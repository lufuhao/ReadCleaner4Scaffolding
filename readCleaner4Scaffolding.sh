#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20180322

Requirements:
    Linux:  perl
    Script: picardrmdup
            bowtie_dedup_mates_separately.sh
            bam_extract_readname_using_region.pl
            bam_filter_by_readname_file.pl
            bam_BamKeepBothMates.pl
            SizeCollectBin_luf.pl
            list2_compar.pl
            bam_exclude_reads_by_windows_size_and_numpairs.pl
            bam_scaffolding_separately.stats_noBioDBSAM.pl
            sspace_seqid_conversion.pl
            sspace_link_scaffolding_ID_to_contig_IDs.pl
    Perl Modules: 
            FuhaoPerl5Lib
    Program: SSPACE, Bowtie, SAMtools, bedtools, bamaddrg, seqtk

Descriptions:

  First Part: Get mapped reads, deduplication and clean
    1. bowtie_dedup_mates_separately.sh for each library
        Mapping each mate (-1 or -2) to reference (-r / -x)
        Merge both mates
        picardrmdup remove duplicates
        Merge libaries
    2. picardrmdup again
    3. SAMtools depth calculate average depth and STDEV
        Find 
        You need to define depth threshold to remove reads in repeat region
      Note: picardrmdup will deduplicate each mate separately,
        so, some reads may have two duplicates

  Second Part: extract reads, remap, and deduplicate
    4. Remove high-depth reads
    5. Extract non-high-depth reads
    6. Bowtie mapping
    7. Clean high depth again
        Keep both mate mapped reads
    8. Extract reads again [1]
        Bowtie map as mate
        rmdup
        Get read names with FLAG 1024
        Exclude these reads from [1]
            * For further analysis
        bowtie map as single again
            * For further analysis
        merge BAMs
            * For further analysis
    9. Extract reads and reference seq for SSPACE
    10. run SSPACE

Options:
  -h    Print this help message
  -1    Fastq R1 files
  -2    Fastq R2 files
  -l    Libraries names
  -x    Bowtie index
  -r    Reference fasta
  -e    Maximum depth threshold;
  -m    Maximum insert size between mates for bowtie -X, default: 300
  -p    Output file prefix, without path
  -d    Delete temporary files
  -t    Number of threads, default: 1
  -sspace path to SSPACE

Example:
  $0 -i ./chr1.fa -t 10

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_q1=''
opt_q2=''
opt_lib=''
opt_ref=''
opt_x=''
opt_t=1
opt_p='MyOutput'
opt_e=0
opt_m=300
opt_d=0
###For debug
startstep=0
endstep=10
path_sspace=''
mean_insertsize=38000
stdev_insertsize='0.2'
strand_orientation='FR'
window_size1=10000
window_size2=20000
num_read_in_window=5
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) opt_q1=$2;shift 2;;
    -2) opt_q2=$2;shift 2;;
    -l) opt_lib=$2;shift 2;;
    -r) opt_ref=$2;shift 2;;
    -x) opt_x=$2;shift 2;;
    -p) opt_p=$2;shift 2;;
    -e) opt_e=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    -d) opt_d=1;shift 1;;
    -a) startstep=$2;shift 2;;
    -b) endstep=$2;shift 2;;
    -sspace) path_sspace=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0
  else
#    echo "I require $1 but it's not installed.  Aborting." >&2
    echo 1
  fi
#  local cmd=$1
#  if command -v $cmd >/dev/null 2>&1;then
#    echo >&2 $cmd "  :  "`command -v $cmd`
#    exit 0
#  else
#    echo >&2 "Error: require $cmd but it's not installed.  Exiting..."
#    exit 1
#  fi
}

###Usage: array=(`split delimiter string`)
SplitStr () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}

#Usage: string=$(join delimiter array)
JoinArr () {
	local separator=$1
	shift 1
	local -a array=(`echo $@`)
	local returnstr=$(printf "$separator%s" "${array[@]}")
	returnstr=${returnstr:1}
	echo $returnstr
}
#abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
CleanTemp () {
	local indtemp;
	
	for indtemp in ${tempfiles[@]}; do
		if [ -e $indtemp ]; then
			echo "SUB(CleanTemp): deleteing $indtemp";
			rm $indtemp > /dev/null 2> /dev/null
		fi
	done
	
	tempfiles=();
	
	return 0;
}

### index BAMs using samtools index
IndexBam () {
	local -a IBbamin=("$@");
	local IBsubinfo='SUB(IndexBam)';
	local IBindbam;
	
	for IBindbam in ${IBbamin[@]}; do
		if [ -z "$IBindbam" ] || [ ! -s "$IBindbam" ]; then
			echo "${IBsubinfo}Error: invalid BAM file: $IBindbam"
			return 1;
		fi
		samtools index $IBindbam > /dev/null 2>&1
		if [ $? -ne 0 ] || [ ! -s $IBindbam.bai ]; then
			echo "${IBsubinfo}Error: samtools index $IBindbam" >&2
			return 1;
		fi
	done
	
	return 0;
}



### BowtieMappingSingle(index, fq, option, outbam)
BowtieMappingSingle () {
	local BMSindex=$1;
	local BMSfq=$2
	local BMSoptions=$3;
	local BMSbamout=$4;
	local BMSsubinfo='BowtieMappingSingle';
	bowtie $BMSoptions $BMSindex $BMSfq  | samtools view -b -h -S -F 4 - > $BMSbamout.notsort.bam
	if [ $? -ne 0 ] || [ ! -s $BMSbamout.notsort.bam ]; then
		echo "${BMSsubinfo}Error: bowtie running failed" >&2
		echo "CMD used: bowtie $BMSoptions $BMSindex $BMSfq  | samtools view -b -h -S -F 4 - > $BMSbamout.notsort.bam" >&2
		return 1;
	fi
	samtools sort -f $BMSbamout.notsort.bam $BMSbamout
	if [ $? -ne 0 ] || [ ! -s $BMSbamout ]; then
		echo "${BMSsubinfo}Error: samtools sort running failed" >&2
		echo "CMD used: samtools sort -f $BMSbamout.notsort.bam $BMSbamout" >&2
		return 1;
	else
		rm $BMSbamout.notsort.bam;
	fi
	
	if IndexBam $BMSbamout; then
		echo "${BMSsubinfo}Info: index $BMSbamout"
	else
		echo "${BMSsubinfo}Error: failed to index $BMSbamout" >&2
		return 1;
	fi
	
	return 0;
}

### BowtieMappingBowtieMappingMates(index, fq1, fq2, option, outbam)
BowtieMappingMates () {
	local BMSindex=$1;
	local BMSfq1=$2
	local BMSfq2=$3
	local BMSoptions=$4;
	local BMSbamout=$5;
	local BMSsubinfo='BowtieMappingMates';
	bowtie $BMSoptions $BMSindex -1 $BMSfq1 -2 $BMSfq2 | samtools view -b -h -S -F 4 - > $BMSbamout.notsort.bam
	if [ $? -ne 0 ] || [ ! -s $BMSbamout.notsort.bam ]; then
		echo "${BMSsubinfo}Error: bowtie running failed" >&2
		echo "CMD used: bowtie $BMSoptions $BMSindex -1 $BMSfq1 -2 $BMSfq2 | samtools view -b -h -S -F 4 - > $BMSbamout.notsort.bam" >&2
		return 1;
	fi
	samtools sort -f $BMSbamout.notsort.bam $BMSbamout
	if [ $? -ne 0 ] || [ ! -s $BMSbamout ]; then
		echo "${BMSsubinfo}Error: samtools sort running failed" >&2
		echo "CMD used: samtools sort -f $BMSbamout.notsort.bam $BMSbamout" >&2
		return 1;
	else
		rm $BMSbamout.notsort.bam;
	fi
	
	if IndexBam $BMSbamout; then
		echo "${BMSsubinfo}Info: index $BMSbamout"
	else
		echo "${BMSsubinfo}Error: failed to index $BMSbamout" >&2
		return 1;
	fi
	
	return 0;
}

### Extract fastq by a file of read ID list
SeqTkSubSeqFastq () {
	local SSFfastqin=$1;
	local SSFidlist=$2;
	local SSFfastqout=$3;
	local SSFsubinfo='SUB(SeqTkSubSeqFastq)';
	
	local SSFnumberlist=0
	local SSFnumberout=0
	
	if [ -z "$SSFfastqin" ] || [ ! -s "$SSFfastqin" ]; then
		echo "${SSFsubinfo}Error: invalid input fastq: $SSFfastqin" >&2
		return 1;
	fi
	if [ -z "$SSFidlist" ] || [ ! -s "$SSFidlist" ]; then
		echo "${SSFsubinfo}Error: invalid fastq ID list: $SSFidlist" >&2
		return 1;
	fi
	if [ -z "$SSFfastqout" ]; then
		echo "${SSFsubinfo}Error: invalid fastq ID list: $SSFfastqout" >&2
		return 1;
	fi
	if [ -e $SSFfastqout ]; then
		rm -rf $SSFfastqout >/dev/null 2>/dev/null
	fi
	
	seqtk subseq $SSFfastqin $SSFidlist > $SSFfastqout
	if [ $? -ne 0 ] || [ ! -s $SSFfastqout ]; then
		echo "${SSFsubinfo}Error: seqtk subseq fastq error" >&2
		echo "${SSFsubinfo}CMD used: seqtk subseq $SSFfastqin $SSFidlist > $SSFfastqout" >&2
		return 1;
	fi
	
	SSFnumberlist=$(perl -ne 'BEGIN{$linenum=0;} $linenum++;END{print $linenum, "\n";}' $SSFidlist)
	SSFnumberout=$(perl -ne 'BEGIN{$linnum=0;} $linenum++;END {if ($linenum%4==0) {print $linenum/4, "\n";} else{print "0\n";}}' $SSFfastqout)
	
	if [ -z "$SSFnumberlist" ] || [ $SSFnumberlist -eq 0 ] || [ $SSFnumberlist -ne $SSFnumberout ]; then
		echo "${SSFsubinfo}Error: seqtk subseq partially failed" >&2
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist" >&2
		echo "${SSFsubinfo}        Total number of R1 IDs extracted: $SSFnumberout" >&2
		return 1;
	else
		echo "${SSFsubinfo}Info: seqtk subseq secceeded"
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist"
	fi
	
	return 0;
}
### Extract fasta by a file of read ID list
SeqTkSubSeqFasta () {
	local SSFfastain=$1;
	local SSFidlist=$2;
	local SSFfastaout=$3;
	local SSFsubinfo='SUB(SeqTkSubSeqFasta)';
	
	local SSFnumberlist=0
	local SSFnumberout=0
	
	if [ -z "$SSFfastain" ] || [ ! -s "$SSFfastain" ]; then
		echo "${SSFsubinfo}Error: invalid input fasta: $SSFfastain" >&2
		return 1;
	fi
	if [ -z "$SSFidlist" ] || [ ! -s "$SSFidlist" ]; then
		echo "${SSFsubinfo}Error: invalid fasta ID list: $SSFidlist" >&2
		return 1;
	fi
	if [ -z "$SSFfastaout" ]; then
		echo "${SSFsubinfo}Error: invalid fasta ID list: $SSFfastaout" >&2
		return 1;
	fi
	if [ -e $SSFfastaout ]; then
		rm -rf $SSFfastaout >/dev/null 2>/dev/null
	fi
	
	seqtk subseq -l 70 $SSFfastain $SSFidlist > $SSFfastaout
	if [ $? -ne 0 ] || [ ! -s $SSFfastaout ]; then
		echo "${SSFsubinfo}Error: seqtk subseq fasta error" >&2
		echo "${SSFsubinfo}CMD used: seqtk subseq -l 70 $SSFfastain $SSFidlist > $SSFfastaout" >&2
		return 1;
	fi
	
	SSFnumberlist=$(wc -l $SSFidlist)
	SSFnumberout=$(perl -ne 'BEGIN{$linnum=0;} $linenum++ if (/^>/); END {print $linenum, "\n";}' < $SSFfastaout)

	if [ -z "$SSFnumberlist" ] || [ $SSFnumberlist -eq 0 ] || [ $SSFnumberlist -ne $SSFnumberout ]; then
		echo "${SSFsubinfo}Error: seqtk lineout partially failed" >&2
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist" >&2
		echo "${SSFsubinfo}        Total number of R1 IDs extracted: $SSFnumberout" >&2
		return 1;
	else
		echo "${SSFsubinfo}Info: seqtk subseq secceeded"
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist"
	fi

	return 0;
}

### merge BAMs using bamadd
BamAddGroup () {
	local BAGbams=$1;
	local BAGout=$2;
	local BAGsubinfo='SUB(BamAddGroup)';
	
	if [ -z "$BAGbams" ]; then
		echo "${BAGsubinfo}Error: invalid input for bamaddrg" >&2
		return 1;
	fi
	
	bamaddrg $BAGbams > $BAGout 2> $BAGout.err
	if [ $? -ne 0 ] || [ ! -s $BAGout ]; then
		echo "${BAGsubinfo}Error: bamaddrg merge error" >&2
		echo "${BAGsubinfo}CMD used: bamaddrg $BAGbams > $BAGout 2> $BAGout.err" >&2
		return 1;
	fi
	rm $BAGout.err > /dev/null 2>&1
	if IndexBam $BAGout; then
		echo "${BAGsubinfo}Info: index $BAGout";
	else
		echo "${BAGsubinfo}Error: failed to index $BAGout" >&2
		return 1;
	fi
	
	return 0;
}



### RunSspace allref ref4sspace config $output
### Global: $opt_t=numthresds, $path_sspace
RunSspace () {
	local RSallref=$1
	local RSref4sspace=$2
	local RSconfig4sspace=$3
	local RSfinalfasta=$4
	
#	default
	local RSmax_allowed_gaps=0
	local RSsspace_out_prefix="$opt_p.sspacerun"
	local RSnumlinks=5
	local RSsspace_running_path=$(pwd)
	local RSsubinfo='SUB(RunSspace)';
	local -a formatted_fasta_arr=()
	local formatted_fasta_file=''
	local RSsspace_final_fasta="$RSsspace_running_path/$RSsspace_out_prefix/$RSsspace_out_prefix.final.scaffolds.fasta"
	local RSsspace_final_evidence="$RSsspace_running_path/$RSsspace_out_prefix/$RSsspace_out_prefix.final.evidence"
	local RSfile_id_convert="$RSsspace_running_path/$RSsspace_out_prefix.id.convert"
	
#	$path_sspace -l $RSconfig4sspace -s $RSref4sspace -g $RSmax_allowed_gaps -T $opt_t -b $RSsspace_out_prefix -k $RSnumlinks -p 1 -v 1
	
	if $path_sspace -l $RSconfig4sspace -s $RSref4sspace -g $RSmax_allowed_gaps -T $opt_t -b $RSsspace_out_prefix -k $RSnumlinks -p 1 -v 1 > sspace.log 2>sspace.err; then
		echo "${RSsubinfo}Info: SSPACE running succeeded"
	else
		echo "${RSsubinfo}Error: SSPACE running failed" >&2
		echo "CMD used: $path_sspace -l $RSconfig4sspace -s $RSref4sspace -g $RSmax_allowed_gaps -T $opt_t -b $RSsspace_out_prefix -k $RSnumlinks -p 1 -v 1" >&2
		return 1;
	fi
	if [ ! -s $RSsspace_final_fasta ]; then
		echo "${RSsubinfo}Error: SSPACE output final fasta error" >&2
		return 1;
	fi
	if [ ! -s $RSsspace_final_evidence ]; then
		echo "${RSsubinfo}Error: SSPACE output final evidence error" >&2
		return 1;
	fi
	formatted_fasta_arr=($(find $RSsspace_running_path/$RSsspace_out_prefix/intermediate_results/ -name "$RSsspace_out_prefix.formattedcontigs*.fasta"))
	if [ ${#formatted_fasta_arr[@]} -eq 1 ]; then
		formatted_fasta_file=${formatted_fasta_arr[0]}
	else
		echo "${RSsubinfo}Error: can not find formated fasta file" >&2
		echo "CMD used: find $RSsspace_running_path/$RSsspace_out_prefix/intermediate_results/ -name $RSsspace_out_prefix.formattedcontigs*.fasta" >&2
		return 1;
	fi

	sspace_seqid_conversion.pl $formatted_fasta_file $RSfile_id_convert
	if [ $? -ne 0 ] || [ ! -s $RSfile_id_convert ]; then
		echo "${RSsubinfo}Error: sspace_seqid_conversion.pl runnng failed" >&2
		echo "CMD used: sspace_seqid_conversion.pl $formatted_fasta_file $RSfile_id_convert" >&2
		return 1;
	fi

	sspace_link_scaffolding_ID_to_contig_IDs.pl $RSsspace_final_evidence $RSfile_id_convert "$RSsspace_out_prefix.joined"
### produce
### 	$RSsspace_out_prefix.joined
### 	$RSsspace_out_prefix.joined.contig
### 	$RSsspace_out_prefix.joined.scaffold
	if [ $? -ne 0 ] || [ ! -s "$RSsspace_running_path/$RSsspace_out_prefix.joined.contig" ] || [ ! -s "$RSsspace_running_path/$RSsspace_out_prefix.joined.scaffold" ]; then
		echo "${RSsubinfo}Error: sspace_link_scaffolding_ID_to_contig_IDs.pl runnng finished" >&2
		echo "CMD used: sspace_link_scaffolding_ID_to_contig_IDs.pl $RSsspace_final_evidence $RSfile_id_convert $RSsspace_out_prefix.joined" >&2
		return 1;
	fi
	local RSassemline1=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.joined.contig")
	local RSassemline2=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.joined.scaffold")
	echo "${RSsubinfo}Info: SSPACE joined $RSassemline1 contigs to $RSassemline2 scaffolds"
	
	if SeqTkSubSeqFasta $RSsspace_final_fasta "$RSsspace_running_path/$RSsspace_out_prefix.joined.scaffold" "$RSsspace_running_path/$RSsspace_out_prefix.joined.scaffold.fasta"; then
		echo "${RSsubinfo}Error: extract joined scaffolds finished"
	else
		echo "${RSsubinfo}Error: extract joined scaffolds failed" >&2
		return 1;
	fi

	if [ ! -s "$RSallref.fai" ]; then
		samtools faidx "$RSallref"
		if [ $? -ne 0 ] || [ ! -s "$RSallref.fai" ]; then
			echo "${RSsubinfo}Error: indexing fasta: $RSallref" >&2
			echo "CMD used: samtools faidx $RSallref" >&2
			return 1;
		fi
	fi

	cut -f 1 "$RSallref.fai" > $RSsspace_running_path/$RSsspace_out_prefix.all.contigs
	if [ $? -ne 0 ] || [ ! -s "$RSsspace_running_path/$RSsspace_out_prefix.all.contigs" ]; then
		echo "${RSsubinfo}Error: extract id names: $RSallref.fai" >&2
		echo "CMD used: cut -f 1 $RSallref.fai > $RSsspace_running_path/$RSsspace_out_prefix.all.contigs" >&2
		return 1;
	fi

	local RSnumline1=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.all.contigs")
	echo "${RSsubinfo}Info: total number of contigs:            $RSnumline1"
	local RSnumline2=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.joined.contig")
	echo "${RSsubinfo}Info: total number of joined contigs:     $RSnumline2"
	local RSnumline3=$(($RSnumline1 - $RSnumline2))
	echo "${RSsubinfo}Info: expect number of un-joined contigs: $RSnumline3"
	echo "${RSsubinfo}Info: Extract unjoined contigs names"
	list2_compar.pl "$RSsspace_running_path/$RSsspace_out_prefix.all.contigs" "$RSsspace_running_path/$RSsspace_out_prefix.joined.contig" "$RSsspace_running_path/$RSsspace_out_prefix.shared.contigs" "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs" "$RSsspace_running_path/$RSsspace_out_prefix.un-expected.contigs"
	if [ $? -ne 0 ] || [ ! -s "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs" ]; then
		echo "${RSsubinfo}Error: extract id names: $RSallref.fai" >&2
		echo "CMD used: cut -f 1 $RSallref.fai > $RSsspace_running_path/$RSsspace_out_prefix.all.contigs" >&2
		return 1;
	fi
	echo "${RSsubinfo}Info: list comparison: $RSsspace_running_path/$RSsspace_out_prefix.all.contigs $RSsspace_running_path/$RSsspace_out_prefix.joined.contig"
	local RSnumline5=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.shared.contigs")
	echo "${RSsubinfo}Info: shared :            $RSnumline5"
	local RSnumline6=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs")
	echo "${RSsubinfo}Info: unique 1:           $RSnumline6"
	local RSnumline7=$(wc -l < "$RSsspace_running_path/$RSsspace_out_prefix.un-expected.contigs")
	echo "${RSsubinfo}Info: unique 2:           $RSnumline7"
	if [ $RSnumline2 -eq $RSnumline5 ] && [ $RSnumline6 -eq $RSnumline3 ] && [ $RSnumline7 -eq 0 ]; then
		echo "${RSsubinfo}Info: total number of un-joined contig names:    $RSnumline6"
	else
		echo "${RSsubinfo}Error: failed to extract un-joined contig names" >&2
		return 1;
	fi
	
	tempfiles+=("$RSsspace_running_path/$RSsspace_out_prefix.shared.contigs")
	tempfiles+=("$RSsspace_running_path/$RSsspace_out_prefix.un-expected.contigs")
	if SeqTkSubSeqFasta $RSallref "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs" "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs.fa"; then
		echo "Info: seqtk subseq fasta success"
	else
		echo "Error: seqtk subseq fasta error" >&2
		echo "CMD used: SeqTkSubSeqFasta $opt_ref $finalrefsid $finalrefsfa" >&2
		return 1;
	fi
	cat "$RSsspace_running_path/$RSsspace_out_prefix.joined.scaffold.fasta" "$RSsspace_running_path/$RSsspace_out_prefix.unjoined.contigs.fa" > $RSfinalfasta
	if [ $? -ne 0 ] || [ ! -s $RSfinalfasta ]; then
		echo "Error: seqtk subseq fasta error" >&2
		echo "CMD used: SeqTkSubSeqFasta $opt_ref $finalrefsid $finalrefsfa" >&2
		return 1
	fi
	tempfiles+=("$RSsspace_running_path/$RSsspace_out_prefix.shared.contigs")
}



### RunStatsPipeline $fastaindex $unpairedstats $RSPpairedstats $RSPmeaninsert $RSPstdevinsert $RSPprefix
RunStatsPipeline () {
	local RSPfastaindex=$1
	local RSPunpairedstats=$2;
	local RSPpairedstats=$3;
	local RSPmeaninsert=$4;
	local RSPstdevinsert=$5;
	local RSPprefix=$6;
	
	local RSPinfo="BASH(RunStatsPipeline)";
	local RSPnumpairs=5
	local RSPrundir=$(pwd)
	
	echo "$RSPprefix:"
	echo "        InDex:    $RSPfastaindex"
	echo "        InSert:   $RSPmeaninsert +- $RSPstdevinsert"
	echo "        Unpaired: $RSPunpairedstats"
	echo "        Paired:   $RSPpairedstats"
	
	
	if [ -d $RSPrundir/$RSPprefix ]; then
		rm -rf $RSPrundir/$RSPprefix > /dev/null 2>&1
	fi
	mkdir -p $RSPrundir/$RSPprefix
	cd $RSPrundir/$RSPprefix
	mkdir -p $RSPrundir/$RSPprefix/1.misassembly
	cd $RSPrundir/$RSPprefix/1.misassembly
	
	cat $RSPunpairedstats $RSPpairedstats > $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats
	if [ $? -ne 0 ] || [ ! -s  $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats ]; then
		echo "${RSPinfo}Error: merge stats error" >&2;
		return 1;
	fi


	echo "${RSPinfo}Info: Detecting good and misassembly errors";
	echo "${RSPinfo}Info: Detecting good and misassembly errors" >&2;
	bam_breaks_by_windowsize_and_numpairs.pl $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats $RSPfastaindex $RSPmeaninsert $RSPstdevinsert $RSPnumpairs $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats.misassembly
	if [ $? -ne 0 ]; then
		echo "${RSPinfo}Error: bam_breaks_by_windowsize_and_numpairs.pl running error" >&2;
		echo "CMD used: bam_breaks_by_windowsize_and_numpairs.pl $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats $RSPfastaindex $RSPmeaninsert $RSPstdevinsert $RSPnumpairs $RSPrundir/$RSPprefix/1.misassembly/$RSPprefix.7.8.merged.stats.misassembly" >&2
		return 1;
	fi


	echo "${RSPinfo}Info: Detecting orientation errors";
	echo "${RSPinfo}Info: Detecting orientation errors" >&2;
	mkdir -p $RSPrundir/$RSPprefix/2.orientation
	cd $RSPrundir/$RSPprefix/2.orientation
	bam_orientations_by_windowsize_and_numpairs.pl $RSPpairedstats $fastaindex $RSPmeaninsert $RSPstdevinsert $RSPnumpairs $RSPrundir/$RSPprefix/2.orientation/$RSPprefix.8.orientation
	if [ $? -ne 0 ]; then
		echo "${RSPinfo}Error: bam_breaks_by_windowsize_and_numpairs.pl running error" >&2;
		echo "CMD used: bam_orientations_by_windowsize_and_numpairs.pl $RSPpairedstats $fastaindex $RSPmeaninsert $RSPstdevinsert $RSPnumpairs $RSPrundir/$RSPprefix/2.orientation/$RSPprefix.8.orientation"
		return 1;
	fi


	echo "${RSPinfo}Info: Detecting orientation errors";
	echo "${RSPinfo}Info: Detecting orientation errors" >&2;
	mkdir -p $RSPrundir/$RSPprefix/3.indel
	cd $RSPrundir/$RSPprefix/3.indel
	RSPmininsert=$(($RSPmeaninsert-$RSPstdevinsert-$RSPstdevinsert))
	RSPmaxinsert=$(($RSPmeaninsert+$RSPstdevinsert+$RSPstdevinsert))
	echo "${RSPinfo}Info: Min insert: $RSPmininsert"
	echo "${RSPinfo}Info: Max insert: $RSPmaxinsert"
	bam_improper_insert_size_by_windowsize_and_numpairs.pl $RSPpairedstats $RSPmininsert $RSPmaxinsert $RSPnumpairs $RSPrundir/$RSPprefix/3.indel/$RSPprefix.8.improper
	if [ $? -ne 0 ]; then
		echo "${RSPinfo}Error: bam_improper_insert_size_by_windowsize_and_numpairs.pl running error" >&2;
		echo "CMD used: bam_improper_insert_size_by_windowsize_and_numpairs.pl $RSPpairedstats $RSPmininsert $RSPmaxinsert $RSPnumpairs $RSPrundir/$RSPprefix/3.indel/$RSPprefix.8.improper"
		return 1;
	fi
	
	
	return 0;
}


#################### Command test ###################################
if [ $(CmdExists 'bowtie') -ne 0 ]; then
	echo "Error: CMD 'bowtie' in PROGRAM 'bowtie' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bowtie-build') -ne 0 ]; then
	echo "Error: CMD 'bowtie-build' in PROGRAM 'bowtie' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'picardrmdup') -ne 0 ]; then
	echo "Error: script 'picardrmdup' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'samtools') -ne 0 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bowtie_dedup_mates_separately.sh') -ne 0 ]; then
	echo "Error: script 'bowtie_dedup_mates_separately.sh' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bam_extract_readname_using_region.pl') -ne 0 ]; then
	echo "Error: script 'bam_extract_readname_using_region.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bam_filter_by_readname_file.pl') -ne 0 ]; then
	echo "Error: script 'bam_filter_by_readname_file.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'seqtk') -ne 0 ]; then
	echo "Error: CMD 'seqtk' in PROGRAM 'seqtk' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'SizeCollectBin_luf.pl') -ne 0 ]; then
	echo "Error: script 'SizeCollectBin_luf.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bamaddrg') -ne 0 ]; then
	echo "Error: CMD 'bamaddrg' in PROGRAM 'bamaddrg' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'list2_compar.pl') -ne 0 ]; then
	echo "Error: script 'list2_compar.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bam_scaffolding_separately.stats_noBioDBSAM.pl') -ne 0 ]; then
	echo "Error: script 'bam_scaffolding_separately.stats_noBioDBSAM.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bam_exclude_reads_by_windows_size_and_numpairs.pl') -ne 0 ]; then
	echo "Error: script 'bam_scaffolding_separately.stats_noBioDBSAM.pl' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ -z "$path_sspace" ] || [ ! -s $path_sspace ]; then
	echo "Error: SSPACE if needed but not found" >&2
	exit 127;
fi
if [ $(CmdExists 'sspace_seqid_conversion.pl') -ne 0 ]; then
	echo "Error: script 'sspace_seqid_conversion.pl' is required but not found.  Aborting..." >&2
	exit 127;
fi
if [ $(CmdExists 'sspace_link_scaffolding_ID_to_contig_IDs.pl') -ne 0 ]; then
	echo "Error: script 'sspace_link_scaffolding_ID_to_contig_IDs.pl' is required but not found.  Aborting..." >&2
	exit 127;
fi


#################### Defaults #######################################
RunDir=$RunPath;
threads=$opt_t;
numgap=0;
declare -a fastq1arr=();
declare -a fastq2arr=();
declare -a libraryarr=();
declare -a sub1seq1arr=();
declare -a sub1seq2arr=();
declare -a sub2seq1arr=();
declare -a sub2seq2arr=();
declare -a bam1filesR1arr=();
declare -a bam1filesR2arr=();
declare -a tempfiles=()
declare -a bamaddrgfiles=()

#################### Input and Output ###############################fasta
if [ -z "$opt_q1" ]; then
	echo "Error: invalid fastq R1 [option -1]" >&2
	exit 10;
fi
if [ -z "$opt_q2" ]; then
	echo "Error: invalid fastq R2 [option -2]" >&2
	exit 10;
fi
if [ -z "$opt_lib" ]; then
	echo "Error: invalid library names [option -l]" >&2
	exit 10;
fi
fastq1arr=(${opt_q1//,/ })
fastq2arr=(${opt_q2//,/ })
libraryarr=(${opt_lib//,/ })
#fastq1arr=$(echo $opt_q1 | tr ',' "\n")
#fastq2arr=$(echo $opt_q2 | tr ',' "\n")

if [ ${#fastq1arr[@]} -ne ${#fastq2arr[@]} ] || [ ${#fastq1arr[@]} -ne ${#libraryarr[@]} ] || [ ${#fastq1arr[@]} -lt 1 ]; then
	echo "Error: uneven input fastq/libraryname files" >&2
	exit 10;
fi
echo "##### Libraries SUMMARY #####"
echo "$opt_lib"
echo "${libraryarr[@]}"
for (( i=0; i < ${#libraryarr[@]}; i++ )); do
	if [ -z "${libraryarr[$i]}" ]; then
		echo "Error: invalid individual library name: ${libraryarr[$i]}" >&2
		exit 10;
	fi
	if [ ! -s ${fastq1arr[$i]} ]; then
		echo "Error: invalid individual R1 FastQ: ${fastq1arr[$i]}" >&2
		exit 10;
	fi
	if [ ! -s ${fastq2arr[$i]} ]; then
		echo "Error: invalid individual R2 FastQ: ${fastq2arr[$i]}" >&2
		exit 10;
	fi
	fastq1arr[$i]=$(readlink -f ${fastq1arr[$i]})
	fastq2arr[$i]=$(readlink -f ${fastq2arr[$i]})
	echo "    Library $i ${libraryarr[$i]}"
	echo "        FastQ : R1: ${fastq1arr[$i]}"
	echo "                R2: ${fastq2arr[$i]}"
done
if [ -z "$opt_e" ] || [[ ! "$opt_e" =~ ^[0-9]+$ ]] || [ !  $opt_e -ge 0 ]; then
	echo "Error: depth threshold [-e INT] option error" >&2
	exit 10;
fi
if [ -z "$opt_t" ] || [[ ! "$opt_t" =~ ^[0-9]+$ ]] || [ !  $opt_t -gt 0 ]; then
	echo "Error: threads [-t INT] option error" >&2
	exit 10;
fi



#################### Main ###########################################



### Building index
step=0
cd $RunDir/
if [ -z "$opt_x" ]; then
	echo -e "\n\n#####    Step$step: building index    #####"
	echo -e "\n\n#####    Step$step: building index    #####" >&2
	if [ ! -s $opt_ref ]; then
		echo "Error: invalid reference sequence file" >&2
		exit 10;
	fi
	if [ ! -d $RunDir/$step.index ]; then
		mkdir -p $RunDir/$step.index
	fi
	cd $RunDir/$step.index
	rm $RunDir/$step.index.* > /dev/null 2>&1
	echo "Info: building index using $opt_ref"
	bowtie-build -f $opt_ref BWTindex > bowtie-build.log 2 > bowtie-build.err
	if [ $? -ne 0 ]; then
		echo "Error: bowtie-build index error" >&2
		echo "CMD used: bowtie-build -f $opt_ref BWTindex" >&2
		exit 10;
	fi
	rm bowtie-build.log bowtie-build.err > /dev/null 2>&1
	opt_x=$RunDir/0.index/BWTindex
	echo "bowtie-build -f $opt_ref BWTindex" > bowtie-build.success
else
	opt_x=$(readlink -f $opt_x)
fi



bamfirstmerge="$RunDir/$opt_p.1.merge.bam"
bamfirstrmdup="$RunDir/$opt_p.1.merge.rmdup.bam"
bam2output="$RunDir/$opt_p.2.merge.rmdup.rmdup.bam"
bam_unpairedstats="$RunDir/$opt_p.7.unpaired.stats"
bam_pairedstats="$RunDir/$opt_p.8.paired.stats"
finalcleanR1=$RunDir/$opt_p.8.paired.R1.fastq.gz
finalcleanR2=$RunDir/$opt_p.8.paired.R2.fastq.gz
finalrefsid=$RunDir/$opt_p.9.finalrefsid
finalrefsfa=$RunDir/$opt_p.10.finalrefsid.fasta
finalreadid=$RunDir/$opt_p.9.finalreadid
finalreadfqR1=$RunDir/$opt_p.10.finalreadid.R1.fastq
finalreadfqR2=$RunDir/$opt_p.10.finalreadid.R2.fastq
finalmerge=$RunDir/$opt_p.10.joinedscaffolds.unjoinedcontigs.fasta
### First Part
cd $RunDir/
echo -e "\n\n\n"
echo -e "\n\n\n" >&2
if [ $opt_e -eq 0 ]; then
	echo "### First Part"
	echo "### First Part" >&2


	### Bowtie
	step=1
	echo -e "\n\n#####    Step$step: Running Bowtie    #####"
	echo -e "\n\n#####    Step$step: Running Bowtie    #####" >&2
	cd $RunDir/
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step1"
		exit 0
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.bowtie ]; then
			rm -rf $RunDir/$step.bowtie > /dev/null 2>&1
		fi
		mkdir $RunDir/$step.bowtie
	fi
	if [ ! -d $RunDir/$step.bowtie ]; then
		echo "Error: subfolder not found: $RunDir/$step.bowtie" >&2
		exit 10;
	fi
	cd $RunDir/$step.bowtie
	bamfiles2num=0
	bamfiles2mergestr=' '
	bamfiles3mergestr=' '
	for (( i=0; i < ${#libraryarr[@]}; i++ )); do
		outputprefix="${libraryarr[$i]}.firstBT"
		echo "    Library $i ${libraryarr[$i]}"
		echo "        FastQ : R1: ${fastq1arr[$i]}"
		echo "                R2: ${fastq2arr[$i]}"
		echo "    Library $i ${libraryarr[$i]}" >&2
		echo "        FastQ : R1: ${fastq1arr[$i]}" >&2
		echo "                R2: ${fastq2arr[$i]}" >&2
		if [ $startstep -le $step ]; then
			echo "Info: bowtie mapping and dedup separately"
			bowtie_dedup_mates_separately.sh -1 ${fastq1arr[$i]} -2 ${fastq2arr[$i]} -b " -q -p $threads -v $numgap -a --sam $opt_x " -p $RunDir/$step.bowtie/$outputprefix > $RunDir/$step.bowtie/${libraryarr[$i]}.bowtie_dedup_mates_separately.log 2> $RunDir/$step.bowtie/${libraryarr[$i]}.bowtie_dedup_mates_separately.err
			if [ $? -ne 0 ] || [ ! -s $RunDir/$step.bowtie/$outputprefix.1.2.st.merge.rmdup.bam ]; then
				echo "Error: bowtie running error" >&2
				exit 10;
			fi
			rm bowtie_dedup_mates_separately.log bowtie_dedup_mates_separately.err > /dev/null 2>&1
		fi
		
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.st.bam")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.st.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.2.st.bam")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.2.st.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.bowtie_dedup_mates_separately.log")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.bowtie_dedup_mates_separately.err")
		
		bamfiles2mergestr="$bamfiles2mergestr -b $RunDir/$step.bowtie/$outputprefix.1.2.st.merge.bam -s LIB -r ${libraryarr[$i]}"
		bamfiles3mergestr="$bamfiles3mergestr -b $RunDir/$step.bowtie/$outputprefix.1.2.st.merge.rmdup.bam -s LIB -r ${libraryarr[$i]}"
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.2.st.merge.bam")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.2.st.merge.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.2.st.merge.rmdup.bam")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.2.st.merge.rmdup.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/$outputprefix.1.2.st.merge.rmdup.metrix")
		let bamfiles2num++
	done
	if [ $bamfiles2num -lt 1 ]; then
		echo "Error: bowtie output No BAMs" >&2
		exit 10;
	fi
	if [ $startstep -le $step ]; then
		echo "Info: bamaddrg merge origin BAMs"
		if BamAddGroup "$bamfiles2mergestr" "$bamfirstmerge"; then
			echo "Info: bamaddrg merge success" >&2
		else
			echo "Error: bamaddrg merge error" >&2
			exit 10;
		fi
		echo "Info: bamaddrg merge separately deduplicated BAMs"
		if BamAddGroup "$bamfiles3mergestr" "$bamfirstrmdup"; then
			echo "Info: bamaddrg merge rmdup success" >&2
		else
			echo "Error: bamaddrg merge rmdup error" >&2
			exit 10;
		fi
	fi
	

	### picard rmdup
	step=2
	echo -e "\n\n#####    Step$step: picardrmdup    #####"
	echo -e "\n\n#####    Step$step: picardrmdup    #####" >&2
	cd $RunDir/
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.rmdup ]; then
			rm -rf $RunDir/$step.rmdup > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.rmdup
	fi
	if [ ! -d $RunDir/$step.rmdup ]; then
		echo "Error: subfolder not found: $RunDir/$step.rmdup" >&2
		exit 10;
	fi
	cd $RunDir/$step.rmdup
	bam2out=${bamfirstrmdup##*/}
	bam2out=${bam2out%.bam}
	bam2out="$RunDir/$step.rmdup/$bam2out.rmdup.bam"
	if [ ! -s $bamfirstrmdup ]; then
		echo "Error: last step output not found for Step$step: $bamfirstrmdup" >&2
		exit 10;
	fi
	if [ $startstep -le $step ]; then
		echo "Info: picard MarkDuplicates to dedup"
		picardrmdup -d $bamfirstrmdup > picardrmdup.log 2> picardrmdup.err
		if [ $? -ne 0 ] || [ ! -s "$bam2out" ] || [ ! -s "$bam2out.bai" ]; then
			echo "Error: picardrmdup running error" >&2
			echo "CMD used: picardrmdup -d $bamfirstrmdup" >&2
			exit 10;
		fi
		mv $bam2out $bam2output
		mv "$bam2out.bai" "$bam2output.bai"
		rm picardrmdup.log picardrmdup.err > /dev/null 2>&1
	fi



	### depth
	step=3
	echo -e "\n\n#####    Step$step: depth    #####"
	echo -e "\n\n#####    Step$step: depth    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.depth ]; then
			rm -rf $RunDir/$step.depth > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.depth
	fi
	if [ ! -d $RunDir/$step.depth ]; then
		echo "Error: subfolder not found: $RunDir/$step.depth" >&2
		exit 10;
	fi
	cd $RunDir/$step.depth
#input
	if [ $startstep -le $step ]; then
		if [ ! -s "$bam2output" ]; then
			echo "Error: last step output not found at step$step: $bam2output" >$2
			exit 10;
		fi
	fi
#output

#run
	if [ $startstep -le $step ]; then
		echo "Info: Generating depth info using samtools depth"
		samtools depth $bam2output > $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth
		if [ $? -ne 0 ] || [ ! -s $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth ]; then
			echo "Error: samtools depth" >&2
			echo "CMD used: samtools depth $bam2output > $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth" >&2
			exit 10;
		fi
		echo "Info: calculate average depth and STDEV"
		perl -lane 'print $F[2];' $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth > $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3
		BAMAVEDEPTH=$(perl -ne 'chomp;$sum+=$_;$n++;END{print $sum/$n, "\n";}' $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3)
		export BAMAVEDEPTH
#		echo "Average depth=$BAMAVEDEPTH"
		perl -ne 'BEGIN{$n=0;$sum=0;$max=0;} chomp;$sum= $sum + ($_-$ENV{"BAMAVEDEPTH"})*($_-$ENV{"BAMAVEDEPTH"});if ($_>$max) {$max=$_;} $n++;END{print "\tMax: $max\n\tAVERAGE depth:    ", $ENV{"BAMAVEDEPTH"}, "\n\tTotal:            $n\n\tSum:              $sum\n\tsum/n:            ", $sum/$n, "\n\tSTDEV:            ", sqrt($sum/$n), "\n";}' $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3
		echo "Info: summarizing depth"
		SizeCollectBin_luf.pl $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3 1 > $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3.bin1
		if [ $? -ne 0 ] || [ ! -s $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3.bin1 ]; then
			echo "Error: SizeCollectBin_luf.pl error" >&2
			exit 1;
		fi
	fi
	echo -e "\n\n\n"
	echo "Info: now plot $RunDir/$step.depth/$opt_p.merge.rmdup.bam.depth.col3.bin1"
	echo "Info: decide allowed minimum depth for option -e [INT]"
	echo -e "\n\n\n"
	
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi

elif [ $opt_e -gt 0 ]; then
	step4bamoutput=' '
	step4depthbed=' '
	step6bam1out=' '
	step6bam2out=' '
	step6bammerge=' '
	step7bamout=' '
	step8bamout1=' '
	step8bamout2=' '
	step8bamfinal=' '

#################### remove repeat ###########################################
	step=4
	echo -e "\n\n#####    Step$step: remove high depth reads    #####"
	echo -e "\n\n#####    Step$step: remove high depth reads    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.RMhighdepth ]; then
			rm -rf $RunDir/$step.RMhighdepth > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.RMhighdepth
	fi
	if [ ! -d $RunDir/$step.RMhighdepth ]; then
		echo "Error: subfolder not found: $RunDir/$step.RMhighdepth" >&2
		exit 10;
	fi
	cd $RunDir/$step.RMhighdepth
#input
	if [ $startstep -le $step ]; then
		if [ ! -s $bam2output ]; then
			echo "Error: can not find $bam2output" >&2
			exit 10;
		fi
		if [ ! -s $bamfirstmerge ]; then
			echo "Error: can not find $bamfirstmerge" >&2
			exit 10;
		fi
	fi
#output
	step4bamoutput="$RunDir/$step.RMhighdepth/$opt_p.4.merge.rmdup.rmdup.rmrepeat.bam"
	step4depthbed="$RunDir/$step.RMhighdepth/$opt_p.4.merge.rmdup.bam.genomecov"
#run
	if [ $startstep -le $step ]; then
		DEPTH_THRESHOLD=$opt_e
		export DEPTH_THRESHOLD
		echo "Info: DEPTH_THRESHOLD was set to $DEPTH_THRESHOLD"
		echo "Info: Get high-depth region in BED format"
		bedtools genomecov -bg -ibam $bam2output | perl -ne 'chomp; @arr=split(/\t/);if ($arr[3]>$ENV{"DEPTH_THRESHOLD"}) {print $_, "\n";} else {next;}' | bedtools merge -i - > $step4depthbed
		if [ $? -ne 0 ] || [ ! -s $step4depthbed ]; then
			echo "Error: bedtools genomecov" >&2
			exit 10;
		fi
		echo "Info: extract read names in high-depth region"
		bam_extract_readname_using_region.pl $bamfirstmerge $step4depthbed $step4depthbed.exclude.readname
		if [ $? -ne 0 ] || [ ! -s $step4depthbed.exclude.readname ]; then
			echo "Error: bam_extract_readname_using_region.pl" >&2
			exit 10;
		fi
		echo "Info: exclude extracted read names in high-depth region"
		bam_filter_by_readname_file.pl $bam2output $step4depthbed.exclude.readname 0 $step4bamoutput
		if [ $? -ne 0 ] || [ ! -s "$step4bamoutput" ]; then
			echo "Error: bam_filter_by_readname_file.pl" >&2
			exit 10;
		fi
		if IndexBam "$step4bamoutput"; then
			echo "Info: index BAMs" >&2
		else
			echo "Error: failed to index BAMs" >&2
			exit 10;
		fi
	fi




##################### Extract reads ####################################
	step=5
	echo -e "\n\n#####    Step$step: extract read IDs    #####"
	echo -e "\n\n#####    Step$step: extract read IDs    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0;
	fi
	if [ $startstep -lt $step ]; then
		if [ -d $RunDir/$step.extractReadIDs ]; then
			rm -rf $RunDir/$step.extractReadIDs > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.extractReadIDs
	fi
	if [ ! -d $RunDir/$step.extractReadIDs ]; then
		echo "Error: subfolder not found: $RunDir/$step.extractReadIDs" >&2
		exit 10;
	fi
	cd $RunDir/$step.extractReadIDs
#Default
	declare -a fastqidlist=();
#input
	if [ $startstep -lt $step ]; then
		if [ ! -s $step4bamoutput ]; then
			echo "Error: last step output not found at step$step: $step4bamoutput" >&2
			exit 10;
		fi
	fi
#output

#run
	for (( i=0; i < ${#libraryarr[@]}; i++ )); do
		echo "    Library $i ${libraryarr[$i]}"
		echo "        FastQ : R1: ${fastq1arr[$i]}"
		echo "                R2: ${fastq2arr[$i]}"
		echo "    Library $i ${libraryarr[$i]}" >&2
		echo "        FastQ : R1: ${fastq1arr[$i]}" >&2
		echo "                R2: ${fastq2arr[$i]}" >&2
		if [ $startstep -lt $step ]; then
			echo "Info: samtools view to extract read names"
			samtools view -r "${libraryarr[$i]}" $step4bamoutput | perl -lane 'BEGIN{%idhash=()};unless (exists $idhash{$F[0]}) {print $F[0]; $idhash{$F[0]}++;}' > $RunDir/$step.extractReadIDs/${libraryarr[$i]}.id
			if [ $? -ne 0 ] || [ ! -s $RunDir/$step.extractReadIDs/${libraryarr[$i]}.id ]; then
				echo "Error: extract read IDs $RunDir/$step.extractReadIDs/${libraryarr[$i]}.id" >&2
				exit 10;
			fi
		fi
		if [ $startstep -eq $step ]; then
			echo "Info: check extracted read names"
			if [ ! -s $RunDir/$step.extractReadIDs/${libraryarr[$i]}.id ]; then
				echo "Error: extracted read ID not exists: $RunDir/$step.extractReadIDs/${libraryarr[$i]}.id" >&2
				exit 10
			fi
		fi
		if [ $startstep -le $step ]; then
			echo "Info: Extract read R1"
			if SeqTkSubSeqFastq "${fastq1arr[$i]}" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R1.subseq.fastq"; then
				echo "Info: seqtk subseq R1 fastq success"
			else
				echo "Error: seqtk subseq R1 fastq error" >&2
				exit 10;
			fi
			echo "Info: Extract read R2"
			if SeqTkSubSeqFastq "${fastq2arr[$i]}" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R2.subseq.fastq"; then
				echo "Info: seqtk subseq R2 fastq success"
			else
				echo "Error: seqtk subseq R2 fastq error" >&2
				exit 10;
			fi
		fi
		sub1seq1arr=("${sub1seq1arr[@]}" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R1.subseq.fastq");
		sub1seq2arr=("${sub1seq2arr[@]}" "$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R2.subseq.fastq");
#		tempfiles+=("$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id")
		fastqidlist+=("$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id")
		tempfiles+=("$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R1.subseq.fastq")
		tempfiles+=("$RunDir/$step.extractReadIDs/${libraryarr[$i]}.id.R2.subseq.fastq")
	done
	
	
	for indidfile in ${fastqidlist[@]}; do
		if [ -e "$indidfile" ]; then
			wc -l $indidfile
		fi
	done
	fastqidlist=();
	if [ ${#sub1seq1arr[@]} -ne ${#libraryarr[@]} ] || [ ${#sub1seq2arr[@]} -ne ${#libraryarr[@]} ]; then
		echo "Error: some seqtk subseq might be failed" >&2
		exit 10;
	fi



	### bowtie mapping again
	step=6
	echo -e "\n\n#####    Step$step: bowtie    #####"
	echo -e "\n\n#####    Step$step: bowtie    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0;
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.bowtie ]; then
			rm -rf $RunDir/$step.bowtie > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.bowtie
	fi
	if [ ! -d $RunDir/$step.bowtie ]; then
		echo "Error: subfolder not found: $RunDir/$step.bowtie" >&2
		exit 10;
	fi
	cd $RunDir/$step.bowtie
	for (( i=0; i < ${#libraryarr[@]}; i++ )); do
		if [ $startstep -le $step ]; then
			echo "Library $i ${libraryarr[$i]}"
			echo "Library $i ${libraryarr[$i]}" >&2
			if [ -s "${sub1seq1arr[$i]}" ]; then
				echo "    FastQ : R1: ${sub1seq1arr[$i]}"
				echo "    FastQ : R1: ${sub1seq1arr[$i]}" >&2
			else
				echo "    FastQ : R1: ${sub1seq1arr[$i]} NOT FOUND" >&2
			fi
			if [ -s "${sub1seq1arr[$i]}" ]; then
				echo "            R2: ${sub1seq2arr[$i]}"
				echo "            R2: ${sub1seq2arr[$i]}" >&2
			else
				echo "            R2: ${sub1seq2arr[$i]} NOT FOUND" >&2
			fi
			if [ -s "${sub1seq1arr[$i]}" ] && [ -s "${sub1seq1arr[$i]}" ]; then
				echo "Info: bowtie mapping R1"
				if BowtieMappingSingle "$opt_x" "${sub1seq1arr[$i]}" " -q -p $threads -v $numgap -a --sam" "$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam" > $RunDir/$step.bowtie/${libraryarr[$i]}.R1.bowtie.log 2> $RunDir/$step.bowtie/${libraryarr[$i]}.R1.bowtie.err; then
					echo "Info: bowtie mapping R1 succeeds"
				else
					echo "Error: bowtie mapping R1 error: $RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam" >&2
					exit 10;
				fi
				echo "Info: bowtie mapping R2"
				if BowtieMappingSingle "$opt_x" "${sub1seq2arr[$i]}" " -q -p $threads -v $numgap -a --sam" "$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam"  > $RunDir/$step.bowtie/${libraryarr[$i]}.R2.bowtie.log 2> $RunDir/$step.bowtie/${libraryarr[$i]}.R2.bowtie.err; then
					echo "Info: bowtie mapping R2 succeeds"
				else
					echo "Error: bowtie mapping R2 error: $RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam" >&2
					exit 10;
				fi
			fi
		fi
#		if [ -s "${sub1seq1arr[$i]}" ] && [ -s "${sub1seq1arr[$i]}" ]; then
			bam1filesR1arr+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam")
			bam1filesR2arr+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam")
			bam1filesR1str="$bam1filesR1str -b $RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam -s LIB -r ${libraryarr[$i]}"
			bam1filesR2str="$bam1filesR2str -b $RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam -s LIB -r ${libraryarr[$i]}"
#		fi
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R1.st.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.subseq.R2.st.bam.bai")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.R1.bowtie.log")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.R1.bowtie.err")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.R2.bowtie.log")
		tempfiles+=("$RunDir/$step.bowtie/${libraryarr[$i]}.R2.bowtie.err")
	done
	if [ ${#bam1filesR1arr[@]} -ne ${#libraryarr[@]} ] || [ ${#bam1filesR2arr[@]} -ne ${#libraryarr[@]} ]; then
		echo "Error: some bowtie mapping might be failed" ${#libraryarr[@]} ${#bam1filesR1arr[@]} ${#bam1filesR2arr[@]} >&2
		exit 10;
	else
		bam1filesR1arr=()
		bam1filesR2arr=()
	fi
	
	step6bam1out="$RunDir/$step.bowtie/$opt_p.subseq.R1.st.bam"; 
	step6bam2out="$RunDir/$step.bowtie/$opt_p.subseq.R2.st.bam";
	step6bammerge="$RunDir/$step.bowtie/$opt_p.subseq.R1R2.st.bam"
	if [ $startstep -le $step ]; then
		echo "Info: bamaddrg merge R1"
		if BamAddGroup "$bam1filesR1str" "$step6bam1out"; then
			echo "Info: bamaddrg merge R1 success"
		else
			echo "Error: bamaddrg merge R1 error" >&2
			exit 10;
		fi
		echo "Info: bamaddrg merge R2"
		if BamAddGroup "$bam1filesR2str" "$step6bam2out"; then
			echo "Info: bamaddrg merge R2 success"
		else
			echo "Error: bamaddrg merge R2 error" >&2
			exit 10;
		fi
		echo "Info: samtools merge"
		samtools merge "$step6bammerge" "$step6bam1out" "$step6bam2out"
		if [ $? -ne 0 ] || [ ! -s "$step6bammerge" ]; then
			echo "Error: samtools merge error" >&2
			echo "CMD used: samtools merge $step6bammerge $step6bam1out $step6bam2out"
			exit 10;
		fi
		if IndexBam "$step6bammerge"; then
			echo "Info: index BAMs"
		else
			echo "Error: failed to index BAMs" >&2
			exit 10;
		fi
	fi
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi
	
	
	
	###### RMrepeat
	step=7
	echo -e "\n\n#####    Step$step: RMrepeat    #####"
	echo -e "\n\n#####    Step$step: RMrepeat    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step$step"
		exit 0
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.RMrepeat ]; then
			rm -rf $RunDir/$step.RMrepeat > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.RMrepeat
	fi
	if [ ! -d $RunDir/$step.RMrepeat ]; then
		echo "Error: subfolder not found: $RunDir/$step.RMrepeat" >&2
		exit 10;
	fi
	cd $RunDir/$step.RMrepeat
#input
	if [ ! -s $step4depthbed ]; then
		echo "Error: last step output not found at step$step: $step4depthbed" >&2
		exit 10;
	fi
	if [ ! -s "$step6bammerge" ]; then
		echo "Error: last step output not found at step$step: $step4depthbed" >&2
		exit 10;
	fi
	if [ ! -s "$step6bam1out" ]; then
		echo "Error: last step output not found at step$step: $step6bam1out" >&2
		exit 10;
	fi
	if [ ! -s "$step6bam2out" ]; then
		echo "Error: last step output not found at step$step: $step6bam2out" >&2
		exit 10;
	fi
#output

#run
	if [ $startstep -le $step ]; then
		if [ ! -s "$step6bam1out" ]; then
			echo "Error: not found: input R1.bam at step$step" >&2
			exit 10;
		fi
		if [ ! -s "$step6bam2out" ]; then
			echo "Error: not found: input R2.bam at step$step" >&2
			exit 10;
		fi
		
		echo "Info: re-extract R1R2 readnames (for double check)"
		bam_extract_readname_using_region.pl $step6bammerge $step4depthbed $RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname
		if [ $? -ne 0 ] || [ ! -e $RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname ]; then
			echo "Error: bam_extract_readname_using_region.pl R2" >&2
			exit 10;
		fi
		num2exlude=$(perl -ne 'BEGIN{$linenum=0;} $linenum++;END{print $linenum, "\n";}' "$RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname")
		echo "Info: total number of read pairss to exclude: $num2exlude"
		if [ $num2exlude -gt 0 ]; then
			echo "Info: filter out R1 readnames"
			bam_filter_by_readname_file.pl $step6bam1out $RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname 0 $RunDir/$step.RMrepeat/subseq.R1.st.rmrepeat.bam
			if [ $? -ne 0 ] || [ ! -s $RunDir/$step.RMrepeat/subseq.R1.st.rmrepeat.bam ]; then
				echo "Error: bam_filter_by_readname_file.pl R1" >&2
				exit 10;
			fi
			if IndexBam "$RunDir/$step.RMrepeat/subseq.R1.st.rmrepeat.bam"; then
				echo "Info: index BAMs"
			else
				echo "Error: failed to index BAMs" >&2
				exit 10;
			fi
			
			echo "Info: filter out R2 readnames"
			bam_filter_by_readname_file.pl $step6bam2out $RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname 0 $RunDir/$step.RMrepeat/subseq.R2.st.rmrepeat.bam
			if [ $? -ne 0 ] || [ ! -s $RunDir/$step.RMrepeat/subseq.R2.st.rmrepeat.bam ]; then
				echo "Error: bam_filter_by_readname_file.pl R2" >&2
				exit 10;
			fi
			if [ $? -ne 0 ] || [ ! -s $RunDir/$step.RMrepeat/subseq.R2.st.rmrepeat.bam.bai ]; then
				echo "Error: samtools index R2 error" >&2
				exit 10;
			fi
			if IndexBam "$RunDir/$step.RMrepeat/subseq.R2.st.rmrepeat.bam"; then
				echo "Info: index BAMs"
			else
				echo "Error: failed to index BAMs" >&2
				exit 10;
			fi
			
			step6bam1out=$RunDir/$step.RMrepeat/subseq.R1.st.rmrepeat.bam
			step6bam2out=$RunDir/$step.RMrepeat/subseq.R2.st.rmrepeat.bam
		else
			echo "Info: No read name to exclude"
			if [ -e "$RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname" ]; then
				rm -rf "$RunDir/$step.RMrepeat/subseq.R1R2.st.bam.genomecov.exclude.readname" >/dev/null 2>&1
			fi
		fi
	fi

	if [ $startstep -le $step ]; then
		echo "Info: keep only both mates mapped"
		bam_BamKeepBothMates.pl "$step6bam1out" "$step6bam2out" "$step6bam1out.Paired.bam" "$step6bam2out.Paired.bam" "$step6bam1out.UnPaired.bam" "$step6bam2out.UnPaired.bam";
		if [ $? -ne 0 ] || [ ! -s "$step6bam1out.Paired.bam" ] || [ ! -s "$step6bam2out.Paired.bam" ]; then
			echo "Error: bam_BamKeepBothMates.pl error" >&2
			exit 10;
		fi
		if IndexBam "$step6bam1out.Paired.bam"; then
			echo "Info: index BAMs"
		else
			echo "Error: failed to index BAMs" >&2
			exit 10;
		fi
		if IndexBam "$step6bam2out.Paired.bam"; then
			echo "Info: index BAMs"
		else
			echo "Error: failed to index BAMs" >&2
			exit 10;
		fi
		bam_scaffolding_separately.stats_noBioDBSAM.pl "$step6bam1out.UnPaired.bam" "$step6bam2out.UnPaired.bam" $bam_unpairedstats
		if [ $? -ne 0 ] || [ ! -s "$bam_unpairedstats" ]; then
			echo "Warnings: failed to stats unpaired BAMs" >&2
		fi
	fi
#	tempfiles+=("$step6bam1out")
#	tempfiles+=("$step6bam1out.bai")
#	tempfiles+=("$step6bam2out")
#	tempfiles+=("$step6bam2out.bai")
	step7bamout="$RunDir/$step.RMrepeat/$opt_p.7.subseq.R1R2.st.rmrepeat.paired.bam"
	if [ $startstep -le $step ]; then
		echo "Info: samtools merge paired"
		samtools merge "$step7bamout" "$step6bam1out.Paired.bam" "$step6bam2out.Paired.bam" >/dev/null 2>&1
		if [ $? -ne 0 ] || [ ! -s "$step7bamout" ]; then
			echo "Error: samtools merge error" >&2
			echo "CMD used: samtools merge $step7bamout $step6bam1out.Paired.bam $step6bam2out.Paired.bam" >&2
			exit 10;
		fi
		if IndexBam "$step7bamout"; then
			echo "Info: index BAMs"
		else
			echo "Error: failed to index BAMs" >&2
			exit 10;
		fi
	fi
	tempfiles+=("$step6bam1out.Paired.bam")
	tempfiles+=("$step6bam1out.Paired.bam.bai")
	tempfiles+=("$step6bam2out.Paired.bam")
	tempfiles+=("$step6bam2out.Paired.bam.bai")
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi



### exclude few duplicated reads
	step=8
	echo -e "\n\n#####    Step$step: Exclude duplicates    #####"
	echo -e "\n\n#####    Step$step: Exclude duplicates    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step $step"
		exit 0;
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.duplicates ]; then
			rm -rf $RunDir/$step.duplicates > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.duplicates
	fi
	if [ ! -d $RunDir/$step.duplicates ]; then
		echo "Error: subfolder not found: $RunDir/$step.duplicates" >&2
		exit 10;
	fi
	cd $RunDir/$step.duplicates
#Default
	declare -a listfastq8all=();
	declare -a listfastq8dedup=();
	fastqR1str=''
	fastqR2str=''
#input
	if [ ! -s "$step7bamout" ]; then
		echo "Error: last step output not found at step$step" >&2
		exit 10;
	fi
#output

#run
	bam2merge=''
	
	for (( i=0; i < ${#libraryarr[@]}; i++ )); do
		echo "    Library $i ${libraryarr[$i]}"
		echo "        FastQ : R1: ${fastq1arr[$i]}"
		echo "                R2: ${fastq2arr[$i]}"
		echo "    Library $i ${libraryarr[$i]}" >&2
		echo "        FastQ : R1: ${fastq1arr[$i]}" >&2
		echo "                R2: ${fastq2arr[$i]}" >&2
		if [ $startstep -le $step ]; then
			echo "Info: extract IDs for readgroup ${libraryarr[$i]}"
			samtools view -r "${libraryarr[$i]}" $step7bamout | perl -lane 'BEGIN{%idhash=()};unless (exists $idhash{$F[0]}) {print $F[0]; $idhash{$F[0]}++;}' > "$RunDir/$step.duplicates/${libraryarr[$i]}.id"
			if [ $? -ne 0 ]; then
				echo "Error: extract read IDs $RunDir/$step.duplicates/${libraryarr[$i]}.id" >&2
				exit 10;
			fi
			if [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.id" ]; then
				echo "Warnings: no read IDs $RunDir/$step.duplicates/${libraryarr[$i]}.id" >&2
				continue
			else
				echo "Info: extract R1 reads ${libraryarr[$i]}"
				if SeqTkSubSeqFastq "${fastq1arr[$i]}" "$RunDir/$step.duplicates/${libraryarr[$i]}.id" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.fastq"; then
					echo "Info: seqtk subseq R1 fastq success"
				else
					echo "Error: seqtk subseq R1 fastq error" >&2
					echo "CMD used: seqtk subseq ${fastq1arr[$i]} $RunDir/$step.duplicates/${libraryarr[$i]}.id > $RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.fastq" >&2
					exit 10;
				fi
				perl -ne 'chomp;s/\s+.*$/\/1/;print $_, "\n";$_=<>;print;$_=<>;print;$_=<>;print;' "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.fastq" > "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq"
				if [ $? -ne 0 ] || [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq" ];then
					echo "Error: R1 convert ID failed" >&2
					exit 10;
				fi
				
				echo "Info: extract R2 reads ${libraryarr[$i]}"
				if SeqTkSubSeqFastq "${fastq2arr[$i]}" "$RunDir/$step.duplicates/${libraryarr[$i]}.id" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.fastq"; then
					echo "Info: seqtk subseq R2 fastq success"
				else
					echo "Error: seqtk subseq R2 fastq error" >&2
					echo "CMD used: seqtk subseq ${fastq2arr[$i]} $RunDir/$step.duplicates/${libraryarr[$i]}.id > $RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.fastq" >&2
					exit 10;
				fi
				perl -ne 'chomp;s/\s+.*$/\/2/;print $_, "\n";$_=<>;print;$_=<>;print;$_=<>;print;' "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.fastq" > "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq"
				if [ $? -ne 0 ] || [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq" ];then
					echo "Error: R2 convert ID failed" >&2
					exit 10;
				fi
				
				fastq_checkid.pl "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq"  "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq"  '^\@(\S+)\/[12]{1}$'
				if [ $? -eq 0 ]; then
					echo "${libraryarr[$i]} paired"
#					touch "${libraryarr[$i]}_paired"
				else
					echo "Error: ${libraryarr[$i]} un-paired" >/dev/stderr
					touch "${libraryarr[$i]}_unpaired"
					exit 1;
				fi
				echo "Info: mapping as mates ${libraryarr[$i]}"
				if BowtieMappingMates "$opt_x" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq" " -q -p $threads -v $numgap -a --sam -X $opt_m " "$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.st.bam" > $RunDir/$step.duplicates/${libraryarr[$i]}.bowtiemates.log 2> $RunDir/$step.duplicates/${libraryarr[$i]}.bowtiemates.err; then
					echo "Info: bowtie mapping mates succeeds"
					rm "$RunDir/$step.duplicates/${libraryarr[$i]}.id" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq" "$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq" > / dev/null 2>&1
					tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.rename.fastq")
					tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.rename.fastq")
				else
					echo "Error: bowtie mapping mates error: $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.st.bam" >&2
					exit 10;
				fi
				bam2merge=" $bam2merge -b $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.st.bam -s LIB -r ${libraryarr[$i]} "
			fi
		fi
		listfastq8all+=("$RunDir/$step.duplicates/${libraryarr[$i]}.id")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.id.R1.subseq.fastq")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.id.R2.subseq.fastq")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.bowtiemates.log")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.bowtiemates.err")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.st.bam")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.st.bam.bai")
	done
	step7bamout=' '
	if [ $startstep -le $step ]; then
		echo "Info: bamaddrg merge"
		if BamAddGroup "$bam2merge" "$RunDir/$step.duplicates/subseq.st.merge.bam"; then
			echo "Info: bamaddrg merge success"
		else
			echo "Error: bamaddrg merge error" >&2
			exit 10;
		fi
		echo "Info: picard MarkDuplicates mark duplicated pairs"
		picardrmdup -d -k $RunDir/$step.duplicates/subseq.st.merge.bam > $RunDir/$step.duplicates/subseq.st.merge.rmdup.bam.log 2> $RunDir/$step.duplicates/subseq.st.merge.rmdup.bam.err
		if [ $? -ne 0 ] || [ ! -s $RunDir/$step.duplicates/subseq.st.merge.rmdup.bam ]; then
			echo "Error: picardrmdup: $RunDir/$step.duplicates/subseq.st.merge.bam" >&2
			exit 10;
		fi
	fi
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.bam")
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.bam.bai")
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.rmdup.bam.log")
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.rmdup.bam.err")
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.rmdup.bam")
	tempfiles+=("$RunDir/$step.duplicates/subseq.st.merge.rmdup.bam.bai")
	bam3merge1str=''
	bam3merge2str=''
	bam3merge1arr=()
	bam3merge2arr=()
	for (( i=0; i < ${#libraryarr[@]}; i++ )); do
		outputprefix="${libraryarr[$i]}.fourthBT"
		echo "    Library $i ${libraryarr[$i]}"
		echo "        FastQ : R1: ${fastq1arr[$i]}"
		echo "                R2: ${fastq2arr[$i]}"
		echo "    Library $i ${libraryarr[$i]}" >&2
		echo "        FastQ : R1: ${fastq1arr[$i]}" >&2
		echo "                R2: ${fastq2arr[$i]}" >&2
		if [ $startstep -le $step ]; then
			echo "Info: extract read names with FLAG 1024"
			samtools view -f 1024 -r "${libraryarr[$i]}" $RunDir/$step.duplicates/subseq.st.merge.rmdup.bam | perl -lane 'BEGIN{%idhash=()};unless (exists $idhash{$F[0]}) {print $F[0]; $idhash{$F[0]}++;}' > $RunDir/$step.duplicates/${libraryarr[$i]}.exclude.id
			if [ $? -ne 0 ]; then
				echo "Error: extract read IDs $RunDir/$step.duplicates/${libraryarr[$i]}.exclude.id" >&2
				exit 10;
			fi
			uniquelist=''
			if [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.id" ]; then
				if [ -s $RunDir/$step.duplicates/${libraryarr[$i]}.exclude.id ]; then
					echo "Error: no read IDs but have IDs to exclude" >&2
					exit 10;
				fi
			elif [ -s "$RunDir/$step.duplicates/${libraryarr[$i]}.id" ]; then
				echo "Info: exclude the duplicated read IDs"
				list2_compar.pl "$RunDir/$step.duplicates/${libraryarr[$i]}.id" \
				"$RunDir/$step.duplicates/${libraryarr[$i]}.exclude.id" \
				"$RunDir/$step.duplicates/${libraryarr[$i]}.share" \
				"$RunDir/$step.duplicates/${libraryarr[$i]}.u1" \
				"$RunDir/$step.duplicates/${libraryarr[$i]}.u2"
				if [ $? -ne 0 ] || [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.u1" ]; then
					echo "Warnings: list2_compar.pl $RunDir/$step.duplicates/${libraryarr[$i]}.u1" >&2
#					exit 10;
				fi
				rm $RunDir/$step.duplicates/${libraryarr[$i]}.exclude.id > /dev/null 2>&1
				uniquelist=$RunDir/$step.duplicates/${libraryarr[$i]}.u1
			fi
			if [ ! -z "$uniquelist" ] || [ -s $uniquelist ]; then
				echo "Info: final list: $uniquelist"
				echo "Info: final list: $uniquelist" >&2
				echo "Info: extract unique R1 read IDs"
				seqtk subseq ${fastq1arr[$i]} $uniquelist > $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq
				if [ $? -ne 0 ] || [ ! -s $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq ]; then
					echo "Error2: seqtk subseq R1 fastq error" >&2
					echo "CMD used: seqtk subseq ${fastq1arr[$i]} $RunDir/$step.duplicates/${libraryarr[$i]}.u1 > $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq" >&2
					exit 10;
				fi
				
				echo "Info: extract unique R2 read IDs"
				seqtk subseq ${fastq2arr[$i]} $uniquelist > $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq
				if [ $? -ne 0 ] || [ ! -s $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq ]; then
					echo "Error2: seqtk subseq R2 fastq error" >&2
					echo "CMD used: seqtk subseq ${fastq2arr[$i]} $RunDir/$step.duplicates/${libraryarr[$i]}.u1 > $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq" >&2
					exit 10;
				fi
				
				echo "Info: bowtie mapping unique R1 read IDs"
				if BowtieMappingSingle "$opt_x" "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq" " -q -p $threads -v $numgap -a --sam" "$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R1.final.st.bam" > $RunDir/$step.duplicates/${libraryarr[$i]}.R1.bowtie.log 2> $RunDir/$step.duplicates/${libraryarr[$i]}.R1.bowtie.err; then
					echo "Info: bowtie mapping succeeds"
				else
					echo "Error: bowtie mapping error: $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R1.final.st.bam" >&2
					exit 100;
				fi
				echo "Info: bowtie mapping unique R2 read IDs"
				if BowtieMappingSingle "$opt_x" "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq" " -q -p $threads -v $numgap -a --sam" "$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R2.final.st.bam"  > $RunDir/$step.duplicates/${libraryarr[$i]}.R2.bowtie.log 2> $RunDir/$step.duplicates/${libraryarr[$i]}.R2.bowtie.err; then
					echo "Info: Bowtie mapping succeeds"
				else
					echo "Error: bowtie mapping error: $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R2.final.st.bam" >&2
					exit 100;
				fi
				bam3merge1str="$bam3merge1str -b $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R1.final.st.bam -s LIB -r ${libraryarr[$i]}"
				bam3merge2str="$bam3merge2str -b $RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R2.final.st.bam -s LIB -r ${libraryarr[$i]}"
				gzip -9 "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq"
				if [ $? -ne 0 ] || [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq.gz" ]; then
					echo "Error: Fastq R1 gzip error: $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq" >&2
					exit 100
				fi
				fastqR1str="$fastqR1str $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R1.subseq.fastq.gz "
				gzip -9 "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq"
				if [ $? -ne 0 ] || [ ! -s "$RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq.gz" ]; then
					echo "Error: Fastq R2 gzip error: $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq" >&2
					exit 100;
				fi
				fastqR2str="$fastqR2str $RunDir/$step.duplicates/${libraryarr[$i]}.u1.id.R2.subseq.fastq.gz "
			else
				echo "Warnings: no read IDs to extract" >&2
			fi
		fi
		listfastq8dedup+=("$RunDir/$step.duplicates/${libraryarr[$i]}.u1")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.share")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.u2")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R1.final.st.bam")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R1.final.st.bam.bai")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R2.final.st.bam")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.subseq.R2.final.st.bam.bai")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.R1.bowtie.log")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.R1.bowtie.err")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.R2.bowtie.log")
		tempfiles+=("$RunDir/$step.duplicates/${libraryarr[$i]}.R2.bowtie.err")
	done
	
	step8bamout1="$RunDir/$step.duplicates/$opt_p.8.subseq.R1.st.merge.final.bam"
	step8bamout2="$RunDir/$step.duplicates/$opt_p.8.subseq.R2.st.merge.final.bam"
	step8bamfinal="$RunDir/$step.duplicates/$opt_p.8.subseq.R1R2.st.merge.final.bam"
	
	if [ $startstep -le $step ]; then
		echo "Info: bamaddrg merge BAMs for unique R1 read IDs"
		bamaddrg $bam3merge1str > $step8bamout1
		if [ $? -ne 0 ] || [ ! -s $step8bamout1 ]; then
			echo "Error: bamaddrg merge R1" >&2
			echo "CMD used: bamaddrg $bam3merge1str > $step8bamout1" >&2
			exit 10;
		fi
		if IndexBam "$step8bamout1"; then
			echo "Info: index BAMs"
		else
			echo "Error: Failed to index" >&2
			exit 10;
		fi
		echo "Info: bamaddrg merge BAMs for unique R2 read IDs"
		bamaddrg $bam3merge2str > $step8bamout2
		if [ $? -ne 0 ] || [ ! -s $step8bamout2 ]; then
			echo "Error: bamaddrg merge R2" >&2
			echo "CMD used: bamaddrg $bam3merge2str > $step8bamout2" >&2
			exit 10;
		fi
		if IndexBam "$step8bamout2"; then
			echo "Info: index BAMs"
		else
			echo "Error: Failed to index" >&2
			exit 10;
		fi
#		echo "Info: samtools merge BAMs for unique R1R2 read IDs"
#		samtools merge $step8bamfinal $step8bamout1 $step8bamout2
#		if [ $? -ne 0 ] || [ ! -s $step8bamfinal ]; then
#			echo "Error: samtools merge R1R2" >&2
#			echo "CMD used: samtools merge $step8bamfinal $step8bamout1 $step8bamout2" >&2
#			exit 10;
#		fi
#		if IndexBam "$step8bamfinal"; then
#			echo "Info: index BAMs"
#		else
#			echo "Error: Failed to index" >&2
#			exit 10;
#		fi
		bam_scaffolding_separately.stats_noBioDBSAM.pl $step8bamout1 $step8bamout2 $bam_pairedstats
		if [ $? -ne 0 ] || [ ! -s "$listfastq8dedup" ]; then
			echo "Warnings: failed to stats paired BAMs" >&2
			echo "CMD used: bam_scaffolding_separately.stats_noBioDBSAM.pl $step8bamout1 $step8bamout2 $bam_pairedstats"
		fi
		
		zcat $fastqR1str | gzip -9 > "$finalcleanR1"
		if [ $? -ne 0 ] || [ ! -s "$finalcleanR1" ]; then
			echo "Error: creating final clean fastq R1: $finalcleanR1" >&2
			echo "CMD used: cat $fastqR1str | gzip -9 > $finalcleanR1" >&2
			exit 100;
		fi
		zcat $fastqR2str | gzip -9 > "$finalcleanR2"
		if [ $? -ne 0 ] || [ ! -s "$finalcleanR2" ]; then
			echo "Error: creating final clean fastq R2: $finalcleanR2" >&2
			echo "CMD used: cat $fastqR2str | gzip -9 > $finalcleanR2" >&2
			exit 100;
		fi
	fi
	echo -e "\n\n\n###Read number after pair-up"
	for indlistfile in ${listfastq8all[@]}; do
		if [ -e "$indlistfile" ]; then
			wc -l "$indlistfile"
		fi
	done
	listfastq8all=();
	echo -e "\n\n\n###Read number after pair dedup"
	for indlistfile in ${listfastq8dedup[@]}; do
		if [ -e "$indlistfile" ]; then
			wc -l "$indlistfile"
		fi
	done
	listfastq8dedup=();
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi
	
	

### Step9
	step=9
	echo -e "\n\n#####    Step$step: ReadCleaning    #####"
	echo -e "\n\n#####    Step$step: ReadCleaning    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step $step"
		exit 0;
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.cleanread ]; then
			rm -rf $RunDir/$step.cleanread > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.cleanread
	fi
	if [ ! -d $RunDir/$step.cleanread ]; then
		echo "Error: subfolder not found: $RunDir/$step.duplicates" >&2
		exit 10;
	fi
	cd $RunDir/$step.cleanread
#input
	if [ $startstep -le $step ]; then
		if [ ! -e "$bam_unpairedstats" ]; then
			echo "Error: step 8 output not found at step$step: $bam_unpairedstats" >&2
			exit 10;
		fi

		if [ ! -s "$bam_pairedstats" ]; then
			echo "Error: step 8 output not found at step$step: $bam_pairedstats" >&2
			exit 10;
		fi
	
		if [ ! -s $finalcleanR1 ]; then
			echo "Error: step 8 output not found at step$step: $finalcleanR1" >&2
			exit 10;
		fi
		if [ ! -s $finalcleanR2 ]; then
			echo "Error: step 8 output not found at step$step: $finalcleanR2" >&2
			exit 10;
		fi
	fi
#run
	if [ $startstep -le $step ]; then
		if [ -z "$opt_ref" ] && [ -s "$opt_ref" ]; then
			echo "Error: reference -r not specified"
		
		fi
		if [ ! -s "$opt_ref.fai" ]; then
			samtools faidx "$opt_ref"
			if [ $? -ne 0 ] || [ ! -s "$opt_ref.fai" ]; then
				echo "Error: indexing fasta at step$step: $opt_ref" >&2
				echo "CMD used samtools faidx $opt_ref" >&2
				exit 10;
			fi
		fi
		bam_exclude_reads_by_windows_size_and_numpairs.pl $bam_pairedstats $opt_ref.fai $window_size1 $window_size2 $num_read_in_window $opt_p.9
		if [ $? -ne 0 ] || [ ! -s "$opt_p.9.finalreadids" ] || [ ! -s "$opt_p.9.finalrefsids" ]; then
			echo "Error: bam_exclude_reads_by_windows_size_and_numpairs.pl running error" >&2
			echo "CMD used: bam_exclude_reads_by_windows_size_and_numpairs.pl $bam_pairedstats $opt_ref.fai $window_size1 $window_size2 $num_read_in_window $opt_p.9" >&2
			exit 100;
		fi
		mv "$opt_p.9.finalrefsids" $finalrefsid
		mv "$opt_p.9.finalreadids" $finalreadid
	fi
#	if RunStatsPipeline $opt_ref.fai $bam_unpairedstats $bam_pairedstats $meaninsert $stdevinsert $prefix; then
#		echo "$prefix succeeded"
#	else
#		echo "Error: $prefix failed"
#	fi
	if [ $startstep -le $step ]; then
		if SeqTkSubSeqFasta $opt_ref $finalrefsid $finalrefsfa; then
			echo "Info: seqtk subseq fasta success"
		else
			echo "Error: seqtk subseq fasta error" >&2
			echo "CMD used: SeqTkSubSeqFasta $opt_ref $finalrefsid $finalrefsfa" >&2
			exit 10;
		fi
		if SeqTkSubSeqFastq $finalcleanR1 $finalreadid $finalreadfqR1; then
			echo "Info: seqtk subseq R1 fastq success"
		else
			echo "Error: seqtk subseq R1 fastq error" >&2
			exit 10;
		fi
		if SeqTkSubSeqFastq $finalcleanR2 $finalreadid $finalreadfqR2; then
			echo "Info: seqtk subseq R2 fastq success"
		else
			echo "Error: seqtk subseq R2 fastq error" >&2
			exit 10;
		fi
	fi
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi



### step 10
	step=10
	echo -e "\n\n#####    Step$step: sspace    #####"
	echo -e "\n\n#####    Step$step: sspace    #####" >&2
	cd $RunDir
	if [ $endstep -lt $step ]; then
		echo "Info: debug mode ends before step $step"
		exit 0;
	fi
	if [ $startstep -le $step ]; then
		if [ -d $RunDir/$step.sspace ]; then
			rm -rf $RunDir/$step.sspace > /dev/null 2>&1
		fi
		mkdir -p $RunDir/$step.sspace
	fi
	if [ ! -d $RunDir/$step.sspace ]; then
		echo "Error: subfolder not found: $RunDir/$step.duplicates" >&2
		exit 10;
	fi
	cd $RunDir/$step.sspace
#input
	if [ $startstep -le $step ]; then
		if [ ! -s $finalreadfqR1 ]; then
			echo "Error: step 8 output not found at step$step: $finalreadfqR1" >&2
			exit 10;
		fi
		if [ ! -s $finalreadfqR2 ]; then
			echo "Error: step 8 output not found at step$step: $finalreadfqR2" >&2
			exit 10;
		fi
		if [ ! -s $finalrefsfa ]; then
			echo "Error: step 8 output not found at step$step: $finalrefsfa" >&2
			exit 10;
		fi
	fi
#output
	sspace_config="$RunDir/$step.sspace/$opt_p.sspace.config"
### Run
	if [ $startstep -le $step ]; then
		echo -e "Lib1\tbowtie\t$finalreadfqR1\t$finalreadfqR2\t$mean_insertsize\t$stdev_insertsize\t$strand_orientation" > $sspace_config
		if [ $? -ne 0 ] || [ ! -s "$sspace_config" ]; then
			echo "Error: creating SSPACE config failed" >&2
			exit 10;
		fi
		if RunSspace "$opt_ref" $finalrefsfa $sspace_config $finalmerge; then
			echo "Info: SSPACE running succeeded"
		else
			echo "Error: SSPACE running failed"
			exit 100;
		fi
	fi
	
	if [ $opt_d -eq 1 ]; then
		CleanTemp
	fi
fi



exit 0;
