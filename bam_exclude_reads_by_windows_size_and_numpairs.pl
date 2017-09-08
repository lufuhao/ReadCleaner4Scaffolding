#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/sum min max/;
use Data::Dumper qw/Dumper/;
use constant USAGE =><<EOH;

usage: $0    bam.stats    fasta.fai    window_size1    window_size2    num_pairs    output_prefix

v20170228

Requirements
    List::Util; Data::Dumper

Input and Output:
    bam.stats:       bam_scaffolding.stats_noBioDBSAM.pl output_prefix
    fasta.fai:       fasta index
    length/pairs:    how many pairs in such a length for a scaffolding
    output:          readname.output  .oneonly .window.pairs .filterout

Examples: 

bam_exclude_reads_by_windows_size_and_numpairs.pl ../MP40Kb7lib.subseq.R1R2.st.merge.final.bam.stats ../../ta3bAllScaffoldsV443.genom.fa.fai 10000 5 testout

# Window distance
cut -f 11,13 testout.distance | perl -ne 'chomp; \@arr=split(/\t/);next unless (scalar(\@arr)==2 and \$arr[0] =~/^\\d+\$/); print \$arr[0], "\\t", \$arr[1],"\\n";' > log

# strand
grep -v ^'#' testout.distance | cut -f 2,5,8 | sort -u > log
perl -ne 'chomp; \@arr=split(/\\t/);if (exists \$hash{\$arr[0]} and \$hash{\$arr[0]}{\$arr[1]}) {die "Error: repeat: \$_ \\n";} if (exists \$hash{\$arr[1]} and \$hash{\$arr[1]}{\$arr[0]}) {next;} \@arr2=split(/\\\//, \$arr[2]); print \$arr2[0], "\\t", \$arr2[1], "\\n"; \$hash{\$arr[0]}{\$arr[1]}++;' log

EOH
die USAGE if (scalar(@ARGV) !=6 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $bamstats=$ARGV[0];
my $fastaindex=$ARGV[1];
my $windows1=$ARGV[2];
my $windows2=$ARGV[3];
my $numpairs=$ARGV[4];
my $output=$ARGV[5];



### Configurable start ###
my $maxmappingtimes=1;
my $ratiothreshold=5;
my $min_allowed_insert=20000;
my $max_allowed_insert=60000;
### Configurable end ###



die "Error: invalid bam.stats input\n" unless (defined $bamstats and -s $bamstats);

die "Error: invalid fasta.fai input\n" unless (defined $fastaindex and -s $fastaindex);

die "Error: invalid window_size\n" unless (defined $windows1);

if ($windows1=~/^\d+[kK]$/) {
	print "Info: window_size in Kilobase accepted\n";
	$windows1=~s/[kK]$/000/;
}
elsif ($windows1=~/^\d+[mM]$/) {
	print "Info: window_size in megabase accepted\n";
	$windows1=~s/[mM]$/000000/;
}
elsif ($windows1=~/^\d+[gG]$/) {
	print "Info: window_size in gigabase accepted\n";
	$windows1=~s/[gG]$/000000000/;
}
die "Error: invalid window_size1: $windows1\n" unless ($windows1 =~/^\d+$/);
if ($windows2=~/^\d+[kK]$/) {
	print "Info: window_size in Kilobase accepted\n";
	$windows2=~s/[kK]$/000/;
}
elsif ($windows2=~/^\d+[mM]$/) {
	print "Info: window_size in megabase accepted\n";
	$windows2=~s/[mM]$/000000/;
}
elsif ($windows2=~/^\d+[gG]$/) {
	print "Info: window_size in gigabase accepted\n";
	$windows2=~s/[gG]$/000000000/;
}
die "Error: invalid window_size2: $windows2\n" unless ($windows2 =~/^\d+$/);
die "Error: invalid num_pairs [only numbers 0-9 are accepted, integer]\n" unless (defined $numpairs and $numpairs=~/^\d+$/);

die "Error: invalid output\n" unless (defined $output);

unlink "$output.oneonly" if (-e "$output.oneonly");
unlink "$output.window$windows1-$windows2.pairs$numpairs" if (-e "$output.window$windows1-$windows2.pairs$numpairs");
unlink "$output.filterout" if (-e "$output.filterout");
unlink "$output.distance" if (-e "$output.distance");
unlink "$output.overlap" if (-e "$output.overlap");
unlink "$output.scaffold.pairs" if (-e "$output.scaffold.pairs");
unlink "$output.finalreadids" if (-e "$output.finalreadids");
unlink "$output.finalrefsids" if (-e "$output.finalrefsids");



### Default
my %strand=();
my $strness='NaN';
my $ratio='NaN';
my %finalarray=();
my %idhash=();
my $linenum=0;
my %seqlength=();
my $testwindows=0;
my $test_overlap=0;
my @insertsum=();
my $insertnum=0;
my $estimatedinsertsize=0;
my $estimatedgap=0;
my $maxinsert=0;
my %read2gap=();
my $readproblematic=0;
my $printthislink=0;
my %printreads=();
my %printrefs=();

### Read sequence length
open (FASTAINDEX, " < $fastaindex") || die "Error: can not open fasta.fai\n";
while (my $line=<FASTAINDEX>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	unless ($arr[0]=~/^\S+$/ and $arr[1]=~/^\d+$/) {
		die "Error: invalid fasta line: $line\n";
	}
	if (exists $seqlength{$arr[0]}) {
		die "Error: repeated sequence ID: $arr[0]\n";
	}
	$seqlength{$arr[0]}=$arr[1];
}
close FASTAINDEX;
print "\n##### SUMMARY 1 #####\n";
print "    Read total fasta.fai lines: $linenum\n";
print "    Total seqeunces:            ", scalar(keys %seqlength), "\n\n";




open (BAMSTATSIN, " < $bamstats") || die "Error: can not open bam.stats: $bamstats\n";
open (FILTERED, " > $output.window$windows1-$windows2.pairs$numpairs") || die "Error: can not write filtered reads: $output.window$windows1-$windows2.pairs$numpairs\n";
open (FILTEROUT, " > $output.filterout") || die "Error: can not write filter_out reads: $output.filterout\n";
open (DISTANCE, " > $output.distance") || die "Error: can not write distance reads: $output.distance\n";
print DISTANCE "#No/Total-NumReads\tRef1name\tRef1len\tRef1std\tRef2name\tRef2len\tRef2std\tS/D\tratio\tAssemblyStd\tDist1\tDist1Detail\tDist2\tDist2Detail\tTest_overlap\tMeaninternalSize\tDetail\n";
open (OVERLAP, " > $output.overlap") || die "Error: can not write contig overlap: $output.overlap\n";
open (SCAFPAIR, " > $output.scaffold.pairs") || die "Error: can not write scaffold.pairs: $output.scaffold.pairs\n";
print SCAFPAIR "Scaffold1\tScaffold2\tST/RC\tDirection\tNum_Fit_Grooup\tNum_Group\tPerc_Fit\tOverlap\tMeanDist\tStdevDist\tMeanGap\tStdevGap\n";
open (FINALREADLIST, " > $output.finalreadids") || die "Error: can not write finalreadids: $output.finalreadids\n";
open (FINALREFSLIST, " > $output.finalrefsids") || die "Error: can not write finalrefsids: $output.finalrefsids\n";

$linenum=0;
while (my $line=<BAMSTATSIN>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
#	$arr[0]=quotemeta($arr[0]);
	if (exists $idhash{$arr[0]}{$arr[2]}{$arr[1]}) {
		$idhash{$arr[0]}{'EXCLUDEDPLZ'}++;
		next;
	}
	if ($arr[5] eq '+') {
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'pos'}=$arr[3];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'str'}=$arr[5];
	}
	elsif ($arr[5] eq '-') {
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'pos'}=$arr[4];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'str'}=$arr[5];
	}
	else {
		die "Error: strand error at line($linenum): $line\n";
	}
}
close BAMSTATSIN;



open (ONLYONE, "> $output.oneonly") || die "Error: can not write onlyonly: $output.oneonly";
my %links=();
my %read2links=();
my @allreads=();
my $unique_mapped_read_num=0;
@allreads=keys (%idhash);

foreach my $indread (@allreads) {
	my @arr=();
	if (exists $idhash{$indread}{'EXCLUDEDPLZ'}) {
		print FILTEROUT $indread, "\tMappOneReferenceTwice\n";
		delete $idhash{$indread};
		next;
	}
	@arr=keys %{$idhash{$indread}};
	
	
	my $excludethisread=0;
	foreach my $indctg (@arr) {### Control read mapped to different strand and -pos is downstream of +pos
		if (exists $idhash{$indread}{$indctg}{1} and exists $idhash{$indread}{$indctg}{2}) {
			if ($idhash{$indread}{$indctg}{1}{'str'} eq '+' and $idhash{$indread}{$indctg}{2}{'str'} eq '-') {
				unless ($idhash{$indread}{$indctg}{2}{'pos'}> $idhash{$indread}{$indctg}{1}{'pos'}) {
					$excludethisread=1;
					last;
				}
			}
			elsif ($idhash{$indread}{$indctg}{1}{'str'} eq '-' and $idhash{$indread}{$indctg}{2}{'str'} eq '+') {
				unless ($idhash{$indread}{$indctg}{1}{'pos'} > $idhash{$indread}{$indctg}{2}{'pos'}) {
					$excludethisread=1;
					last;
				}
			}
			else {
				$excludethisread=1;
				last;
			}
		}
	}
	if ($excludethisread==1) {
		print FILTEROUT $indread, "\tMappingProblem\n";
		delete $idhash{$indread};
		next;
	}
	
	if (scalar(@arr)==1) {### Control mapped to only one reference
		print ONLYONE $indread, "\t";
		if (exists $idhash{$indread}{$arr[0]} and exists $idhash{$indread}{$arr[0]}{2} and exists $idhash{$indread}{$arr[0]}{1}) {
			my $insetsize=abs($idhash{$indread}{$arr[0]}{2}{'pos'}-$idhash{$indread}{$arr[0]}{1}{'pos'})+1;
			if ($insetsize>=$min_allowed_insert and $insetsize<=$max_allowed_insert) {
				push (@insertsum, $insetsize);
				$insertnum++;
				print ONLYONE $insetsize, "\n";
				print FINALREADLIST $indread, "\n";
				$unique_mapped_read_num++;
			}
		}
		else {
			print ONLYONE "NaN\n";
		}
		delete $idhash{$indread};
		next;
	}
	
	my @arr2=();
	my @arr3=();
	foreach my $contig (@arr) {
		push (@arr2, $contig) if (exists $idhash{$indread}{$contig}{1});
		push (@arr3, $contig) if (exists $idhash{$indread}{$contig}{2});
	}
	if (scalar(@arr2)<1 or scalar(@arr3)<1) {### Remove Only one mate mapped
		print FILTEROUT $indread, "\tOnlyOneMateMapped\n";
		delete $idhash{$indread};
		next;
	}
	my $maxmapping=max (scalar(@arr2), scalar(@arr3));
	if ($maxmapping>$maxmappingtimes) {### control mapping times
		print FILTEROUT $indread, "\tMappedMoreThan", $maxmappingtimes, "Times$maxmapping\n";
		delete $idhash{$indread};
		next;
	}
	foreach my $first (@arr2) {
		foreach my $second (@arr3) {
			next if ($first eq $second);
#			next if (exists $links{$second} and exists $links{$second}{$first});
			$links{$first}{$second}{$indread}++;
			$read2links{$indread}{$first}{$second}++;
			$links{$second}{$first}{$indread}++;
			$read2links{$indread}{$second}{$first}++;
		}
	}
}

#print "Test: \%links\n"; print Dumper \%links; print "\n"; exit 0; ### For test ###
#print "Test: \%read2links\n"; print Dumper \%read2links; print "\n"; exit 0; ### For test ###

close ONLYONE;
print "### Summary ###\n";
print "Num: $insertnum\n";
$estimatedinsertsize=int (&CalAverage(@insertsum));
my $stdev_insertsize=int (CalStdev(@insertsum));
print "Estimated insert size average: ", $estimatedinsertsize, "\n";
print "Estimated insert size stdev : ",$stdev_insertsize , "\n";
$maxinsert=$estimatedinsertsize+$stdev_insertsize;
print "Maximum insert size : ", $maxinsert, "\n";
@insertsum=();
@allreads=();


my @firsts=();
my %scafstrandness=();
@firsts=keys %links;
#print "Test: \n"; print "Total number first strands: ", scalar(@firsts), "\n\n"; ### For test ###
foreach my $first (@firsts) {
	next unless (exists $links{$first});
	my @seconds=keys %{$links{$first}};
	delete $links{$first} unless (scalar(@seconds)>0);
#	print "SCAFFOLD1 $first\n"; ### For test ###
	foreach my $second (@seconds) {
#		print "\tSCAFFOLD2 $second\n"; ### For test ###
		next unless (exists $links{$first}{$second});
		my @readid=();
		@readid=keys %{$links{$first}{$second}};
#		print "SCAFFOLD1 $first SCAFFOLD2 $second Read number: ", scalar(@readid), "\n"; ### For test ###
		unless (scalar(@readid)>=$numpairs) {
#			print STDERR "Info: pairs not have enough evidence: SCAFFOLD1 $first SCAFFOLD2 $second\n"; ### For test ###
			delete $links{$first}{$second} if (exists $links{$first} and exists $links{$first}{$second});
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} unless (scalar(keys %{$links{$second}})>0);
			next;
		}
#		print "SCAFFOLD1 $first SCAFFOLD2 $second Read number: ", scalar(@readid), "\n"; ### For test ###
		
		
		
		### Decide strand
		%strand=('s' => 0, 'd' => 0);
		$strness='NaN';
		$ratio='NaN';
		my %readstrand=();
		foreach my $indread (@readid) {
			my @firstpos=();
			my @secondpos=();
			unless (exists $idhash{$indread}) {
				&RemoveRead($indread);
				next;
			}
			$readproblematic=0;
			&EstimateOverlap($indread, $first, $second);
			if ($readproblematic==1) {
				print FILTEROUT $indread, "\tProblematicMapping\n";
				&RemoveRead($indread);
				next;
			}
			
			my $firstmate=0;
			my $secondmate=0;
			if (exists $idhash{$indread}{$first}) {
				if (exists $idhash{$indread}{$first}{1}) {
					push (@firstpos, 1);
					push (@firstpos, $idhash{$indread}{$first}{1}{'pos'});
					push (@firstpos, $idhash{$indread}{$first}{1}{'str'});
					$firstmate++;
				}
				if (exists $idhash{$indread}{$first}{2}) {
					push (@secondpos, 2);
					push (@secondpos, $idhash{$indread}{$first}{2}{'pos'});
					push (@secondpos, $idhash{$indread}{$first}{2}{'str'});
					$secondmate++;
				}
			}
			if (exists $idhash{$indread}{$second}) {
				if (exists $idhash{$indread}{$second}{1}) {
					push (@secondpos, 1);
					push (@secondpos, $idhash{$indread}{$second}{1}{'pos'});
					push (@secondpos, $idhash{$indread}{$second}{1}{'str'});
					$secondmate++;
				}
				if (exists $idhash{$indread}{$second}{2}) {
					push (@firstpos, 2);
					push (@firstpos, $idhash{$indread}{$second}{2}{'pos'});
					push (@firstpos, $idhash{$indread}{$second}{2}{'str'});
					$firstmate++;
				}
			}
#			print "Read: $indread\t$first vs $second\n"; ### For test ###
#			print "Test: \@firstpos: ", join ("\t", @firstpos), "\n"; ### For test ###
#			print "Test: \@secondpos: ", join ("\t", @secondpos), "\n"; ### For test ###
			if ($firstmate==2 or $secondmate==2) {
				if ($firstmate==2 and $secondmate==2) {
					if ($firstpos[2] eq $firstpos[5]) {
						$strand{'s'}++; $readstrand{$indread}{'s'}++;
					}
					else {
						$strand{'d'}++;$readstrand{$indread}{'d'}++;
					}
					if ($secondpos[2] eq $secondpos[5]) {
						$strand{'s'}++;$readstrand{$indread}{'s'}++;
					}
					else {
						$strand{'d'}++;$readstrand{$indread}{'d'}++;
					}
				}
				elsif ($firstmate==2 and $secondmate<2) {
					if ($firstpos[2] eq $firstpos[5]) {
						$strand{'s'}++;$readstrand{$indread}{'s'}++;
					}
					else {
						$strand{'d'}++;$readstrand{$indread}{'d'}++;
					}
				}
				elsif ($firstmate<2 and $secondmate==2) {
					if ($secondpos[2] eq $secondpos[5]) {
						$strand{'s'}++;$readstrand{$indread}{'s'}++;
					}
					else {
						$strand{'d'}++;$readstrand{$indread}{'d'}++;
					}
				}
				else  {
					print "Warnings: can not understand: $first ($firstmate) vs $second ($secondmate)\tRead: $indread\n\n";
				}
			}
			else {
				print STDERR "Warnings: not mapped to any: READ: $indread SCAFFOLD1 $first SCAFFOLD2 $second\n";
				next;
			}
		}
#		print "Test: FIRST $first SECOND $second\n"; print Dumper \%strand; print "\n"; ### for test ###




###checking strand
		if (($strand{'d'}+$strand{'s'})<$numpairs) {
			delete $links{$first}{$second} if (exists $links{$first} and exists $links{$first}{$second});
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} unless (scalar(keys %{$links{$second}})>0);
			next;
		}
		my $finalstrand='';
		if ($strand{'d'}==0 and $strand{'s'}>=$numpairs) {
			$strness='Reverse';
			$ratio=$strand{'s'};
		}
		elsif ($strand{'s'}==0 and $strand{'d'}>=$numpairs) {
			$strness='Same';
			$ratio=$strand{'d'};
		}
		elsif ($strand{'s'}>=$strand{'d'} and $strand{'d'}>0) {
			 {
				$strness='Reverse';
				$ratio=$strand{'s'}/$strand{'d'};
			}
		}
		elsif ($strand{'s'}<$strand{'d'} and $strand{'s'}>0) {
			if ($strand{'d'}/$strand{'s'}>=$ratiothreshold) {
				$strness='Same';
				$ratio=$strand{'d'}/$strand{'s'};
			}
		}
		else {
			print STDERR "Warnings: unknown strand count: FIRST $first SECOND $second Same:", $strand{'s'}, "  Diff", $strand{'d'}, "\n";
			print Dumper $links{$first}{$second};
		}
		if ($ratio>=$ratiothreshold) {
			if ($strand{'s'}>$strand{'d'}) {$finalstrand='s';} 
			if ($strand{'s'}<$strand{'d'}) {$finalstrand='d';}
		}
		unless ($finalstrand=~/^[sd]{1,1}$/) {
			delete $links{$first}{$second} if (exists $links{$first} and exists $links{$first}{$second});
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} unless (scalar(keys %{$links{$second}})>0);
			next;
		}
#		print "SCAFFOLD1 $first SCAFFOLD2 $second STRAND $finalstrand\n";### For test ###
		
		foreach my $indread (keys %readstrand) {
			if (exists $readstrand{$indread}{'s'} and exists $readstrand{$indread}{'d'}) {
				print FILTEROUT "$indread\tControversyStrandingbySelf\n";
				&RemoveRead($indread);
				next;
			}
			unless (exists $readstrand{$indread} and exists $readstrand{$indread}{$finalstrand}) {
				print FILTEROUT "$indread\tControversyStrandingWithOthers\n";
				&RemoveRead($indread);
				next;
			}
		}
		unless (exists $scafstrandness{$first} and exists $scafstrandness{$first}{$second}) {
			$scafstrandness{$first}{$second}=$finalstrand;
#			$scafstrandness{$second}{$first}=$finalstrand;
		}
		else {
			die "Error: repeated strand: SCAFFOLD1 $first SCAFFOLD2 $second STRAND $finalstrand\n";
		}
	}
}
#print "\n\n\n"; print Dumper %scafstrandness; print "\n\n\n"; exit 0;### For test ###



###Double check
@firsts=();
@firsts=keys %links;
foreach my $first (@firsts) {
	next unless (exists $links{$first});
	my @seconds=keys %{$links{$first}};
	delete $links{$first} unless (scalar(@seconds)>0);
	foreach my $second (@seconds) {
		next unless (exists $links{$first}{$second});
		my @readid=();
		@readid=keys %{$links{$first}{$second}};
		unless (scalar(@readid)>0) {
			delete $links{$first}{$second};
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} if (exists $links{$second} and scalar(keys %{$links{$second}})==0);
		}
		my $numlinks=0;
		foreach my $indread (@readid) {
			if (exists $idhash{$indread}) {
				$numlinks++;
			}
			else {
				delete $links{$first}{$second}{$indread};
			}
		}
		
		unless ($numlinks>=$numpairs) {
			print STDERR "Info: pairs not have enough evidence2: SCAFFOLD1 $first SCAFFOLD2 $second\n";
			delete $links{$first}{$second} if (exists $links{$first} and exists $links{$first}{$second});
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} if (exists $links{$second} and scalar(keys %{$links{$second}})==0);
			next;
		}
	}
}



@firsts=();
@firsts=keys %links;
foreach my $first (@firsts) {
	next unless (exists $links{$first});
	my @seconds=keys %{$links{$first}};
	delete $links{$first} unless (scalar(@seconds)>0);
	foreach my $second (@seconds) {
		next unless (exists $links{$first}{$second});
		my @readid=();
		@readid=keys %{$links{$first}{$second}};
		unless (scalar(@readid)>0) {
			delete $links{$first}{$second};
			delete $links{$first} unless (scalar(keys %{$links{$first}})>0);
			delete $links{$second}{$first} if (exists $links{$second} and exists $links{$second}{$first});
			delete $links{$second} if (exists $links{$second} and scalar(keys %{$links{$second}})==0);
		}
		$testwindows=0;
		unless (scalar(@readid)>=$numpairs) {
			print STDERR "Info: pairs not have enough evidence\n";
			next;
		}
		my $plusstart=0; my %finalplus=();
		my $minusstart=0; my %finalminus=();
		my $arrstart=0; %finalarray=();
		%strand=('s' => 0, 'd' => 0);
		$strness='NaN';
		$ratio='NaN';
		$test_overlap=0;
		$printthislink=0;
		%read2gap=();
		my @disance=();
		my @gaps=();
		foreach my $indread (@readid) {
			my @firstpos=();
			my @secondpos=();
			die "Error: Non-existed2 readname: $indread\n" unless (exists $idhash{$indread});
			my $firstmate=0;
			my $secondmate=0;
			if (exists $idhash{$indread}{$first}) {
				if (exists $idhash{$indread}{$first}{1}) {
					print FILTERED "$indread\t1\t$first\t", $idhash{$indread}{$first}{1}{'pos'}, "\t", $idhash{$indread}{$first}{1}{'str'}, "\n";
					push (@firstpos, 1);
					push (@firstpos, $idhash{$indread}{$first}{1}{'pos'});
					push (@firstpos, $idhash{$indread}{$first}{1}{'str'});
					$firstmate++;
				}
				if (exists $idhash{$indread}{$first}{2}) {
					print FILTERED "$indread\t2\t$first\t", $idhash{$indread}{$first}{2}{'pos'}, "\t", $idhash{$indread}{$first}{2}{'str'}, "\n";
					push (@secondpos, 2);
					push (@secondpos, $idhash{$indread}{$first}{2}{'pos'});
					push (@secondpos, $idhash{$indread}{$first}{2}{'str'});
					$secondmate++;
				}
			}
			if (exists $idhash{$indread}{$second}) {
				if (exists $idhash{$indread}{$second}{1}) {
					print FILTERED "$indread\t1\t$second\t", $idhash{$indread}{$second}{1}{'pos'}, "\t", $idhash{$indread}{$second}{1}{'str'}, "\n";
					push (@secondpos, 1);
					push (@secondpos, $idhash{$indread}{$second}{1}{'pos'});
					push (@secondpos, $idhash{$indread}{$second}{1}{'str'});
					$secondmate++;
				}
				if (exists $idhash{$indread}{$second}{2}) {
					print FILTERED "$indread\t2\t$second\t", $idhash{$indread}{$second}{2}{'pos'}, "\t", $idhash{$indread}{$second}{2}{'str'}, "\n";
					push (@firstpos, 2);
					push (@firstpos, $idhash{$indread}{$second}{2}{'pos'});
					push (@firstpos, $idhash{$indread}{$second}{2}{'str'});
					$firstmate++;
				}
			}
			else {
				die "Error: Non-existed2 readname:contig $indread\n";
			}
#			print "Read: $indread\t$first vs $second\n";
#			print "Test: \@firstpos: ", join ("\t", @firstpos), "\n"; ### For test ###
#			print "Test: \@secondpos: ", join ("\t", @secondpos), "\n"; ### For test ###
			if ($firstmate==2 or $secondmate==2) {
				if ($firstmate==2 and $secondmate==2) {
					print OVERLAP "Info: defined overlap: $first vs $second\tRead: $indread\n";
					
					$finalarray{$arrstart}{'readid'}=$indread;
					@{$finalarray{$arrstart++}{'arr'}}=@firstpos;
					if ($firstpos[2] eq '+') {
						$finalplus{$plusstart}{'readid'}=$indread;
						@{$finalplus{$plusstart++}{'arr'}}=@firstpos;
					}
					elsif ($firstpos[2] eq '-') {
						$finalminus{$minusstart}{'readid'}=$indread;
						@{$finalminus{$minusstart++}{'arr'}}=@firstpos;
					}
					if ($firstpos[2] eq $firstpos[5]) {
						$strand{'s'}++;
					}
					else {
						$strand{'d'}++;
					}
					$finalarray{$arrstart}{'readid'}=$indread;
					@{$finalarray{$arrstart++}{'arr'}}=@secondpos;
					if ($secondpos[2] eq '+') {
						$finalplus{$plusstart}{'readid'}=$indread;
						@{$finalplus{$plusstart++}{'arr'}}=@secondpos;
					}
					elsif ($secondpos[2] eq '-') {
						$finalminus{$minusstart}{'readid'}=$indread;
						@{$finalminus{$minusstart++}{'arr'}}=@secondpos;
					}
					
					if ($secondpos[2] eq $secondpos[5]) {
						$strand{'s'}++;
					}
					else {
						$strand{'d'}++;
					}
				}
				elsif ($firstmate==2 and $secondmate<2) {
					print OVERLAP "Info: possible overlap: $first vs $second\tRead: $indread\n" if ($secondmate==1);
					$finalarray{$arrstart}{'readid'}=$indread;
#						print join("\t", @firstpos), "\n"; ### For test ###
					@{$finalarray{$arrstart++}{'arr'}}=@firstpos;
					if ($firstpos[2] eq '+') {
						$finalplus{$plusstart}{'readid'}=$indread;
						@{$finalplus{$plusstart++}{'arr'}}=@firstpos;
					}
					elsif ($firstpos[2] eq '-') {
						$finalminus{$minusstart}{'readid'}=$indread;
						@{$finalminus{$minusstart++}{'arr'}}=@firstpos;
					}
					if ($firstpos[2] eq $firstpos[5]) {
						$strand{'s'}++;
					}
					else {
						$strand{'d'}++;
					}
				}
				elsif ($firstmate<2 and $secondmate==2) {
					print OVERLAP "Info: possible overlap: $first vs $second\tRead: $indread\n" if ($firstmate==1);
#						print join("\t", @secondpos), "\n"; ### For test ###
					$finalarray{$arrstart}{'readid'}=$indread;
					@{$finalarray{$arrstart++}{'arr'}}=@secondpos;
					if ($secondpos[2] eq '+') {
						$finalplus{$plusstart}{'readid'}=$indread;
						@{$finalplus{$plusstart++}{'arr'}}=@secondpos;
					}
					elsif ($secondpos[2] eq '-') {
						$finalminus{$minusstart}{'readid'}=$indread;
						@{$finalminus{$minusstart++}{'arr'}}=@secondpos;
					}
					if ($secondpos[2] eq $secondpos[5]) {
						$strand{'s'}++;
					}
					else {
						$strand{'d'}++;
					}
				}
				else  {
					print "Warnings: can not understand: $first ($firstmate) vs $second ($secondmate)\tRead: $indread\n\n";
				}
				
				
				unless (&EstimateOverlap($indread, $first, $second)) {
					die "Error: SUB(EstimateOverlap) returns error\n";
				}
				
				if (0) {### For test ###
					print "\%read2gap\n";
					print $indread, "\t$first: ", $seqlength{$first}, "\t$second: ", $seqlength{$second}, "\n";
					print Dumper $idhash{$indread}{$first};
					print Dumper $idhash{$indread}{$second};
					print Dumper $read2gap{$indread}; 
					print "\n"; 
				}
				if (exists $read2gap{$indread}and exists $read2gap{$indread}{'dist'}) {
					push (@disance, $read2gap{$indread}{'dist'});
				}
				if (exists $read2gap{$indread}and exists $read2gap{$indread}{'gap'}) {
					push (@gaps, $read2gap{$indread}{'gap'});
				}
				
			}
			else {
				next;
			}
		}
		
		if (0) {### For test ###
			print "\%read2gap\n"; 
			print Dumper \%read2gap; 
			print "\n"; 
			exit 0;
		}
		my $meandist=0;
		my $stdevdist=0;
		$meandist=int(&CalAverage(@disance));
		$stdevdist=int(&CalStdev(@disance));
		my $meangap=0;
		my $stdevgap=0;
		$meangap=int(&CalAverage(@gaps));
		$stdevgap=int(&CalStdev(@gaps));
		if (abs($meandist)<=$maxinsert) {# and abs($meangap)<=$maxinsert) {### Limit Mean dist and gap
			if ($maxmappingtimes==1) {
				if (abs($meangap)<=$maxinsert) {
					$printthislink=1;
				}
				else {
					$printthislink=0;
				}
			}
			else {
				$printthislink=1;
			}
		}
		else {
			$printthislink=0;
		}
		
		if ($strand{'d'}==0 and $strand{'s'}>=$numpairs) {
			$strness='RC';
			$ratio=$strand{'s'};
		}
		elsif ($strand{'s'}==0 and $strand{'d'}>=$numpairs) {
			$strness='ST';
			$ratio=$strand{'d'};
		}
		else {
			print STDERR "Warnings: unknown strand count: Same:",$strand{'s'},"  Diff",$strand{'d'},"\n";
		}

#			print "Test: \%finalarray\n";print Dumper \%finalarray;print "\n"; ### For test ###
###check distance
		my $total_windows=0;
		if (scalar (keys %finalplus)>=$numpairs) {
			unless(&GetDistance($first, $second, \%finalplus)) {
				print STDERR "Error: \%finalplus\n";
				die "\n";
			}
			$total_windows=scalar (keys %finalplus)-$numpairs+1+$total_windows;
		}
		else {
			print DISTANCE "#NNNplus\t$first\t+\t$second\n";
		}
		if (scalar (keys %finalminus)>=$numpairs) {
			unless(&GetDistance($first, $second, \%finalminus)) {
				print STDERR "Error: \%finalminus\n";
				die "\n";
			}
			$total_windows=scalar (keys %finalminus)-$numpairs+1+$total_windows;
		}
		else {
			print DISTANCE "#NNNminus\t$first\t-\t$second\n";
		}
		
		if ($testwindows>0 and $printthislink==1) {
			print SCAFPAIR $first, "\t", $seqlength{$first}, "\t", $second, "\t", $seqlength{$second}, "\t", $strand{'d'}, '/', $strand{'s'}, "\t", $ratio, "\t$strness\t", $testwindows, "\t", $total_windows, "\t", $testwindows/$total_windows, "\t$test_overlap\t$meandist\t$stdevdist\t$meangap\t$stdevgap\n";
#			print "### Test ###\n"; print Dumper $links{$first}{$second}; exit 0;### For test ###
			foreach my $readidv (@readid) {
				print FINALREADLIST $readidv, "\n" unless (exists $printreads{$readidv});
				$printreads{$readidv}++;
			}
			print FINALREFSLIST $first, "\n" unless (exists $printrefs{$first});
			$printrefs{$first}++;
			print FINALREFSLIST $second, "\n" unless (exists $printrefs{$second});
			$printrefs{$second}++;
			
		}
	}
}
close FILTEROUT;
close FILTERED;
close DISTANCE;
close OVERLAP;
close SCAFPAIR;
close FINALREFSLIST;
close FINALREADLIST;

print "Output unique reads: $unique_mapped_read_num\n";
print "Output joinable reads: ", scalar(keys %printreads), "\n";
print "Output references: ", scalar(keys %printrefs), "\n";


### Global: $numpairs
sub GetDistance {
	my ($GDfirst, $GDsecond, $GDhash)=@_;
	
	my $GDsubinfo='SUB(GetDistance)';
	
	my ($GDtest_sort, $GDreadorder)=&SortIdHash($GDhash);
	unless ($GDtest_sort) {
		print STDERR $GDsubinfo, "Error: SortIdHash failed\n";
		return 0;
	}

	for (my $GDi=0; $GDi<=(scalar(@{$GDreadorder})-$numpairs); $GDi++) {
		my @GDtemparrpos1=();
		my @GDtempstrand1=();
		my @GDtemparrpos2=();
		my @GDtempstrand2=();
		my @GDtempfix=();
		my @GDtempflanking=();
		for (my $GDj=0; $GDj<$numpairs; $GDj++) {
			my $GDorder=${$GDreadorder}[$GDi+$GDj];
#			print " $GDi $GDj ", ${$GDhash}{$GDorder}{'readid'}, " $GDfirst vs $GDsecond\n";### For test ###
			unless (exists ${$GDhash}{$GDorder} and exists ${$GDhash}{$GDorder}{'arr'} and defined ${$GDhash}{$GDorder}{'arr'}[1] and ${$GDhash}{$GDorder}{'arr'}[1]=~/^\d+$/ and defined ${$GDhash}{$GDorder}{'arr'}[4] and ${$GDhash}{$GDorder}{'arr'}[4]=~/^\d+$/) {
				print STDERR $GDsubinfo, "Error: second pos: $GDfirst vs $GDsecond: \n";
				return 0;
			}
			my ($GDmatenum1, $GDmatepos1, $GDmatestr1, $GDmatenum2, $GDmatepos2, $GDmatestr2)=@{${$GDhash}{$GDorder}{'arr'}};
			push (@GDtemparrpos1, $GDmatepos1);
			push (@GDtempstrand1, $GDmatestr1);
			push (@GDtemparrpos2, $GDmatepos2);
			push (@GDtempstrand2, $GDmatestr2);
			
			my ($GDtest, $GDfixed, $GDgapsize)= &EstimateGap($GDfirst, $GDsecond, ${$GDhash}{$GDorder});
			unless ($GDtest) {
				print STDERR $GDsubinfo, "Error: EstimateGap failed: $GDfirst vs $GDsecond\n";
				print Dumper ${$GDhash}{$GDorder};
				return 0;
			}
			
#			push (@GDtempfix, $GDfixed);
			push (@GDtempflanking, $GDgapsize);
			
#			print " $GDi $GDj ", ${$GDhash}{$GDorder}{'readid'}, " $GDfirst vs $GDsecond $GDgapsize\n"; ### For test ###
		}
		my @GDnewtemparrpos2=sort {$a<=>$b} @GDtemparrpos2;
		my $GDmaxdistance1=$GDtemparrpos1[-1]-$GDtemparrpos1[0];
		my $GDmaxdistance2=$GDnewtemparrpos2[-1]-$GDnewtemparrpos2[0];
		
#		$test_overlap=0;
#		if (sum(@GDtempfix) >1) {$test_overlap=1;} ### would defined overlap if 1 of then show overlap
		my $GDaveragedistancebetween=&CalAverage(@GDtempflanking);
#		print DISTANCE"No/Total-NumReads\tRef1name\tRef1len\tRef1std\tRef2name\tRef2len\tRef2std\tS/D\tratio\tAssemblyStd\tDist1\tDist1Detail\tDist2\tDist2Detail\tTest_overlap\tMeaninternalSize\tDetail\n";
		print DISTANCE $GDi,'/',scalar(keys %{$GDhash}),"-$numpairs\t", "$GDfirst\t", $seqlength{$GDfirst}, "\t", join ('', @GDtempstrand1),"\t", $GDsecond, "\t", $seqlength{$GDsecond}, "\t", join ('', @GDtempstrand2), "\t", $strand{'s'}, '/', $strand{'d'}, "\t$ratio\t$strness\t", "$GDmaxdistance1\t(", join (",", @GDtemparrpos1), ")\t$GDmaxdistance2\t(", join (',', @GDtemparrpos2), ")\t$test_overlap\t", $GDaveragedistancebetween, "\t(", join(',', @GDtempflanking), ")\n";
		
		if ($windows2!=0 and $GDmaxdistance1<=$windows1 and $GDmaxdistance2<=$windows2 and $GDaveragedistancebetween<=$maxinsert) {
			$testwindows++;
		}
	}
	
	
	return 1;
}




### Global: %idhash
sub EstimateGap {
	my ($EGref1, $EGref2, $EGhashref)=@_;
	
	my $EGsubinfo='SUB(EstimateGap)';
	my @EGfixedgap=();
	my $EGestgapsize=0;
	my $EGfinalsize=0;
	my $EGfinalfix=0;
	my $EGflank1=0;
	my $EGflank2=0;
	
	my $EGreadid=${$EGhashref}{'readid'};
	my ($EGmatenum1, $EGmatepos1, $EGmatestr1, $EGmatenum2, $EGmatepos2, $EGmatestr2)=@{${$EGhashref}{'arr'}};
	unless (defined $EGmatenum1 and $EGmatenum1=~/^[12]{1}$/ and defined $EGmatenum2 and $EGmatenum2=~/^[12]{1}$/ and ($EGmatenum1 != $EGmatenum2)) {
		print STDERR $EGsubinfo, "Error: read $EGreadid array: $EGmatenum1, $EGmatepos1, $EGmatestr1, $EGmatenum2, $EGmatepos2, $EGmatestr2\n";
		return 0;
	}
	unless (exists $idhash{$EGreadid} and exists $idhash{$EGreadid}{$EGref1}) {
		print STDERR $EGsubinfo, "Error: idhash not existed for read $EGreadid\n";
		return 0;
	}
	unless (exists $idhash{$EGreadid}{$EGref1}{$EGmatenum1}) {
		print STDERR $EGsubinfo, "Error: READ $EGreadid R$EGmatenum1 REF $EGref1 ARRAY $EGmatenum1, $EGmatepos1, $EGmatestr1\n";
		return 0;
	}
	unless (exists $idhash{$EGreadid}{$EGref2}{$EGmatenum2}) {
		print STDERR $EGsubinfo, "Error: READ $EGreadid R$EGmatenum2 REF $EGref1 ARRAY $EGmatenum2, $EGmatepos2, $EGmatestr2\n";
		return 0;
	}
	if (exists $idhash{$EGreadid}{$EGref1}{$EGmatenum2}) {
		my $EGgap1=abs($idhash{$EGreadid}{$EGref1}{$EGmatenum2}{'pos'}-$EGmatepos1)+1;
		push (@EGfixedgap, $EGgap1);
	}
	if (exists $idhash{$EGreadid}{$EGref2}{$EGmatenum1}) {
		my $EGgap2=abs($idhash{$EGreadid}{$EGref2}{$EGmatenum1}{'pos'}-$EGmatepos2)+1;
		push (@EGfixedgap, $EGgap2);
	}
	
	if (scalar(@EGfixedgap)==2) {
		$EGfinalsize= max @EGfixedgap;
		$EGfinalfix=1;
	}
	elsif (scalar(@EGfixedgap)==1) {
		$EGfinalsize=$EGfixedgap[0];
		$EGfinalfix=1;
	}
	else {
		unless (exists $seqlength{$EGref1} and $seqlength{$EGref1}>=$EGmatepos1) {
			print STDERR $EGsubinfo, "Error: unknown length: $EGref1\n";
			return 0;
		}
		unless (exists $seqlength{$EGref2} and $seqlength{$EGref2}>=$EGmatepos2) {
			print STDERR $EGsubinfo, "Error: unknown length: $EGref2\n";
			return 0;
		}
		$EGflank1=0;
		$EGflank2=0;
		if ($EGmatestr1 eq '+') {
			$EGflank1=$seqlength{$EGref1}-$EGmatepos1+1;
		}
		elsif ($EGmatestr1 eq '-') {
			$EGflank1=$EGmatepos1;
		}
		else {
			print STDERR $EGsubinfo, "Error: unknown R$EGmatenum1 strand: READID $EGreadid REF $EGref1 STRAND $EGmatestr1\n";
			return 0;
		}
		
		if ($EGmatestr2 eq '+') {
			$EGflank2=$seqlength{$EGref2}-$EGmatepos2+1;
		}
		elsif ($EGmatestr2 eq '-') {
			$EGflank2=$EGmatepos2;
		}
		else {
			print STDERR $EGsubinfo, "Error: unknown R$EGmatenum2 strand: READID $EGreadid REF $EGref2 STRAND $EGmatestr2\n";
			return 0;
		}
		
		if (defined $EGflank1 and $EGflank1=~/^\d+$/ and defined $EGflank2 and $EGflank2=~/^\d+$/) {
			$EGfinalsize=$EGflank1+$EGflank2;
		}
		else {
			print STDERR $EGsubinfo, "Error: unknown R$EGmatenum1 flanking size: READID $EGreadid REF $EGref1 STRAND $EGmatestr1\n" unless (defined $EGflank1 and $EGflank1=~/^\d+$/);
			print STDERR $EGsubinfo, "Error: unknown R$EGmatenum2 flanking size: READID $EGreadid REF $EGref2 STRAND $EGmatestr2\n" unless (defined $EGflank2 and $EGflank2=~/^\d+$/);
			return 0;
		}
	}
	
#	print "READID $EGreadid REF $EGref1 ($seqlength{$EGref1}) POS1 $EGmatepos1 STRAND $EGmatestr1 REF $EGref2($seqlength{$EGref2}) POS2 $EGmatepos2 STRAND $EGmatestr2 Fix $EGfinalfix FLANK1 $EGflank1 FLANK2 $EGflank2\n"; ### for test ###
	
	return (1, $EGfinalfix, $EGfinalsize);
}




sub CalAverage {
	return sum(@_)/@_;
}


sub CalStdev {
	my $CSsum=0;
	my $CSaverage=&CalAverage(@_);
	
	foreach (@_) {
		$CSsum=$CSsum+($_-$CSaverage) ** 2;
	}
	$CSsum = $CSsum / scalar(@_);
	$CSsum= sqrt ($CSsum);
	
	return $CSsum;
}

### Global: $numpairs
### Dependency: 
### Note: 
sub SortIdHash {
	my $SIHhash=shift;
	
	my @SIHarr=();
	my $SIHsubinfo='SUB(SortIdHash)';
	my %SIHtemp=();
	
	if (scalar(keys %{$SIHhash})<$numpairs) {
		print STDERR $SIHsubinfo, "Error: invalid input hash: ", scalar(keys %{$SIHhash}), "<5 keys\n";
		return 0;
	}
	
	foreach my $SIHnum (keys %{$SIHhash}) {
		if (exists ${$SIHhash}{$SIHnum} and exists ${$SIHhash}{$SIHnum}{'arr'}) {
			if (scalar(@{${$SIHhash}{$SIHnum}{'arr'}})==6) {
				$SIHtemp{${$SIHhash}{$SIHnum}{'arr'}[1]}{${$SIHhash}{$SIHnum}{'arr'}[4]}{$SIHnum}++;
			}
			else {
				print STDERR $SIHsubinfo, 'Warnings: invalid array element number: Read ', ${$SIHhash}{$SIHnum}{'readid'},"\t", join (" ", @{${$SIHhash}{$SIHnum}{'arr'}}), "\n";
				return 0;
			}
		}
		else {
			print STDERR $SIHsubinfo, "Warnings: none existing array($SIHnum)\n";
			return 0;
		}
	}
	
	foreach my $SIHnum (sort {$a<=>$b} keys %SIHtemp) {
		foreach my $SIHnum2 (sort {$a<=>$b} keys %{$SIHtemp{$SIHnum}}) {
			my @SIGarr2=();
			@SIGarr2=keys %{$SIHtemp{$SIHnum}{$SIHnum2}};
			foreach my $SIHnum3 (@SIGarr2) {
				push (@SIHarr, $SIHnum3);
			}
		}
	}
	
	unless (scalar(@SIHarr)>=$numpairs and scalar(@SIHarr)==scalar(keys %{$SIHhash})) {
		print STDERR $SIHsubinfo, "Error: invalid output array: <5 keys\n";
		return 0;
	}
	
	%SIHtemp=();
	
	return (1, \@SIHarr);
}




### Global: %read2links, %links, %idhash
sub RemoveRead {
	my $RRreadname=shift;
	
	if (exists $idhash{$RRreadname}) {
		delete $idhash{$RRreadname};
	}
	if (exists $read2links{$RRreadname}) {
		my @RRarr1=keys %{$read2links{$RRreadname}};
		foreach my $RRele1 (@RRarr1) {
			my @RRarr2=keys %{$read2links{$RRreadname}{$RRele1}};
			foreach my $RRele2 (@RRarr2) {
				if (exists $links{$RRele1} and exists $links{$RRele1}{$RRele2} and exists $links{$RRele1}{$RRele2}{$RRreadname}) {
					delete $links{$RRele1}{$RRele2}{$RRreadname};
					delete $links{$RRele1}{$RRele2} if (scalar(keys %{$links{$RRele1}{$RRele2}})>0);
					delete $links{$RRele1} if (scalar(keys %{$links{$RRele1}})>0);
				}
				if (exists $links{$RRele2} and exists $links{$RRele2}{$RRele1} and exists $links{$RRele2}{$RRele1}{$RRreadname}) {
					delete $links{$RRele2}{$RRele1}{$RRreadname};
					delete $links{$RRele2}{$RRele1} if (scalar(keys %{$links{$RRele2}{$RRele1}})>0);
					delete $links{$RRele2} if (scalar(keys %{$links{$RRele2}})>0);
				}
			}
		}
	}
	
	return 1;
}





### Global: %idhash, %seqlength, $estimatedinsertsize, %read2gap
sub EstimateOverlap {
	my ($EOreadid, $EOfirst, $EOsecond)=@_;
	
	my $EOsubinfo='SUB(EstimateOverlap)';
	my @EOarr1=(); my @EOarr2=();
	my $EOnum1=0;my $EOnum2=0;
	
	unless (exists $idhash{$EOreadid} and exists $idhash{$EOreadid}{$EOfirst} and exists $idhash{$EOreadid}{$EOsecond}) {
		print STDERR $EOsubinfo, "Error: read hash exists: $EOreadid\n";
		return 0;
	}
	unless (exists $seqlength{$EOfirst}) {
		print STDERR $EOsubinfo, "Error: seq length NOT exists: $EOfirst\n";
		return 0;
	}
	unless (exists $seqlength{$EOsecond}) {
		print STDERR $EOsubinfo, "Error: seq length NOT exists: $EOsecond\n";
		return 0;
	}
	
	if (exists $idhash{$EOreadid}{$EOfirst}{1} and exists $idhash{$EOreadid}{$EOsecond}{2}) {
		$EOnum1=1; $EOnum2=2;
	}
	elsif (exists $idhash{$EOreadid}{$EOfirst}{2} and exists $idhash{$EOreadid}{$EOsecond}{1}) {
		$EOnum1=2; $EOnum2=1;
	}
	else {
		print STDERR $EOsubinfo, "Error: read not paired: $EOreadid\n";
		return 0;
	}
	
	my $EOpairednum1=0;
	if (exists $idhash{$EOreadid}{$EOfirst}{$EOnum2}) {
		$EOpairednum1=1;
	}
	my $EOpairednum2=0;
	if (exists $idhash{$EOreadid}{$EOsecond}{$EOnum1}) {
		$EOpairednum2=1;
	}
	
#	print "Info: READ $EOreadid First $EOfirst Mate $EOnum1 Paired? $EOpairednum1 SECOND $EOsecond Mate $EOnum2 Paired? $EOpairednum2\n"; ### For test ###
	if ($EOpairednum1==0 and $EOpairednum2==0) {
		my $EOd1=0; my $EOd2=0;
		if ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '+') {
			$EOd1=$seqlength{$EOfirst}-$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '-') {
			$EOd1=$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid strand: READ $EOreadid MATE $EOnum1 REFERENCE $EOfirst\n";
			return 0;
		}
		if ($idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '+') {
			$EOd2=$seqlength{$EOsecond}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '-') {
			$EOd2=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid strand: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			return 0;
		}
		$read2gap{$EOreadid}{'dist'}=$EOd1+$EOd2;
		$read2gap{$EOreadid}{'gap'}=$estimatedinsertsize-($EOd1+$EOd2);
	}
	elsif ($EOpairednum1==1 and $EOpairednum2==0) {
		$test_overlap=1;
		my $EOd1=0; my $EOd2=0;
		if ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '+' and $idhash{$EOreadid}{$EOfirst}{$EOnum2}{'str'} eq '-') {
			$read2gap{$EOreadid}{'dist'}=$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}-$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}+1;
			$EOd1=$seqlength{$EOfirst}-$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '-' and $idhash{$EOreadid}{$EOfirst}{$EOnum2}{'str'} eq '+') {
			$read2gap{$EOreadid}{'dist'}=$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}-$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}+1;
			$EOd1=$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read mapping1: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			print Dumper $idhash{$EOreadid}{$EOfirst};
			print Dumper $idhash{$EOreadid}{$EOsecond};
			return 0;
		}
		if ($idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '+') {
			$EOd2=$seqlength{$EOsecond}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '-') {
			$EOd2=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read strand1: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			return 0;
		}
		if (($EOd2> (max $read2gap{$EOreadid}{'dist'}, $maxinsert))  or ($read2gap{$EOreadid}{'dist'}<=0) ) {
			$readproblematic=1;
			return 1;
		}
		$read2gap{$EOreadid}{'gap'}=-($EOd1+$EOd2-1);
	}
	elsif ($EOpairednum1==0 and $EOpairednum2==1) {
		$test_overlap=1;
		my $EOd1=0; my $EOd2=0;
		if ($idhash{$EOreadid}{$EOsecond}{$EOnum1}{'str'} eq '+' and $idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '-') {
			$read2gap{$EOreadid}{'dist'}=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}-$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'}+1;
			$EOd1=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'};
		}
		elsif ($idhash{$EOreadid}{$EOsecond}{$EOnum1}{'str'} eq '-' and $idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '+') {
			$read2gap{$EOreadid}{'dist'}=$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
			$EOd1=$seqlength{$EOsecond}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read mapping2: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			print Dumper $idhash{$EOreadid}{$EOsecond};
			print Dumper $idhash{$EOreadid}{$EOsecond};
			return 0;
		}
		if ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '+') {
			$EOd2=$seqlength{$EOfirst}-$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '-') {
			$EOd2=$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read strand2: READ $EOreadid MATE $EOnum1 REFERENCE $EOfirst\n";
			return 0;
		}
		if (($EOd2> (max $read2gap{$EOreadid}{'dist'}, $maxinsert)) or ($read2gap{$EOreadid}{'dist'}<=0) ) {
			$readproblematic=1;
			return 1;
		}
		$read2gap{$EOreadid}{'gap'}=-($EOd1+$EOd2-1);
	}
	elsif ($EOpairednum1==1 and $EOpairednum2==1) {
		$test_overlap=1;
		my $EOd1=0; my $EOd2=0; my $EOd3=0; my $EOd4=0; my $EOinsert1=0; my $EOinsert2=0;
		if ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '+' and $idhash{$EOreadid}{$EOfirst}{$EOnum2}{'str'} eq '-') {
			$EOinsert1=$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}-$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}+1;
			$EOd1=$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'};
			$EOd2=$seqlength{$EOfirst}-$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOfirst}{$EOnum1}{'str'} eq '-' and $idhash{$EOreadid}{$EOfirst}{$EOnum2}{'str'} eq '+') {
			$EOinsert1=$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}-$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'}+1;
			$EOd1=$seqlength{$EOfirst}-$idhash{$EOreadid}{$EOfirst}{$EOnum1}{'pos'}+1;
			$EOd2=$idhash{$EOreadid}{$EOfirst}{$EOnum2}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read mapping3: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			print Dumper $idhash{$EOreadid}{$EOfirst};
			print Dumper $idhash{$EOreadid}{$EOsecond};
			return 0;
		}
		
		if ($idhash{$EOreadid}{$EOsecond}{$EOnum1}{'str'} eq '+' and $idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '-') {
			$EOinsert2=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}-$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'}+1;
			$EOd3=$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'};
			$EOd4=$seqlength{$EOsecond}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
		}
		elsif ($idhash{$EOreadid}{$EOsecond}{$EOnum1}{'str'} eq '-' and $idhash{$EOreadid}{$EOsecond}{$EOnum2}{'str'} eq '+') {
			$read2gap{$EOreadid}{'dist'}=$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'}-$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'}+1;
			$EOd3=$seqlength{$EOsecond}-$idhash{$EOreadid}{$EOsecond}{$EOnum1}{'pos'}+1;
			$EOd4=$idhash{$EOreadid}{$EOsecond}{$EOnum2}{'pos'};
		}
		else {
			print STDERR $EOsubinfo, "Error: invalid read mapping4: READ $EOreadid MATE $EOnum2 REFERENCE $EOsecond\n";
			print Dumper $idhash{$EOreadid}{$EOsecond};
			print Dumper $idhash{$EOreadid}{$EOsecond};
			return 0;
		}
		
		$read2gap{$EOreadid}{'dist'}= max($EOinsert1, $EOinsert2);

		$read2gap{$EOreadid}{'gap'}=-($read2gap{$EOreadid}{'dist'} + (min($EOd1, $EOd3)) + min ($EOd2, $EOd4) -2);
	}
	
	return 1;
}


__END__
