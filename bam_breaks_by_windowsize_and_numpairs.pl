#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/sum min max/;
use Data::Dumper; ### For test ###
use constant USAGE =><<EOH;

usage: $0 in.stats fasta.fai insert_mean insert_stdev Num_Pairs out.breaks

This script is used to detect breaks problem

  in.stats: bam_scaffolding_separately.stats_noBioDBSAM.pl output
  fasta.fai: index for fasta
  insert_mean: Mean insert size
  insert_stdev: StDev insert size
  Num_Pairs: INT pairs would be count as 1 evidence
  out.breaks: output

v20170417

EOH

die USAGE if (scalar(@ARGV) !=6 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $input=shift @ARGV;
my $fastaindex=shift @ARGV;
my $distave=shift @ARGV;
my $disstdev=shift @ARGV;
my $numpairs=shift @ARGV;
my $outfile_breaks=shift @ARGV;


### Configurable
my $maxmappingtimes=1;
my $opt_disance=10000;
my $debug=0;

### Default
my $linenum=0;
my %idhash=();
my %allrefs=();
my %seqlength=();
my $numContigs=0;
my $total_bases=0;
my $total_errors=0;



#input output
unless (defined $input and -s $input) {
	die "Error: invalid input.stats\n";
}
unless (defined $numpairs and $numpairs=~/^\d+/) {
	die "Error: invalid Num_Pairs, should be INT\n";
}
unless (defined $outfile_breaks and $outfile_breaks=~/^[-\.\/\w_]+$/) {
	die "Error: invalid out.breaks\n";
}
unlink $outfile_breaks if (-e $outfile_breaks);



$linenum=0;
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
print "    STATS:                      $input\n";
print "    Index:                      $fastaindex\n";
print "    Read total fasta.fai lines: $linenum\n";
print "    Total seqeunces:            ", scalar(keys %seqlength), "\n\n";



$linenum=0;
open (STATSIN, " < $input") || die "Error: can not open in.stats: $input\n";
while (my $line=<STATSIN>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	if (exists $idhash{$arr[0]}{$arr[2]}{$arr[1]}) {
		$idhash{$arr[0]}{'EXCLUDEDPLZ'}++;
		next;
	}
	if ($arr[5] eq '+') {
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'pos'}=$arr[3];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'end'}=$arr[4];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'str'}=$arr[5];
	}
	elsif ($arr[5] eq '-') {
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'pos'}=$arr[4];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'end'}=$arr[3];
		$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'str'}=$arr[5];
	}
	else {
		die "Error: strand error at line($linenum): $line\n";
	}
}
close STATSIN;



my @allreads=();
@allreads=keys (%idhash);
foreach my $indread (@allreads) {
	if (exists $idhash{$indread}{'EXCLUDEDPLZ'}) {
		delete $idhash{$indread};
		next;
	}
	
	my @contigs=keys %{$idhash{$indread}};
	my @arr2=();
	my @arr3=();
	foreach my $contig (@contigs) {
		push (@arr2, $contig) if (exists $idhash{$indread}{$contig}{1});
		push (@arr3, $contig) if (exists $idhash{$indread}{$contig}{2});
	}
	my $maxmapping= max (scalar(@arr2), scalar(@arr3));
	if ($maxmapping>$maxmappingtimes) {### control mapping times
		delete $idhash{$indread};
		next;
	}
	
	foreach (@contigs) {
		$allrefs{$_}{$indread}++;
	}
	
}


my %goodregions=();
my %forwardstrand=();
my %reversestrand=();
my %badregions=();
my $total_good_length=0;
my $total_good_num=0;
my $contigs;
open (BREAKOUT, " > $outfile_breaks") || die "Error: can not write $outfile_breaks\n";
foreach $contigs (sort keys %allrefs) {
	my @readnames=keys %{$allrefs{$contigs}};
	unless (scalar(@readnames)>=$numpairs) {
		print STDERR "Warnings: less than $numpairs pairs on Seq $contigs\n" if ($debug);
		next;
	}
	my $numunpaired=0;
	my $numpaired=0;
	my %readinfo=();
	my $numproperpaired=0;
	%goodregions=();
	%forwardstrand=();
	%reversestrand=();
	%badregions=();
	
	### detect paired or unpaired
	foreach my $readid (@readnames) {
		if (exists $idhash{$readid}{$contigs}{1} and $idhash{$readid}{$contigs}{2}) {
			$numpaired++;
			if (($idhash{$readid}{$contigs}{1}{'str'} eq '+') and ($idhash{$readid}{$contigs}{2}{'str'} eq '-')) {
				my $distance=$idhash{$readid}{$contigs}{2}{'pos'}-$idhash{$readid}{$contigs}{1}{'pos'};
				if ($distance>=($distave-$disstdev) and $distance<=($distave+$disstdev)) {
					if (exists $goodregions{$idhash{$readid}{$contigs}{1}{'pos'}}) {
						if ($idhash{$readid}{$contigs}{2}{'pos'}>$goodregions{$idhash{$readid}{$contigs}{1}{'pos'}}) {
							$goodregions{$idhash{$readid}{$contigs}{1}{'pos'}}=$idhash{$readid}{$contigs}{2}{'pos'};
						}
					}
					else {
						$goodregions{$idhash{$readid}{$contigs}{1}{'pos'}}=$idhash{$readid}{$contigs}{2}{'pos'};
					}
					$numproperpaired++;
					$readinfo{$readid}{'paired'}++;
				}
				
			}
			elsif (($idhash{$readid}{$contigs}{1}{'str'} eq '-') and ($idhash{$readid}{$contigs}{2}{'str'} eq '+')) {
				my $distance=$idhash{$readid}{$contigs}{1}{'pos'}-$idhash{$readid}{$contigs}{2}{'pos'};
				if ($distance>=($distave-$disstdev) and $distance<=($distave+$disstdev)) {
					if (exists $goodregions{$idhash{$readid}{$contigs}{2}{'pos'}}) {
						if ($idhash{$readid}{$contigs}{1}{'pos'}>$goodregions{$idhash{$readid}{$contigs}{2}{'pos'}}) {
							$goodregions{$idhash{$readid}{$contigs}{2}{'pos'}}=$idhash{$readid}{$contigs}{1}{'pos'};
						}
					}
					else {
						$goodregions{$idhash{$readid}{$contigs}{2}{'pos'}}=$idhash{$readid}{$contigs}{1}{'pos'};
					}
					$numproperpaired++;
					$readinfo{$readid}{'paired'}++;
				}
			}
		}
		else {
			$numunpaired++;
			$readinfo{$readid}{'unpaired'}++;
		}
	}

	
	
	my $temphash1=&MergeRegion(\%goodregions);
#	print "Test: $contigs Total ",scalar(@readnames) ," Unpaired $numunpaired Paired $numpaired Properpaired $numproperpaired\n";
#	print Dumper $temphash1; ### For test ###
	if (scalar(keys %{$temphash1})>0) {
		$total_good_num++;
		foreach my $goodstart (sort {$a<=>$b} keys %{$temphash1}) {
			$total_good_length=$total_good_length+(${$temphash1}{$goodstart}-$goodstart+1);
		}
	}
	
	unless ($numunpaired>=$numpairs) {
		print STDERR "Warnings: less than $numpairs unpaired on Seq $contigs\n" if ($debug);
		next;
	}
	
	### group
	my $temp_i=0;
	foreach my $readid (@readnames) {
		my @matenum=keys %{$idhash{$readid}{$contigs}};
		foreach my $whichmate (@matenum) {
			if ($idhash{$readid}{$contigs}{$whichmate}{'str'} eq '+') {
				$forwardstrand{$temp_i++}=[$readid, $whichmate, $idhash{$readid}{$contigs}{$whichmate}{'pos'}];
			}
			else {
				$reversestrand{$temp_i++}=[$readid, $whichmate, $idhash{$readid}{$contigs}{$whichmate}{'pos'}];
			}
		}
	}
	my $forwref1=&SortHash(\%forwardstrand);
	my $revsref2=&SortHash(\%reversestrand);
	
	for (my $temp_x=0; $temp_x<(scalar(@{$forwref1})-$numpairs); $temp_x++) {
		my @temparr=();
		foreach (my $temp_y=$temp_x; $temp_y<($temp_x+$numpairs); $temp_y++) {
			push (@temparr, ${$forwref1}[$temp_y]);
		}
		next unless (($forwardstrand{$temparr[-1]}[2]-$forwardstrand{$temparr[0]}[2])<=$opt_disance);
		my $total_evidence=0;
		foreach my $temp2 (@temparr) {
			my $temp3=$forwardstrand{$temp2}[0];
			if (exists $readinfo{$temp3} and exists $readinfo{$temp3}{'unpaired'}) {
				$total_evidence++;
			}
		}
		if ($total_evidence==$numpairs) {
			my $temp4=$forwardstrand{$temparr[-1]}[2]+($distave+$disstdev);
			if ($temp4 < ($seqlength{$contigs}-$disstdev)) {
				my $temp5=$idhash{$forwardstrand{$temparr[-1]}[0]}{$contigs}{$forwardstrand{$temparr[-1]}[1]}{'end'}+1;
				unless (exists $badregions{$temp5} and $temp4<=$badregions{$temp5}) {
					$badregions{$temp5}=$temp4;
				}
			}
		}
	}
#	print "Test: Bad region1: \%badregions\n"; print Dumper \%badregions; print "\n"; ### For test ###
	for (my $temp_x=0; $temp_x<(scalar(@{$revsref2})-$numpairs); $temp_x++) {
		my @temparr=();
		foreach (my $temp_y=$temp_x; $temp_y<($temp_x+$numpairs); $temp_y++) {
			push (@temparr, ${$revsref2}[$temp_y]);
		}
		next unless (($reversestrand{$temparr[-1]}[2]-$reversestrand{$temparr[0]}[2])<=$opt_disance);
		my $total_evidence=0;
		foreach my $temp2 (@temparr) {
			my $temp3=$reversestrand{$temp2}[0];
			if (exists $readinfo{$temp3} and exists $readinfo{$temp3}{'unpaired'}) {
				$total_evidence++;
			}
		}
		if ($total_evidence==$numpairs) {
			my $temp4=$reversestrand{$temparr[0]}[2] - ($distave+$disstdev);
			if ($temp4 > $disstdev) {
				my $temp5=$idhash{$reversestrand{$temparr[0]}[0]}{$contigs}{$reversestrand{$temparr[0]}[1]}{'end'}-1;
				unless (exists $badregions{$temp4} and $temp5<=$badregions{$temp4}) {
					$badregions{$temp4}=$temp5;
				}
			}
		}
	}
	next unless (scalar(keys %badregions)>0);
#	print "Test: Bad region2: \%badregions\n"; print Dumper \%badregions; print "\n"; ### For test ###
	
	my $temphash2={};
	$temphash2=&MergeRegion(\%badregions);
	
#	print "Test: Good region: \$temphash1\n"; print Dumper $temphash1; print "\n"; ### For test ###
#	print "Test: Bad region: \$temphash2\n"; print Dumper $temphash2; print "\n"; ### For test ###
	
	my $finalhash={};
	$finalhash=&BadDeductGood($temphash2, $temphash1);
#	print "Test: $contigs length $seqlength{$contigs} Problematic region: \$finalhash\n"; print Dumper $finalhash; print "\n"; ### For test ###
	next unless (scalar(keys %{$finalhash})>0);
	$numContigs++;

	foreach (sort {$a<=>$b} keys %{$finalhash}) {
		print BREAKOUT $contigs, ':', $_, '-', ${$finalhash}{$_}, "\n";
		$total_bases=$total_bases+${$finalhash}{$_}-$_+1;
		$total_errors++;
	}
}
close BREAKOUT;



print "##### SUMMARY 2 #####\n";
print "    Good:\n";
print "        Total number:      $total_good_num\n";
print "        Total bases:       $total_good_length\n";
print "    Error:\n";
print "        Num errors:        $total_errors\n";
print "        Num contigs:       $numContigs\n";
print "        Num bases:         $total_bases\n";



###
sub MergeRegion {
	my $MRorginalhash=shift;

	my %MRhash=();
	my $MRtest=0;
	my $laststart=0;
	my $lastend=0;
	foreach my $MRnum1 (sort {$a<=>$b} keys %{$MRorginalhash}) {
		if ($MRtest==0) {
			$laststart=$MRnum1;
			$lastend=${$MRorginalhash}{$MRnum1};
			$MRtest++;
		}
		else {
			if ($MRnum1<=($lastend+1)) {
				$lastend=${$MRorginalhash}{$MRnum1};
			}
			elsif ($MRnum1>($lastend+1)) {
				$MRhash{$laststart}=$lastend;
				$laststart=$MRnum1;
				$lastend=${$MRorginalhash}{$MRnum1};
			}
		}
	}
	$MRhash{$laststart}=$lastend;
	
	return \%MRhash;
}


### Global: $contigs
sub SortHash {
	my $SHhash=shift;
	
#	print "Test: Sort Hash\n"; print Dumper $SHhash; print "\n"; ### For test ###
	
	my @SHarr=();
	my %MRtemphash=();
	
	foreach my $SHarrref (keys %{$SHhash}) {
		$MRtemphash{${$SHhash}{$SHarrref}[2]}{$SHarrref}++;
	}
	foreach (sort {$a<=>$b} keys %MRtemphash) {
		###Control repeats; TO BE CONTINUED
#		my @SHarr2=sort {$a<=>$b} keys %{$MRtemphash{$_}};
#		
#		$contigs
		foreach my $MRnum (keys %{$MRtemphash{$_}}) {
			push (@SHarr, $MRnum);
		}
	}
	
#	print "Test: Order \n"; print Dumper \@SHarr; print "\n"; ### For test ###
	return \@SHarr;
}



sub BadDeductGood {
	my ($BDGhash1, $BDGhash2)=@_;
	
	my $BDGtest=0;
	my $BDGsubinfo='SUB(BadDeductGood)';
	
	BDGLOOP: while ($BDGtest==0) {
		my $BDGoverlapnum=0;
#		print $BDGsubinfo, "Test: \$BDGhash1\n"; print Dumper $BDGhash1; print "\n";
		my @BDGarr1=sort {$a<=>$b} keys %{$BDGhash1};
		foreach my $BDGbadstart (@BDGarr1) {
			my $BDGbadend=${$BDGhash1}{$BDGbadstart};
			foreach my $BDGgoodstart (sort {$a<=>$b} keys %{$BDGhash2}) {
				my $BDGgoodend=${$BDGhash2}{$BDGgoodstart};
				if ($BDGgoodend<$BDGbadstart or $BDGgoodstart>$BDGbadend) {
					next;
				}
				else {
					if ($BDGgoodstart<=$BDGbadstart and $BDGgoodend>=$BDGbadend) {
						delete ${$BDGhash1}{$BDGbadstart}
					}
					elsif ($BDGgoodend<=$BDGbadend and $BDGgoodend>=$BDGbadstart) {
						if ($BDGgoodstart>$BDGbadstart) {
							${$BDGhash1}{$BDGbadstart}=$BDGgoodstart-1;
						}
						else {
							delete ${$BDGhash1}{$BDGbadstart};
						}
#						if ($BDGgoodend<$BDGbadend) {
#							if (exists ${$BDGhash1}{$BDGgoodend+1}) {
#								die "($BDGsubinfo)Error: 111\n";
#							}
#							${$BDGhash1}{$BDGgoodend+1}=$BDGbadend;
#						}
						
					}
					elsif ($BDGgoodstart>=$BDGbadstart and $BDGgoodstart<=$BDGbadend) {
						if ($BDGgoodstart>$BDGbadstart) {
							${$BDGhash1}{$BDGbadstart}=$BDGgoodstart-1;
						}
						else {
							delete ${$BDGhash1}{$BDGbadstart};
						}
#						if ($BDGgoodend<$BDGbadend) {
#							if (exists ${$BDGhash1}{$BDGgoodend+1}) {
#								die "($BDGsubinfo)Error: 222\n";
#							}
#							${$BDGhash1}{$BDGgoodend+1}=$BDGbadend;
#						}
					}
					else {
						die "($BDGsubinfo)Error: unknow range Good $BDGgoodstart-$BDGgoodend Bad $BDGbadstart-$BDGbadend\n";
					}
					$BDGoverlapnum++;
					next BDGLOOP;
				}
			}
		}
		if ($BDGoverlapnum==0) {
			$BDGtest=1;
		}
		else {
			$BDGtest=0;
		}
	}
	return $BDGhash1;
}
