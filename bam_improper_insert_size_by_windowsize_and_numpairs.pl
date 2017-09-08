#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/sum min max/;
use Data::Dumper; ### For test ###
use constant USAGE =><<EOH;

usage: $0 in.stats insert_min insert_max Num_Pairs out.improper

This script is used to detect breaks and orientation problem

  in.stats: bam_scaffolding_separately.stats_noBioDBSAM.pl output
  insert_min: minimum insert size: mean - 2 * Stdev
  insert_max: maximum insert size: mean + 2 * Stdev
  Num_Pairs: INT pairs would be count as 1 evidence
  out.improper: output
    out.improper.insertion
    out.improper.deletion

v20170707

EOH

die USAGE if (scalar(@ARGV) !=5 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $input=shift @ARGV;
my $min_allowed_insertsize=shift @ARGV;
my $max_allowed_insertsize=shift @ARGV;
my $numpairs=shift @ARGV;
my $outfile_breaks=shift @ARGV;


### Configurable
### Max mapping times for each read
my $maxmappingtimes=1;
### window maximum size for a group of neighbouring 5 reads
my $opt_disance=10000; 
### window maximum size for the mates of a group of neighbouring 5 reads
my $opt_mdisance=20000;



### Default
my $linenum=0;
my %idhash=();
my %allrefs=();
my %seqlength=();
my $numContigs=0;
my $total_bases=0;
my $total_errors=0;
my $num_ins_Contigs=0;
my $total_ins_bases=0;
my $total_ins_errors=0;
my $num_del_Contigs=0;
my $total_del_bases=0;
my $total_del_errors=0;

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
	my $paired_on_unique=0;
	
	my $uniquecontig=$contigs[0];
	my @arr2=();
	my @arr3=();
	foreach my $contig (@contigs) {
		push (@arr2, $contig) if (exists $idhash{$indread}{$contig}{1});
		push (@arr3, $contig) if (exists $idhash{$indread}{$contig}{2});
		$paired_on_unique=1 if (exists $idhash{$indread}{$contig}{1} and exists $idhash{$indread}{$contig}{2});
	}
	my $maxmapping= max (scalar(@arr2), scalar(@arr3));
	if ($maxmapping>$maxmappingtimes) {### control mapping times
		delete $idhash{$indread};
		next;
	}
	
	unless ($paired_on_unique>0) {
		delete $idhash{$indread};
		next;
	}
	foreach (@contigs) {
		$allrefs{$_}{$indread}++ if (exists $idhash{$indread}{$_}{1} and exists $idhash{$indread}{$_}{2});
	}
}



my %goodregions=();
my %forwardstrand=();
my %reversestrand=();
my %badregions=();
my %insert_region=();
my %deletion_region=();
open (BREAKOUT, " > $outfile_breaks") || die "Error: can not write $outfile_breaks\n";
open (BREAKINS, " > $outfile_breaks.insertion") || die "Error: can not write $outfile_breaks.insertion\n";
open (BREAKDEL, " > $outfile_breaks.deletion") || die "Error: can not write $outfile_breaks.deletion\n";
foreach my $contigs (sort keys %allrefs) {
	my @readnames=keys %{$allrefs{$contigs}};
	next unless (scalar(@readnames)>=$numpairs);
	### Default
	my $numpaired=0; my $numunpaired=0; my $numproperpaired=0; my $orientationpaired=0; my $numimproperpaired=0;
	%goodregions=();
	%forwardstrand=();
	%reversestrand=();
	%badregions=();
	%insert_region=();
	%deletion_region=();
	my %readinfo=();
	my $temp_i=1;
	my $temp_a=1;
	my $temp_b=1;
	
	### detect paired or unpaired
	foreach my $readid (@readnames) {
		if (exists $idhash{$readid}{$contigs}{1} and $idhash{$readid}{$contigs}{2}) {
			$numpaired++;
			if (($idhash{$readid}{$contigs}{1}{'str'} eq '+') and ($idhash{$readid}{$contigs}{2}{'str'} eq '-')) {
				my $distance=$idhash{$readid}{$contigs}{2}{'pos'}-$idhash{$readid}{$contigs}{1}{'pos'};
				if ($distance>=$min_allowed_insertsize and $distance<=$max_allowed_insertsize) {
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
				else {
					if ($distance<$min_allowed_insertsize) {
						$readinfo{$readid}{'del'}++;
					}
					elsif ($distance>$max_allowed_insertsize) {
						$readinfo{$readid}{'ins'}++;
					}
					 $numimproperpaired++;
					 $readinfo{$readid}{'improper'}++;
				}
				$forwardstrand{$temp_i++}=[$readid, 1, $idhash{$readid}{$contigs}{1}{'pos'}];
				$reversestrand{$temp_i++}=[$readid, 2, $idhash{$readid}{$contigs}{2}{'pos'}];
			}
			elsif (($idhash{$readid}{$contigs}{1}{'str'} eq '-') and ($idhash{$readid}{$contigs}{2}{'str'} eq '+')) {
				my $distance=$idhash{$readid}{$contigs}{1}{'pos'}-$idhash{$readid}{$contigs}{2}{'pos'};
				if ($distance>=$min_allowed_insertsize and $distance<=$max_allowed_insertsize) {
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
				else {
					if ($distance<$min_allowed_insertsize) {
						$readinfo{$readid}{'del'}++;
					}
					elsif ($distance>$max_allowed_insertsize) {
						$readinfo{$readid}{'ins'}++;
					}
					$numimproperpaired++;
					$readinfo{$readid}{'improper'}++;
				}
				$forwardstrand{$temp_i++}=[$readid, 2, $idhash{$readid}{$contigs}{2}{'pos'}];
				$reversestrand{$temp_i++}=[$readid, 1, $idhash{$readid}{$contigs}{1}{'pos'}];
			}
			elsif ($idhash{$readid}{$contigs}{1}{'str'} eq $idhash{$readid}{$contigs}{2}{'str'}) {
				$orientationpaired++;
				$readinfo{$readid}{'orientation'}++;
				if (($idhash{$readid}{$contigs}{1}{'str'} eq '+') and ($idhash{$readid}{$contigs}{2}{'str'} eq '+')) {
					$forwardstrand{$temp_i++}=[$readid, 1, $idhash{$readid}{$contigs}{1}{'pos'}];
					$forwardstrand{$temp_i++}=[$readid, 2, $idhash{$readid}{$contigs}{2}{'pos'}];
				}
				elsif (($idhash{$readid}{$contigs}{1}{'str'} eq '-') and ($idhash{$readid}{$contigs}{2}{'str'} eq '-')) {
					$reversestrand{$temp_i++}=[$readid, 1, $idhash{$readid}{$contigs}{2}{'pos'}];
					$reversestrand{$temp_i++}=[$readid, 2, $idhash{$readid}{$contigs}{2}{'pos'}];
				}
				else {
					die "Error: strand\n";
				}
			}
		}
		else {
			$numunpaired++;
			$readinfo{$readid}{'unpaired'}++;
		}
	}
	next unless ($numimproperpaired>=$numpairs);
	my $temphash1=&MergeRegion(\%goodregions);
#	print "Test: $contigs Total ",scalar(@readnames) ," Unpaired $numunpaired Paired $numpaired Properpaired $numproperpaired OrientationProblems $orientationpaired\n";
	
	my $forwref1=&SortHash(\%forwardstrand);
	my $revsref2=&SortHash(\%reversestrand);
	
	for (my $temp_x=0; $temp_x<(scalar(@{$forwref1})-$numpairs); $temp_x++) {
		my @temparr=();
		foreach (my $temp_y=$temp_x; $temp_y<($temp_x+$numpairs); $temp_y++) {
			push (@temparr, ${$forwref1}[$temp_y]);
		}
		next unless (($forwardstrand{$temparr[-1]}[2]-$forwardstrand{$temparr[0]}[2])<=$opt_disance);
		my $total_evidence=0;
		my $total_insertion=0;
		my $total_deletion=0;
		my @matepos=();
		foreach my $temp2 (@temparr) {
			my $thisreadid=$forwardstrand{$temp2}[0];
			my $thismate=3-$forwardstrand{$temp2}[1];
			if (exists $readinfo{$thisreadid} and exists $readinfo{$thisreadid}{'improper'}) {
				$total_evidence++;
				if (exists $idhash{$thisreadid} and exists $idhash{$thisreadid}{$contigs} and exists $idhash{$thisreadid}{$contigs}{$thismate} and ($idhash{$thisreadid}{$contigs}{$thismate}{'str'} eq '-')) {
					push (@matepos, $idhash{$thisreadid}{$contigs}{$thismate}{'pos'});
				}
				else {
					die "Error3\n";
				}
				if (exists $readinfo{$thisreadid}{'ins'}) {
					$total_insertion++;
				}
				if (exists $readinfo{$thisreadid}{'del'}) {
					$total_deletion++;
				}
			}
		}
		if (($total_evidence==$numpairs) and ((($total_insertion==$numpairs) and ($total_deletion==0)) or (($total_insertion==0) and ($total_deletion==$numpairs)))) {
			unless (scalar(@matepos)==$numpairs) {
				die "Error4\n";
			}
			@matepos=sort {$a<=>$b} @matepos;
			if (($matepos[-1]>$matepos[0]) and (($matepos[-1]-$matepos[0])<=$opt_mdisance)) {
				if ($forwardstrand{$temparr[-1]}[2]< $matepos[0]) {
					$badregions{$temp_i++}=[$forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
					if (($total_insertion==$numpairs) and ($total_deletion==0)) {
						$insert_region{$temp_a++}=[$forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
						
					}
					if (($total_insertion==0) and ($total_deletion==$numpairs)) {
						$deletion_region{$temp_b++}=[$forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
					}
				}
				elsif ($matepos[-1]<$forwardstrand{$temparr[0]}[2]) {
					$badregions{$temp_i++}=[$matepos[0], $matepos[-1], $forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2]];
					if (($total_insertion==$numpairs) and ($total_deletion==0)) {
						$insert_region{$temp_a++}=[$matepos[0], $matepos[-1], $forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2]];
						
					}
					if (($total_insertion==0) and ($total_deletion==$numpairs)) {
						$deletion_region{$temp_b++}=[$matepos[0], $matepos[-1], $forwardstrand{$temparr[0]}[2], $forwardstrand{$temparr[-1]}[2]];
					}
				}
				else {
					print STDERR "Warnings: $contigs ", $forwardstrand{$temparr[0]}[2], '-', $forwardstrand{$temparr[-1]}[2], ' vs ', $matepos[0],'-', $matepos[-1], "\n";
				}
			}
		}
	}
	
	for (my $temp_x=0; $temp_x<(scalar(@{$revsref2})-$numpairs); $temp_x++) {
		my @temparr=();
		foreach (my $temp_y=$temp_x; $temp_y<($temp_x+$numpairs); $temp_y++) {
			push (@temparr, ${$revsref2}[$temp_y]);
		}
		next unless (($reversestrand{$temparr[-1]}[2]-$reversestrand{$temparr[0]}[2])<=$opt_disance);
		my $total_evidence=0;
		my $total_insertion=0;
		my $total_deletion=0;
		my @matepos=();
		foreach my $temp2 (@temparr) {
			my $thisreadid=$reversestrand{$temp2}[0];
			my $thismate=3-$reversestrand{$temp2}[1];
			if (exists $readinfo{$thisreadid} and exists $readinfo{$thisreadid}{'improper'}) {
				$total_evidence++;
				if (exists $idhash{$thisreadid} and exists $idhash{$thisreadid}{$contigs} and exists $idhash{$thisreadid}{$contigs}{$thismate} and ($idhash{$thisreadid}{$contigs}{$thismate}{'str'} eq '+')) {
					push (@matepos, $idhash{$thisreadid}{$contigs}{$thismate}{'pos'});
				}
				else {
					die "Error5\n";
				}
				if (exists $readinfo{$thisreadid}{'ins'}) {
					$total_insertion++;
				}
				if (exists $readinfo{$thisreadid}{'del'}) {
					$total_deletion++;
				}
			}
		}
		if (($total_evidence==$numpairs) and ((($total_insertion==$numpairs) and ($total_deletion==0)) or (($total_insertion==0) and ($total_deletion==$numpairs)))) {
			unless (scalar(@matepos)==$numpairs) {
				die "Error6\n";
			}
			@matepos=sort {$a<=>$b} @matepos;
			if (($matepos[-1]>$matepos[0]) and (($matepos[-1]-$matepos[0])<=$opt_mdisance)) {
				if ($reversestrand{$temparr[-1]}[2]< $matepos[0]) {
					$badregions{$temp_i++}=[$reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
					if (($total_insertion==$numpairs) and ($total_deletion==0)) {
						$insert_region{$temp_a++}=[$reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
						
					}
					if (($total_insertion==0) and ($total_deletion==$numpairs)) {
						$deletion_region{$temp_b++}=[$reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2], $matepos[0], $matepos[-1]];
					}
				}
				elsif ($matepos[-1]<$reversestrand{$temparr[0]}[2]) {
					$badregions{$temp_i++}=[$matepos[0], $matepos[-1], $reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2]];
					if (($total_insertion==$numpairs) and ($total_deletion==0)) {
						$insert_region{$temp_a++}=[$matepos[0], $matepos[-1], $reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2]];
						
					}
					if (($total_insertion==0) and ($total_deletion==$numpairs)) {
						$deletion_region{$temp_b++}=[$matepos[0], $matepos[-1], $reversestrand{$temparr[0]}[2], $reversestrand{$temparr[-1]}[2]];
					}
				}
				else {
					print STDERR "Warnings: $contigs ", $reversestrand{$temparr[0]}[2], '-', $reversestrand{$temparr[-1]}[2], ' vs ', $matepos[0],'-', $matepos[-1], "\n";
				}
			}
		}
	}
	
	
	my $mergedbad={};
	$mergedbad=&MergeBadRegions(\%badregions);
	my $mergeinsert={};
	$mergeinsert=&MergeBadRegions(\%insert_region);
	my $mergedeletion={};
	$mergedeletion=&MergeBadRegions(\%deletion_region);
	
	
#	print "Test: \%badregions\n"; print Dumper \%badregions; print "\n"; ### for test ###
#	print "Test: \$mergedbad\n"; foreach (sort {$a<=>$b} keys %{$mergedbad}) {print Dumper ${$mergedbad}{$_};} print "\n";### for test ###
#	exit 0 if (scalar(keys %{$mergedbad})>0);### for test ###
	
	if (scalar(keys %{$mergedbad})>0) {
		$numContigs++;
	}
	else {
		next;
	}
#	print "mergedbad: $contigs\n";### For test ###
	$total_bases+= &CountTotalRegion($mergedbad);
	foreach (sort {$a<=>$b} keys %{$mergedbad}) {
		print BREAKOUT $contigs, "\t", ${$mergedbad}{$_}[0], '-', ${$mergedbad}{$_}[1], ' vs ', ${$mergedbad}{$_}[2],'-', ${$mergedbad}{$_}[3], "\n";
		$total_errors++;
	}
#	print "mergeinsert: $contigs\n";### For test ###
	$num_ins_Contigs++ if (scalar(keys %{$mergeinsert})>0);
	$total_ins_bases+= &CountTotalRegion($mergeinsert);
	foreach (sort {$a<=>$b} keys %{$mergeinsert}) {
		print BREAKINS $contigs, "\t", ${$mergeinsert}{$_}[0], '-', ${$mergeinsert}{$_}[1], ' vs ', ${$mergeinsert}{$_}[2],'-', ${$mergeinsert}{$_}[3], "\n";
		$total_ins_errors++;
#		$total_ins_bases=$total_ins_bases + abs(${$mergeinsert}{$_}[1]-${$mergeinsert}{$_}[0]) + 1;
#		$total_ins_bases=$total_ins_bases + abs(${$mergeinsert}{$_}[3]-${$mergeinsert}{$_}[2]) + 1;
	}
#	print "mergedeletion:$contigs\n";### For test ###
	$num_del_Contigs++ if (scalar(keys %{$mergedeletion})>0);
	$total_del_bases+=&CountTotalRegion($mergedeletion);
	foreach (sort {$a<=>$b} keys %{$mergedeletion}) {
		print BREAKDEL $contigs, "\t", ${$mergedeletion}{$_}[0], '-', ${$mergedeletion}{$_}[1], ' vs ', ${$mergedeletion}{$_}[2],'-', ${$mergedeletion}{$_}[3], "\n";
		$total_del_errors++;
	}
}
close BREAKOUT;
close BREAKINS;
close BREAKDEL;
print "##### SUMMARY 2 #####\n";
print "    Number of problematic errors: $total_errors\n";
print "    Number of problematic contigs: $numContigs\n";
print "    Number of problematic bases: $total_bases\n\n";

print "    Number of insertion errors: $total_ins_errors\n";
print "    Number of insertion contigs: $num_ins_Contigs\n";
print "    Number of insertion bases: $total_ins_bases\n\n";

print "    Number of deletion errors: $total_del_errors\n";
print "    Number of deletion contigs: $num_del_Contigs\n";
print "    Number of deletion bases: $total_del_bases\n\n";

#####################################################################
########################### SUB functions ###########################
#####################################################################



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



sub MergeBadRegions {
	my $MBRhash=shift;
	
	my $MBRrethash={};
	my %MBRtemphash=();
	my @MBRarr=();
	my $MBRindexnum=1;
	my $MBRtest=0;
	
	foreach my $MBRnum (sort {$a <=> $b} keys %{$MBRhash}) {
		$MBRtemphash{${$MBRhash}{$MBRnum}[0]}{$MBRnum}++;
	}
	foreach my $MBRnum (sort {$a <=> $b} keys %MBRtemphash) {
		foreach my $MBRnum2 (sort {$a <=> $b} keys %{$MBRtemphash{$MBRnum}}) {
			push (@MBRarr, $MBRnum2);
		}
	}

	foreach my $MBRnum3 (@MBRarr) {
		my $testoverlap3=0;
		my $MBRtest_overlap1=0;my $MBRtest_overlap2=0;
		my ($MBRtemp1, $MBRtemp2, $MBRtemp3, $MBRtemp4);
		my @MBRarr2=sort {$a <=> $b} keys %{$MBRrethash};
#		print "Test1:  ", ${$MBRhash}{$MBRnum3}[0], '-', ${$MBRhash}{$MBRnum3}[1], ' vs ', ${$MBRhash}{$MBRnum3}[2], '-', ${$MBRhash}{$MBRnum3}[3], "\n";### For test ###
		foreach my $MBRnum4 (@MBRarr2) {
#			print "    $MBRnum4: ", ${$MBRrethash}{$MBRnum4}[0], '-', ${$MBRrethash}{$MBRnum4}[1], ' vs ', ${$MBRrethash}{$MBRnum4}[2], '-', ${$MBRrethash}{$MBRnum4}[3], "\t"; ### For test ###
			($MBRtest_overlap1, $MBRtemp1, $MBRtemp2)=&TestOverlap (${$MBRrethash}{$MBRnum4}[0], ${$MBRrethash}{$MBRnum4}[1], ${$MBRhash}{$MBRnum3}[0], ${$MBRhash}{$MBRnum3}[1]);
###			next if ($MBRtest_overlap1==1);
			($MBRtest_overlap2, $MBRtemp3, $MBRtemp4)=&TestOverlap (${$MBRrethash}{$MBRnum4}[2], ${$MBRrethash}{$MBRnum4}[3], ${$MBRhash}{$MBRnum3}[2], ${$MBRhash}{$MBRnum3}[3]);
###			next if ($MBRtest_overlap2==1);
			if ($MBRtest_overlap1==1 and $MBRtest_overlap2==1) {
				${$MBRrethash}{$MBRnum4}=[$MBRtemp1, $MBRtemp2, $MBRtemp3, $MBRtemp4];
				$testoverlap3++;
				$MBRtest++;
#				print "1\n";### For test ###
				last;
			}
#			else { print "0\n";}### For test ###
		}
		unless ($testoverlap3>0) {
			${$MBRrethash}{$MBRindexnum++}=[${$MBRhash}{$MBRnum3}[0], ${$MBRhash}{$MBRnum3}[1], ${$MBRhash}{$MBRnum3}[2], ${$MBRhash}{$MBRnum3}[3]];
		}
	}
	
	unless ($MBRtest==0) {
		$MBRrethash=&MergeBadRegions($MBRrethash);
	}
	
	return $MBRrethash;
}



sub TestOverlap {
	my ($TOs1, $TOs2, $TOq1, $TOq2) =@_;

	my ($TOn1, $TOn2);
	
#	print "Test2: ", $TOs1, '-', $TOs2, ' vs ', $TOq1, '-', $TOq2; ### For test ###
	
	if (($TOs1>=$TOq1 and $TOs1<=$TOq2) or ($TOs2>=$TOq1 and $TOs2<=$TOq2) or ($TOq1>=$TOs1 and $TOq1<=$TOs2) or ($TOq2>=$TOs1 and $TOq2<=$TOs2)) {
		my @MBRarr=sort {$a<=>$b} ($TOs1, $TOs2, $TOq1, $TOq2);
#		print "\t1\n"; ### For test ###
		return (1, $MBRarr[0], $MBRarr[-1]);
	}
	else {
#		print "\t0\n"; ### For test ###
		return 0;
	}
}



sub CountTotalRegion {
	my $CTRorihash=shift;
	
	my $CTRi=0;
	my @CTRarr=();
	my %CTRtemphash=();
	my $CTRsubinfo='SUB(CountTotalRegion)';
	my $CTSstart=0;
	my $CTRend=0;
	my $CTRret_sum=0;
#	print Dumper $CTRorihash; ### for test ###
	foreach (sort {$a<=>$b} keys %{$CTRorihash}) {
		unless ((defined ${$CTRorihash}{$_}[0]) and (${$CTRorihash}{$_}[0]=~/^\d+$/) and defined ${$CTRorihash}{$_}[1] and (${$CTRorihash}{$_}[1]=~/^\d+$/) and (${$CTRorihash}{$_}[0] <= ${$CTRorihash}{$_}[1])) {
			die $CTRsubinfo, "Error1\n";
		}
		unless ((defined ${$CTRorihash}{$_}[2]) and (${$CTRorihash}{$_}[3]=~/^\d+$/) and (defined ${$CTRorihash}{$_}[3]) and (${$CTRorihash}{$_}[3]=~/^\d+$/) and (${$CTRorihash}{$_}[2] <= ${$CTRorihash}{$_}[3])) {
			die $CTRsubinfo, "Error2\n";
		}
		
		unless (exists $CTRtemphash{${$CTRorihash}{$_}[0]} and $CTRtemphash{${$CTRorihash}{$_}[0]}>=${$CTRorihash}{$_}[1]) {
			$CTRtemphash{${$CTRorihash}{$_}[0]}=${$CTRorihash}{$_}[1];
		}
		unless (exists $CTRtemphash{${$CTRorihash}{$_}[2]} and $CTRtemphash{${$CTRorihash}{$_}[2]}>=${$CTRorihash}{$_}[3]) {
			$CTRtemphash{${$CTRorihash}{$_}[2]}=${$CTRorihash}{$_}[3];
		}
	}
	foreach my $CTRx (sort {$a<=>$b} keys %CTRtemphash) {
		if ($CTRi==0) {
			$CTSstart=$CTRx;
			$CTRend=$CTRtemphash{$CTRx};
			$CTRi++;
		}
		else {
			if ($CTRx>$CTSstart and ($CTRx<=($CTRend+1))) {
				my @CTRarr=($CTSstart, $CTRend, $CTRx, $CTRtemphash{$CTRx});
				@CTRarr=sort {$a<=>$b} @CTRarr;
				$CTRend=$CTRarr[-1];
			}
			elsif ($CTRx>($CTRend+1)) {
				$CTRret_sum=$CTRret_sum + abs($CTRend-$CTSstart) +1;
			}
		}
	}
	$CTRret_sum=$CTRret_sum + abs($CTRend-$CTSstart) +1;
	
	return $CTRret_sum;
}
