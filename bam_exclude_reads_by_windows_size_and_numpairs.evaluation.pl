#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw /Dumper/;
use constant USAGE =><<EOH;

usage: $0 input order output

v20161205

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $input=$ARGV[0];
my $order=$ARGV[1];
my $output=$ARGV[2];


my %idhash=();


open (CHRORDER, "< $order") || die "Error: can not open order\n";
while (my $line=<CHRORDER>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	unless (exists $idhash{$arr[0]}) {
		$idhash{$arr[0]}{'beg'}=$arr[4];
		$idhash{$arr[0]}{'end'}=$arr[5];
		$idhash{$arr[0]}{'str'}=$arr[7];
	}
	else {
		die "Error: repeated ID:$arr[0]\n"
	}
}
close CHRORDER;

#print Dumper \%idhash;exit 0; ### For test ###

open (LINKINPUT, " < $input") || die "Error: can not open links input\n";
open (LINKOUTPUT, " > $output ") || die "Error: can not write oputput\n";
while (my $line=<LINKINPUT>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	my $strandol='';
	if (exists $idhash{$arr[0]} and exists $idhash{$arr[2]}) {
		
		my $distance=0;
		if ($idhash{$arr[0]}{'str'} eq '+') {
			$strandol.='+>';
		}
		else {
			$strandol.='<-';
		}
		if ($idhash{$arr[2]}{'str'} eq '+') {
			$strandol.='+>';
		}
		else {
			$strandol.='<-';
		}
		if ($idhash{$arr[0]}{'beg'}>$idhash{$arr[2]}{'end'} or $idhash{$arr[2]}{'beg'}>$idhash{$arr[0]}{'end'}) {
			if ($idhash{$arr[0]}{'beg'}>$idhash{$arr[2]}{'end'}) {
				$distance=$idhash{$arr[0]}{'beg'}-$idhash{$arr[2]}{'end'};
			}
			if ($idhash{$arr[2]}{'beg'}>$idhash{$arr[0]}{'end'}) {
				$distance=$idhash{$arr[2]}{'beg'}-$idhash{$arr[0]}{'end'}-1;
			}
		}
		elsif ($idhash{$arr[0]}{'beg'}<=$idhash{$arr[2]}{'end'} and $idhash{$arr[0]}{'beg'}>=$idhash{$arr[2]}{'beg'}) {
			if ($idhash{$arr[0]}{'end'} > $idhash{$arr[2]}{'end'}) {
				$distance=-($idhash{$arr[2]}{'end'}-$idhash{$arr[0]}{'beg'}+1);
			}
			else {
				$distance=-($idhash{$arr[0]}{'end'}-$idhash{$arr[0]}{'beg'}+1);
			}
		}
		elsif ($idhash{$arr[2]}{'beg'}<=$idhash{$arr[0]}{'end'} and $idhash{$arr[2]}{'beg'}>=$idhash{$arr[0]}{'beg'}) {
			if ($idhash{$arr[2]}{'end'} > $idhash{$arr[0]}{'end'}) {
				$distance=-($idhash{$arr[0]}{'end'}-$idhash{$arr[2]}{'beg'}+1);
			}
			else {
				$distance=-($idhash{$arr[2]}{'end'}-$idhash{$arr[2]}{'beg'}+1);
			}
		}
		$line= $line."\t$strandol\t$distance";
	}
	else {
		if (exists $idhash{$arr[0]}) {
			if ($idhash{$arr[0]}{'str'} eq '+') {
				$strandol.='+>';
			}
			else {
				$strandol.='<-';
			}
		}
		else {
			$strandol.='XX';
		}
		if ($idhash{$arr[2]}) {
			if ($idhash{$arr[2]}{'str'} eq '+') {
				$strandol.='+>';
			}
			else {
				$strandol.='<-';
			}
		}
		else {
			$strandol.='XX';
		}
		$line= $line."\t$strandol\tNaN";
	}
	
	if (exists $idhash{$arr[0]}) {
		$line=$line."\t".$idhash{$arr[0]}{'beg'}."\t".$idhash{$arr[0]}{'end'}."\t".$idhash{$arr[0]}{'str'};
	}
	else {
		$line=$line."\tNaN\tNaN\tNaN";
	}
	if (exists $idhash{$arr[2]}) {
		$line=$line."\t".$idhash{$arr[2]}{'beg'}."\t".$idhash{$arr[2]}{'end'}."\t".$idhash{$arr[2]}{'str'};
	}
	else {
		$line=$line."\tNaN\tNaN\tNaN";
	}
	print LINKOUTPUT $line, "\n";
}
close LINKINPUT;
close LINKOUTPUT;
