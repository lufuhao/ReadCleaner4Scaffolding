#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/VerifyCigarLength CalCigarRefLength/;
use constant USAGE =><<EOH;

usage: $0 file.bam out.file

    Evaluate BAM files before scaffolding

    v20160919

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



if ($ARGV[0]=~/\.sam$/i) {
	open (BAMIN, "samtools view -h -S $ARGV[0] | ") || die "Error: can not open SAM file: $ARGV[0]\n";
}
elsif ($ARGV[0]=~/\.bam$/i) {
	open (BAMIN, "samtools view -h $ARGV[0] | ") || die "Error: can not open SAM file: $ARGV[0]\n";
}
else {
	die "Error: can not open guess BAM file format from file.suffix: .sam or .bam?\n";
}

#die "Error: invalid reference.fasta file\n" unless (defined $ARGV[1] and -s $ARGV[1]);

die "Error: invalid out.file\n" unless (defined $ARGV[1] and $ARGV[1]=~/^\S+$/);
die "Error: out.file existed: $ARGV[1]\n" if (-e $ARGV[1]);


open (OUTPUT, "> $ARGV[1]") || die "Error: can not write out.file\n";
my $linenum=0;
while (my $line=<BAMIN>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^\@/);
	my @arr=split(/\t/, $line);
	my $readname=$arr[0];
	my $flag=$arr[1];
	my $seqid=$arr[2];
	my $start=$arr[3];
	my $cigar=$arr[5];
	my $queryseq=$arr[9];
	
	unless (VerifyCigarLength($cigar, length($queryseq))) {
		print OUTPUT "Error:cigar_length_error at line($linenum): $line\n";
		next;
	}
	my $mate=0;
	if ($flag & 0x0040) {
		$mate=1;
	}
	elsif ($flag & 0x0080) {
		$mate=2;
	}
	elsif ($readname=~/^s\d+a$/) {###For SOPRA
		$mate=1;
	}
	elsif ($readname=~/^s\d+b$/) {###For SOPRA
		$mate=2;
	}
	my $strand='NaN';
	if ($flag & 0x0010) {
		$strand='-';
	}
	else {
		$strand='+';
	}
	my $endlength= CalCigarRefLength($cigar);
	my $end=$start+$endlength-1;
	print OUTPUT $readname, "\t", $mate, "\t", $seqid, "\t", $start, "\t", $end, "\t", $strand, "\t", $cigar, "\n";
}
close OUTPUT;
close BAMIN;
