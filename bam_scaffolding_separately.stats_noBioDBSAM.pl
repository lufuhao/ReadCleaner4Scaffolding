#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/VerifyCigarLength CalCigarRefLength/;
use constant USAGE =><<EOH;

usage: $0 R1.bam R2.bam out.file

    Evaluate BAM files before scaffolding

    v20160919

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



###INPUT
die "Error: invalid R1.BAM file\n" unless (defined $ARGV[0] and -s $ARGV[0]);
die "Error: invalid R2.BAM file\n" unless (defined $ARGV[1] and -s $ARGV[1]);



###OUTPUT
die "Error: invalid output.file\n" unless (defined $ARGV[2] and $ARGV[2]=~/^\S+$/);
unlink $ARGV[2] if (-e $ARGV[2]);



### Default
my $path_samtools='samtools';
my $bamin1=$ARGV[0];
my $bamin2=$ARGV[1];
my $statout=$ARGV[2];


### OPEN
if ($bamin1=~/\.sam$/i) {
	open (BAMINR1, "$path_samtools view -S $bamin1 | ") || die "Error: can not open SAM file: $bamin1\n";
}
elsif ($bamin1=~/\.bam$/i) {
	open (BAMINR1, "$path_samtools view $bamin1 | ") || die "Error: can not open BAM file: $bamin1\n";
}
else {
	die "Error: can not guess R1.BAM file format from file.suffix: .sam or .bam?\n";
}
die "Error: can not guess R2.BAM file format from file.suffix: .sam or .bam?\n" unless ($bamin2=~/\.bam$/i or $bamin2=~/\.bam$/i);

open (OUTPUT, "> $statout") || die "Error: can not write out.file\n";


### R1
my $linenum=0;
my $mate=1;
while (my $line=<BAMINR1>) {
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
close BAMINR1;



###R2
if ($bamin2=~/\.sam$/i) {
	open (BAMINR2, "$path_samtools view -h -S $bamin2 | ") || die "Error: can not open R2.SAM file: $bamin2\n";
}
elsif ($bamin2=~/\.bam$/i) {
	open (BAMINR2, "$path_samtools view -h $bamin2 | ") || die "Error: can not open R2.BAM file: $bamin2\n";
}
else {
	die "Error: can not guess R2.BAM file format from file.suffix: .sam or .bam?\n";
}
$linenum=0;
$mate=2;
while (my $line=<BAMINR2>) {
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

close BAMINR2;
close OUTPUT;
