#!/usr/bin/env perl
use strict;
use warnings;
use Bio::DB::Sam;
use constant USAGE =><<EOH;

usage: $0 file.bam reference.fasta out.file

    Evaluate BAM files before scaffolding

    v20160919

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


die "Error: invalid BAM/SAM file\n" unless (defined $ARGV[0] and -s $ARGV[0]);
die "Error: invalid reference.fasta file\n" unless (defined $ARGV[1] and -s $ARGV[1]);
die "Error: invalid out.file\n" unless (defined $ARGV[2] and $ARGV[2]=~/^\S+$/);
die "Error: out.file existed: $ARGV[2]\n" if (-e $ARGV[2]);

my $bamobj=Bio::DB::Sam -> new (-bam => "$ARGV[0]", -fasta => "$ARGV[1]");

my @seqids=$bamobj->seq_ids;
open (OUTPUT, "> $ARGV[2]") || die "Error: can not write out.file\n";
foreach my $seq (@seqids) {
	my @alignment=$bamobj->get_features_by_location(-seq_id => "$seq");
	foreach my $idvalign (@alignment) {
		my $query_name=$idvalign->name;
		my $query_start = $idvalign->query->start;     
		my $query_end   = $idvalign->query->end;
		my $seqid  = $idvalign->seq_id;
		my $start  = $idvalign->start;
		my $end    = $idvalign->end;
		my $strand = $idvalign->strand;
		my $cigar  = $idvalign->cigar_str;
		
		print OUTPUT $query_name, "\t", $query_start, "\t", $query_end, "\t", $seqid, "\t", $start, "\t", $end, "\t", $strand, "\t", $cigar, "\n";
	}
}
close OUTPUT;
