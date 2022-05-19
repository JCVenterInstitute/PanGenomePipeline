#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $dirname = dirname(__FILE__);

my ($mash_exec, $help, $out_dir, $kmer, $size, $input_file, $id_file);
my $cwd = getcwd;
my @genome_ids; #array of genome identifiers used for column labels
$size = 10000;
$kmer = 17;
$out_dir = $cwd;
GetOptions("mash_exec|M=s"=>\$mash_exec, "size|s=s"=>\$size, "kmer|k=s"=>\$kmer, "out_dir|o=s"=>\$out_dir, "help|h|?"=>\$help, "input_file|f=s"=>\$input_file, "id_file|i=s"=>\$id_file);
if ($help || !$input_file || !$mash_exec || !$id_file) {
    if (!$help) {
	if (!$input_file) {
	    print STDERR "Must specify an input file of genome paths!\n\n";
	}
	if (!$mash_exec) {
	    print STDERR "Must specify a path to the MASH executable!\n\n";
	}
	if (!$id_file) {
	    print STDERR "Must specify an identifier file of genome identifiers!\n\n";
	}
    }
    print STDERR "Creates MASH sketch files in out_dir(-o) for a set of genome fasta files specified by the input_file(-f)\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome fasta file paths, one per line (required)\n";
    print STDERR "	-i input file of genome idenitifiers, one per line, in the same order and number of identifiers as the genome fasta file paths (required)\n";
    print STDERR "	-o output directory (required if you do not want your MASH sketch files in the current directory)\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-k kmer size. Default is 17\n";	
    print STDERR "	-s sketch size. Default is 10000\n";	
    print STDERR "	-? or -h help\n";
		
    exit(1);
}

if ($size =~ /\D/) {
    die ("ERROR: $size the sketch size is not a nonnegative integer\n");
}
if (($size < 1000) || ($size > 100000)) {
    die ("ERROR: $size is outside of the expected sketch size range of 1000-100000\n");
}

if ($kmer =~ /\D/) {
    die ("ERROR: $kmer the k-mer size is not a nonnegative integer\n");
}
if (($kmer < 7) || ($kmer > 33)) {
    die ("ERROR: $kmer is outside of the expected k-mer size range of 7-33\n");
}

if (!(-e $mash_exec)) {
    die ("$mash_exec mash executable file does not exist!\n");
}
if (-d $mash_exec) {
    die ("$mash_exec mash executable file is a directory!\n");
}
if (!(-x $mash_exec)) {
    die ("$mash_exec mash executable file is not executable!\n");
}
if (!(-e $id_file)) {
    die ( "$id_file the genome identifers input file does not exist!\n");
}
if (-d $id_file) {
    die ( "$id_file the genome identifiers input file is a directory!\n");
}
if (!(-e $input_file)) {
    die ( "$input_file - the genome paths input file does not exist!\n");
}
if (-d $input_file) {
    die ( "$input_file - the genome pathss input file is a directory!\n");
}
if (!(-e $out_dir)) {
    die ( "$out_dir - the output directory for the sketch files does not exist!\n");
}
if (!(-d $out_dir)) {
    die ( "$out_dir - the output directory for the sketch files is not a directory!\n");
}
if (substr($out_dir, -1, 1) ne "/") {
    $out_dir .= "/";
}
if (substr($out_dir, 0, 1) ne "/") {
    $out_dir = $cwd . "/$out_dir";
}
my $id_fh;
unless (open ($id_fh, "<", $id_file) )  {
    die ("ERROR: Cannot open genome identifiers input file $id_file!\n");
}
while (my $id = <$id_fh>) {
    chomp ($id); #remove newline
    push @genome_ids, $id;
}
my $genome_ids_size = @genome_ids; #size of the array / number of genome identifiers
my $genome_paths_size = `wc -l < $input_file`;
chomp($genome_paths_size);
if ($genome_paths_size != $genome_ids_size) {
    die ("ERROR: The number of genome identifiers ($genome_ids_size) is not equal to the number of genome MASH sketch paths ($genome_paths_size)!\n");
}
my $input_fh;
unless (open ($input_fh, "<", $input_file) )  {
    die ("ERROR: Cannot open genome paths input file $input_file!\n");
}
my $index = 0;
while (my $genome_path = <$input_fh>) {
    chomp ($genome_path); #remove newline
    if (!(-e $genome_path)) {
	die ( "$genome_path the genome path file does not exist!\n");
    }
    if (-d $genome_path) {
	die ( "$genome_path the genome path file is a directory!\n");
    }
    my ($gp_id, $gp_dir, $gp_suf) = fileparse($genome_path, qr/\.[^.]*/);
    my $mash_sketch_out = $out_dir . $gp_id . ".mash_sketch_error";
    my $out = $out_dir . $gp_id;
    print STDERR "Executing command:\n$mash_exec sketch -k $kmer -s $size -o $out -I $genome_ids[$index] $genome_path >& $mash_sketch_out\n";
    `$mash_exec sketch -k $kmer -s $size -o $out -I $genome_ids[$index] $genome_path >& $mash_sketch_out`;
    my $mash_file = $out . ".msh";
    if (!(-e $mash_file) || !(-s $mash_file)) {
	die ("mash had a problem quitting:\nNo or zero size .msh sketch file created\nSee file: $mash_sketch_out for complete mash sketch output\n");
    }
    my $sketch_fh;
    unless (open ($sketch_fh, "<", $mash_sketch_out) )  {
	die ("ERROR: Cannot open mash sketch output file $mash_sketch_out!\n");
    }
    while (my $mash_output = <$sketch_fh>) {
	if (($mash_output =~ /ERROR/) || ($mash_output =~ /WARNING/)) {
	    die ("mash had a problem quitting:\n$mash_output\nSee file: $mash_sketch_out for complete mash sketch output\n");
	}
    }
    close($sketch_fh);
    `rm $mash_sketch_out`;
    $index++;
}
exit (0);
