#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
use File::Basename;

my $dirname = dirname(__FILE__);

my ($mash_file, $mash_exec, $help, $cutoff, $out, $kmer, $size, $input_file, $id_file, $type_strain);
my $cwd = getcwd;
my @genome_ids; #array of genome identifiers used for row labels
my @print_ids; #array parallel to @genome_ids which indicates the row/col should be printed (1) or skipped (0) based on ANI cutoff to the type strain.
#my $MAX = 1000;
$size = 10000;
$kmer = 17;
$cutoff = "";
GetOptions("mash_file|m=s"=>\$mash_file, "mash_exec|M=s"=>\$mash_exec, "size|s=s"=>\$size, "kmer|k=s"=>\$kmer, "out_prefix|o=s"=>\$out, "help|h|?"=>\$help, "cutoff|c=s"=>\$cutoff, "input_file|f=s"=>\$input_file, "id_file|i=s"=>\$id_file, "type_strain|t"=>\$type_strain);
if ($help || !($mash_file || $input_file) || !$id_file || !$mash_exec || ($mash_file && $input_file) || ($cutoff && !$type_strain)) {
    if ($mash_file && $input_file) {
	print STDERR "Cannot specify both a  mask sketch file and an input file of genome paths!\n\n";
    }
    if ($cutoff && !$type_strain) {
	print STDERR "Cannot specify a percent identity cutoff without specifyin a type strain!\n\n";
    }
    print STDERR "Runs an all versus all MASH of a set of genome fasta files specified by the input_file(-f) or a set of MASH sketches specified by the mash_file(-m)\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome fasta file paths, one per line (required or -m mash sketch file)\n";
    print STDERR "	-i input file of genome idenitifiers, one per line, in the same order and number of identifiers as the genome paths file or mash sketch file (required)\n";
    print STDERR "	-o output file prefix (required if you dont want some generic silly name for your files)\n";
    print STDERR "	-c cutoff value to use to ignore genomes below this ANI level to the type strain genome if desired must only be used with -t\n";
    print STDERR "	-m mash sketch file to compare to (required or -f input file of genome paths)\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-t flag to indicate that the first line of the genome fasta paths file and genome identifiers file shoud be treated as the type strain for cutoff filtering(optional)\n";
    print STDERR "	-k kmer size. Default is 17\n";	
    print STDERR "	-s sketch size. Default is 10000\n";	
    print STDERR "	-? or -h help\n";
		
    exit(0);
}

if (!$out) {
    $out = $cwd . "/temp_mash_out";
}

if ($size =~ /\D/) {
    die ("ERROR: $size the sketch size is not a nonnegative integer\n");
}
if (($size < 1000) || ($size > 100000)) {
    die ("ERROR: $size is outside of the expected sketch size range of 1000-100000\n");
}

if ($cutoff ne "") { 
    if (!looks_like_number($cutoff)) {
	die ("ERROR: $cutoff does not look like a number\n");
    }
    if (($cutoff < 0) || ($cutoff > 100)) {
	die ("ERROR: $cutoff is outside of the expected ANI cutoff range of 0-100\n");
    }
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
my $id_fh;
unless (open ($id_fh, "<", $id_file) )  {
    die ("ERROR: Cannot open genome identifiers input file $id_file!\n");
}
while (my $id = <$id_fh>) {
    chomp ($id); #remove newline
    push @genome_ids, $id;
}
my $genome_ids_size = @genome_ids; #size of the array / number of genome identifiers
if ($input_file) {
    if (!(-e $input_file)) {
	die ( "$input_file the genomes paths input file does not exist!\n");
    }
    if (-d $input_file) {
	die ( "$input_file the genomes paths input file is a directory!\n");
    }
    my $genome_paths_size = `wc -l < $input_file`;
    chomp($genome_paths_size);
    if ($genome_paths_size != $genome_ids_size) {
	die ("ERROR: The number of genome identifiers ($genome_ids_size) is not equal to the number of genome paths ($genome_paths_size)!\n");
    }
    print STDERR "Executing command:\n$mash_exec sketch -k $kmer -s $size -o $out -l $input_file 2>&1";
    my $mash_output = `$mash_exec sketch -k $kmer -s $size -o $out -l $input_file 2>&1`;
    $mash_file = $out . ".msh";
    if (($mash_output =~ /ERROR/) || ($mash_output =~ /WARNING/) || !(-e $mash_file) || !(-s $mash_file)) {
	die ("mash had a problem quitting:\n$mash_output");
    }
}

print STDERR "Executing command:\n$mash_exec dist $mash_file $mash_file | cut -f 3";
my $list = `$mash_exec dist $mash_file $mash_file | cut -f 3`;
my $row_count = 0;
my $col_count = 0;
print STDOUT "ID";
if ($cutoff) {
    while ($list =~ /([^\n\r]+)([\n\r])/g) { #process the first line to output the column headers when cutoff filtering
	my $ani_est = 100 * (1 - $1);
	if ($ani_est >= $cutoff) {
	    $print_ids[$col_count] = 1;
	    print STDOUT "\t$genome_ids[$col_count]";
	} else {
	    $print_ids[$col_count] = 0;
	}
	$col_count++;
	if ($col_count >= $genome_ids_size) {
	    $col_count = 0;
	    print STDOUT "\n";
	    last;
	}
    }
} else {
    foreach my $id (@genome_ids) {
	print STDOUT "\t$id";
    }
    print STDOUT "\n";
}

while ($list =~ /([^\n\r]+)([\n\r])/g) {
    my $ani_est = 100 * (1 - $1);
    if ($col_count == 0) {
	if ($row_count >= $genome_ids_size) {
	    die ("ERROR: too many mash distances returned for the number of genome identifiers ($genome_ids_size)\n");
	}
	if (!$cutoff || $print_ids[$row_count]) {
	    print STDOUT "$genome_ids[$row_count]";
	}
	$row_count++;
    }
    if (!$cutoff || ($print_ids[$col_count] && $print_ids[$row_count])) {
	print STDOUT "\t$ani_est";
    }
    $col_count++;
    if ($col_count >= $genome_ids_size) {
	$col_count = 0;
	print STDOUT "\n";
    }
}
exit (1);
