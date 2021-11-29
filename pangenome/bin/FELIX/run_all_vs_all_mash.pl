#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $dirname = dirname(__FILE__);

my ($mash_file, $mash_exec, $help, $cutoff, $out, $kmer, $size, $input_file, $id_file, $type_strain);
my $cwd = getcwd;
my @genome_ids; #array of genome identifiers used for row labels
my @print_ids; #array parallel to @genome_ids which indicates the row/col should be printed (1) or skipped (0) based on ANI cutoff to the type strain.
my @sum_distances; #sum of distances to other genomes used to determine medoid
#my $MAX = 1000;
my $mean_all;
my $num_all = 0;
my $total_all = 0;
my $min_all = 101;
my $max_all = -1;
my $mean_kept;
my $num_kept = 0;
my $total_kept = 0;
my $min_kept = 101;
my $max_kept = -1;
my $mean_discard;
my $num_discard = 0;
my $total_discard = 0;
my $min_discard = 101;
my $max_discard = -1;
my $mean_cross;
my $num_cross = 0;
my $total_cross = 0;
my $min_cross = 101;
my $max_cross = -1;
my $genomes_kept = 0;
my $genomes_discard = 0;
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
    print STDERR "Runs an all versus all MASH of a set of genome fasta files specified by the input_file(-f) or a set of MASH sketches specified by the mash_file(-m)\nOutputs an ANI matrix to STDOUT and various information to STDERR\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome fasta file paths, one per line (required or -m mash sketch file)\n";
    print STDERR "	-i input file of genome idenitifiers, one per line, in the same order and number of identifiers as the genome paths file or mash sketch file (required)\n";
    print STDERR "	-o output file prefix (required if you dont want some generic silly name for your mash sketch file)\n";
    print STDERR "	-c cutoff value to use to ignore genomes below this ANI level to the type strain genome if desired must only be used with -t\n";
    print STDERR "	-m mash sketch file to compare to (required or -f input file of genome paths)\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-t flag to indicate that the first line of the genome fasta paths file and genome identifiers file shoud be treated as the type strain for cutoff filtering(optional)\n";
    print STDERR "	-k kmer size. Default is 17\n";	
    print STDERR "	-s sketch size. Default is 10000\n";	
    print STDERR "	-? or -h help\n";
		
    exit(1);
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
my $row_count = 0;
while (my $id = <$id_fh>) {
    chomp ($id); #remove newline
    push @genome_ids, $id;
    $sum_distances[$row_count] = 0;
    $row_count++;
}
$row_count = 0;
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
    my $mash_sketch_out = $out . ".mash_sketch_error";
    print STDERR "Executing command:\n$mash_exec sketch -k $kmer -s $size -o $out -l $input_file >& $mash_sketch_out\n";
    `$mash_exec sketch -k $kmer -s $size -o $out -l $input_file >& $mash_sketch_out`;
    $mash_file = $out . ".msh";
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
}

my $mash_out_file = $out . ".mash_dist_out";
print STDERR "Executing command:\n$mash_exec dist -t $mash_file $mash_file > $mash_out_file\n";
`$mash_exec dist -t $mash_file $mash_file > $mash_out_file`;
if (!(-e $mash_out_file) || !(-s $mash_out_file)) {
    die ( "ERROR: mash dist command failed: $mash_out_file does not exist or is zero size!\n");
}
my $dist_fh;
unless (open ($dist_fh, "<", $mash_out_file) )  {
    die ("ERROR: Cannot open mash distances file $mash_out_file!\n");
}
my $col_count = 0;
print STDOUT "ID";
if ($cutoff) {
    my $dist = <$dist_fh>; #skip header line
    $dist = <$dist_fh>;
    my $first = 1;
    while ($dist =~ /([^\t\n\r]+)([\t\n\r])/g) { #process the tab delimited distances
	if ($first) {
	    $first = 0;
	    next; #skip header
	}
	if ($col_count >= $genome_ids_size) {
	    die ("ERROR: The number of distances in a single row of the tabular mash output exceeds the number of genome identifiers provided ($genome_ids_size).\n");
	}
	my $ani_est = 100 * (1 - $1);
	if ($ani_est >= $cutoff) {
	    $print_ids[$col_count] = 1;
	    print STDOUT "\t$genome_ids[$col_count]";
	    $genomes_kept++;
	} else {
	    $print_ids[$col_count] = 0;
	    print STDERR "Genome $genome_ids[$col_count] is being filtered out for ANI ($ani_est) below cutoff ($cutoff)\n";
	    $genomes_discard++;
	}
	$col_count++;
    }
    $col_count = 0;
    print STDOUT "\n";
    #pos($dist) = 0; #reset global regex pattern matching to begining of $list
    close($dist_fh);
    unless (open ($dist_fh, "<", $mash_out_file) )  {
	die ("ERROR: Cannot reopen mash distances file $mash_out_file!\n");
    }
} else {
    foreach my $id (@genome_ids) {
	print STDOUT "\t$id";
    }
    print STDOUT "\n";
}

my $dist = <$dist_fh>; #skip header line
while ($dist = <$dist_fh>) {
    if ($row_count >= $genome_ids_size) {
	die ("ERROR: too many mash distances rows for the number of genome identifiers ($genome_ids_size)\n");
    }
    if (!$cutoff || $print_ids[$row_count]) {
	print STDOUT "$genome_ids[$row_count]";
    }
    my $first = 1;
    while ($dist =~ /([^\t\n\r]+)([\t\n\r])/g) { #process the tab delimited distances
	if ($first) {
	    $first = 0;
	    next; #skip header
	}
	if ($col_count >= $genome_ids_size) {
	    die ("ERROR: The number of distances in a single row of the tabular mash output exceeds the number of genome identifiers provided ($genome_ids_size).\n");
	}
	my $ani_est = 100 * (1 - $1);
	if (!$cutoff || ($print_ids[$col_count] && $print_ids[$row_count])) {
	    print STDOUT "\t$ani_est";
	}
	if ($row_count < $col_count) {#don't count self matches in statistics - matrix is symmetric so just use upper half
	    $sum_distances[$row_count] += $1;
	    $sum_distances[$col_count] += $1; #have to do this because we are only doing upper half of symmetric matrix
	    $num_all++;
	    $total_all += $ani_est;
	    if ($ani_est < $min_all) {
		$min_all = $ani_est;
	    }
	    if ($ani_est > $max_all) {
		$max_all = $ani_est;
	    }
	    if ($cutoff) {
		if ($print_ids[$col_count] && $print_ids[$row_count]) {
		    $num_kept++;
		    $total_kept += $ani_est;
		    if ($ani_est < $min_kept) {
			$min_kept = $ani_est;
		    }
		    if ($ani_est > $max_kept) {
			$max_kept = $ani_est;
		    }
		} elsif ($print_ids[$col_count] || $print_ids[$row_count]) {
		    $num_cross++;
		    $total_cross += $ani_est;
		    if ($ani_est < $min_cross) {
			$min_cross = $ani_est;
		    }
		    if ($ani_est > $max_cross) {
			$max_cross = $ani_est;
		    }
		} else {
		    $num_discard++;
		    $total_discard += $ani_est;
		    if ($ani_est < $min_discard) {
			$min_discard = $ani_est;
		    }
		    if ($ani_est > $max_discard) {
			$max_discard = $ani_est;
		    }
		}
	    }
	}
	$col_count++;
    }
    if (!$cutoff || $print_ids[$row_count]) {
	print STDOUT "\n";
    }
    $col_count = 0;
    $row_count++;
}
close($dist_fh);
if ($row_count < $genome_ids_size) {
    die ("ERROR: too few mash distances rows ($row_count) for the number of genome identifiers ($genome_ids_size)\n");
}
`rm $mash_out_file`;
my $index = 0;
my $min_sum_dist = $genome_ids_size + 1;
my $medoid_genome_id;
foreach my $sum_dist (@sum_distances) {
    if ($sum_dist < $min_sum_dist) {
	$min_sum_dist = $sum_dist;
	$medoid_genome_id = $genome_ids[$index];
    }
    $index++;
}
my $medoid_genome_ANI = 100 * (1 - ($min_sum_dist / ($genome_ids_size - 1)));
print STDERR "Medoid genome: $medoid_genome_id\nMean pairwise ANI for medoid genome: $medoid_genome_ANI\n"; 
if ($num_all > 0) {
    $mean_all = $total_all / $num_all;
} else {
    $mean_all = 0;
}
print STDERR "Mean, min, max pairwise ANI for all $genome_ids_size genomes: $mean_all, $min_all, $max_all\n";
if ($cutoff && ($genomes_discard > 0)) {
    if ($num_kept > 0) {
	$mean_kept = $total_kept / $num_kept;
    } else {
	$mean_kept = 0;
    }
    if ($num_discard > 0) {
	$mean_discard = $total_discard / $num_discard;
    } else {
	$mean_discard = 0;
    }
    if ($num_cross > 0) {
	$mean_cross = $total_cross / $num_cross;
    } else {
	$mean_cross = 0;
    }
    print STDERR "For cutoff ($cutoff) and type strain $genome_ids[0]: $genomes_kept kept, $genomes_discard discarded\n";
    print STDERR "Mean, min, max pairwise ANI for just kept genomes: $mean_kept, $min_kept, $max_kept\n";
    print STDERR "Mean, min, max pairwise ANI for just discarded genomes: $mean_discard, $min_discard, $max_discard\n";
    print STDERR "Mean, min, max pairwise ANI between kept and discarded genomes: $mean_cross, $min_cross, $max_cross\n";
}
exit (0);
