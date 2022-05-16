#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $dirname = dirname(__FILE__);

my ($mash_exec, $help, $cutoff, $out, $input_file, $id_file, $type_strain, $type_strain_id);
my $cwd = getcwd;
my @genome_ids; #array of genome identifiers used for column labels
my @print_ids; #array parallel to @genome_ids which indicates the column should be printed (1) or skipped (0) based on ANI cutoff to the type strain.
my @distances; #distances from type strain to other genomes used to determine median and other stats
my @indices; #indices of parallel arrays from 0 -> N-1
my @ordered_indices; #indices of parallel arrays from 0 -> N-1 ordered by distance from type strain smallest first
my @diffs; #differences between ordered ditances
my @ordered_diffs; #ordered differences between ordered ditances smallest first
my $mean_all;
my $num_all = 0;
my $total_all = 0;
my $min_all = 101;
my $max_all = -1;
my $median_all;
my $mean_kept;
my $num_kept = 0;
my $total_kept = 0;
my $min_kept = 101;
my $max_kept = -1;
my $median_kept;
my $mean_discard;
my $num_discard = 0;
my $total_discard = 0;
my $min_discard = 101;
my $max_discard = -1;
my $median_discard;
my $genomes_kept = 0;
my $genomes_discard = 0;
$cutoff = "";
GetOptions("mash_exec|M=s"=>\$mash_exec, "out_prefix|o=s"=>\$out, "help|h|?"=>\$help, "cutoff|c=s"=>\$cutoff, "input_file|f=s"=>\$input_file, "id_file|i=s"=>\$id_file, "type_strain|t=s"=>\$type_strain, "type_strain_id|T=s"=>\$type_strain_id);
if ($help || !$input_file || !$id_file || !$mash_exec || !$type_strain || !$type_strain_id) {
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
	if (!$type_strain) {
	    print STDERR "Must specify a path to the type strain genome MASH sketch!\n\n";
	}
	if (!$type_strain_id) {
	    print STDERR "Must specify an identifier for the type strain genome!\n\n";
	}
    }
    print STDERR "Runs a type strain versus all MASH of a set of genome MASH sketches specified by the input_file(-f)\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome MASH sketch file paths, one per line (required)\n";
    print STDERR "	-i input file of genome idenitifiers, one per line, in the same order and number of identifiers as the genome MASH sketch paths file (required)\n";
    print STDERR "	-o output file prefix (required if you do not want some generic silly name for your output files)\n";
    print STDERR "	-c cutoff value to use to ignore genomes below this ANI level to the type strain genome if desired\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-t genome MASH sketch file path to the type strain\n";
    print STDERR "	-T genome identifier for the type strain\n";
    print STDERR "	-? or -h help\n";
		
    exit(1);
}

if (!$out) {
    $out = $cwd . "/temp_mash_out";
}

if ($cutoff ne "") { 
    if (!looks_like_number($cutoff)) {
	die ("ERROR: the ANI cutoff $cutoff does not look like a number\n");
    }
    if (($cutoff < 0) || ($cutoff > 100)) {
	die ("ERROR: $cutoff is outside of the expected ANI cutoff range of 0-100\n");
    }
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
    die ( "$input_file the genomes MASH sketch paths input file does not exist!\n");
}
if (-d $input_file) {
    die ( "$input_file the genomes MASH sketch paths input file is a directory!\n");
}
if (!(-e $type_strain)) {
    die ( "$type_strain the type strain genome MASH sketch path file does not exist!\n");
}
if (-d $type_strain) {
    die ( "$type_strain the type strain genome MASH sketch path file is a directory!\n");
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

my $mash_out_file = $out . ".mash_dist_out";
print STDERR "Executing command:\n$mash_exec dist -t $type_strain `cat $input_file` > $mash_out_file\n";
`$mash_exec dist -t $type_strain \`cat $input_file\` > $mash_out_file`;
if (!(-e $mash_out_file) || !(-s $mash_out_file)) {
    die ( "ERROR: mash dist command failed: $mash_out_file does not exist or is zero size!\n");
}
my $dist_fh;
unless (open ($dist_fh, "<", $mash_out_file) )  {
    die ("ERROR: Cannot open mash distances file $mash_out_file!\n");
}
print STDOUT "#Type strain $type_strain_id\n";
my $dist = <$dist_fh>; #skip header line
my $row_count = 0;
while ($dist = <$dist_fh>) { #process the MASH lines
    if ($row_count >= $genome_ids_size) {
	die ("ERROR: The number of distances in the tabular mash output exceeds the number of genome identifiers provided ($genome_ids_size).\n");
    }
    my @fields = split(/\s/, $dist);
    my $num_fields = @fields;
    if ($num_fields != 2) {
	die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
    }
    my $ani_est = 100 * (1 - $fields[1]);
    if (!$cutoff || ($ani_est >= $cutoff)) {
	$print_ids[$row_count] = 1;
	print STDOUT "$genome_ids[$row_count]\t$ani_est\n";
	$genomes_kept++;
    } else {
	$print_ids[$row_count] = 0;
	print STDERR "Genome $genome_ids[$row_count] is being filtered out for ANI ($ani_est) below cutoff ($cutoff)\n";
	$genomes_discard++;
    }
    $distances[$row_count] += $fields[1];
    $indices[$row_count] = $row_count;
    $num_all++;
    $total_all += $ani_est;
    if ($ani_est < $min_all) {
	$min_all = $ani_est;
    }
    if ($ani_est > $max_all) {
	$max_all = $ani_est;
    }
    if ($cutoff) {
	if ($print_ids[$row_count]) {
	    $num_kept++;
	    $total_kept += $ani_est;
	    if ($ani_est < $min_kept) {
		$min_kept = $ani_est;
	    }
	    if ($ani_est > $max_kept) {
		$max_kept = $ani_est;
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
    $row_count++;
}
close($dist_fh);
`rm $mash_out_file`;
@ordered_indices = sort { $distances[$a] <=> $distances[$b] || $a <=> $b } @indices; # sort indices from smallest to largest distance to type strain
if ($num_all > 0) {
    $mean_all = $total_all / $num_all;
    $median_all = 100 * (1 - (($num_all % 2) ? $distances[$ordered_indices[($num_all / 2)]] : (($distances[$ordered_indices[(($num_all / 2) - 1)]] + $distances[$ordered_indices[($num_all / 2)]]) / 2)));
} else {
    $mean_all = 0;
    $median_all = 0;
}
print STDERR "Mean, median, min, max pairwise ANI for all $genome_ids_size genomes: $mean_all, $median_all, $min_all, $max_all\n";
if ($cutoff && ($genomes_discard > 0)) {
    if ($num_kept > 0) {
	$mean_kept = $total_kept / $num_kept;
	$median_kept = 100 * (1 - (($num_kept % 2) ? $distances[$ordered_indices[($num_kept / 2)]] : (($distances[$ordered_indices[(($num_kept / 2) - 1)]] + $distances[$ordered_indices[($num_kept / 2)]]) / 2)));
    } else {
	$mean_kept = 0;
	$median_kept = 0;
    }
    if ($num_discard > 0) {
	$mean_discard = $total_discard / $num_discard;
	$median_discard = 100 * (1 - (($num_discard % 2) ? $distances[$ordered_indices[$num_kept + ($num_discard / 2)]] : (($distances[$ordered_indices[$num_kept + (($num_discard / 2) - 1)]] + $distances[$ordered_indices[$num_kept + ($num_discard / 2)]]) / 2)));
    } else {
	$mean_discard = 0;
	$median_discard = 0;
    }
    print STDERR "For cutoff ($cutoff) and type strain $type_strain_id: $genomes_kept kept, $genomes_discard discarded\n";
    print STDERR "Mean, median, min, max pairwise ANI for just kept genomes: $mean_kept, $median_kept, $min_kept, $max_kept\n";
    print STDERR "Mean, median, min, max pairwise ANI for just discarded genomes: $mean_discard, $median_discard, $min_discard, $max_discard\n";
}
if ($num_kept > 0) {
    print STDERR "Percentiles for ANI to the type strain for kept genomes:\n";
    my $slice = $num_kept / 100;
    for (my $i=1; $i <= 100; $i++) {
	my $index = int(($i * $slice) + 0.499);
	my $percentile = 100 * (1 - $distances[$ordered_indices[$index]]);
	print STDERR "$i:$percentile\n";
    }
    for (my $i=1; $i < $num_kept; $i++) {
	$diffs[$i-1] = $distances[$ordered_indices[$i-1]] - $distances[$ordered_indices[$i]];
    }
    my $num_diffs = @diffs;
    if ($num_diffs > 0) {
	my $total_diff = 0;
	my $min_diff = 100;
	my $max_diff = 0;
	my $mean_diff;
	my $median_diff;
	foreach my $diff (@diffs) {
	    my $ani_est = 100 * (1 - $diff);
	    $total_diff += $ani_est;
	    if ($ani_est < $min_diff) {
		$min_diff = $ani_est;
	    }
	    if ($ani_est > $max_diff) {
		$max_diff = $ani_est;
	    }
	}
	@ordered_diffs = sort { $a <=> $b } @diffs; # sort diffs from smallest to largest
	$mean_diff = $total_diff / $num_diffs;
	$median_diff = 100 * (1 - (($num_diffs % 2) ? $ordered_diffs[($num_diffs / 2)] : (($ordered_diffs[(($num_diffs / 2) - 1)] + $ordered_diffs[($num_diffs / 2)]) / 2)));
	print STDERR "Mean, median, min, max pairwise ANI differences for just kept genomes: $mean_diff, $median_diff, $min_diff, $max_diff\n";
	print STDERR "Percentiles for ANI differences for kept genomes:\n";
	my $slice = $num_diffs / 100;
	for (my $i=1; $i <= 100; $i++) {
	    my $index = int(($i * $slice) + 0.499);
	    my $percentile = 100 * (1 - $ordered_diffs[$index]);
	    print STDERR "$i:$percentile\n";
	}
    }
}
exit (0);
