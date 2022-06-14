#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $dirname = dirname(__FILE__);

my ($mash_exec, $help, $cutoff, $redundant, $max_reps, $out, $input_file, $type_strain);
my $cwd = getcwd;
my @genome_ids; #array of genome identifiers used for labels
my @genome_paths; #array of genome sketch paths
my @print_ids; #array parallel to @genome_ids which indicates the genome should be printed (1) or skipped (0) based on ANI cutoff to the type strain.
my @distances; #distances from type strain to other genomes used to determine median and other stats
my @indices; #indices of parallel arrays from 0 -> N-1
my @ordered_indices; #indices of parallel arrays from 0 -> N-1 ordered by distance from type strain smallest first
my @diffs; #differences between ordered ditances
my @ordered_diffs; #ordered differences between ordered ditances smallest first
my $type_strain_id;
my $mean_all;
my $num_all = 0;
my $total_all = 0;
my $min_all = 101;
my $max_all = -1;
my $median_all;
my $stddev_all;
my $sumsquared_all = 0;
my $mean_kept;
my $num_kept = 0;
my $total_kept = 0;
my $min_kept = 101;
my $max_kept = -1;
my $median_kept;
my $stddev_kept;
my $sumsquared_kept = 0;
my $mean_discard;
my $num_discard = 0;
my $total_discard = 0;
my $min_discard = 101;
my $max_discard = -1;
my $median_discard;
my $stddev_discard;
my $sumsquared_discard = 0;
my $mean_redundant;
my $num_redundant = 0;
my $total_redundant = 0;
my $min_redundant = 101;
my $max_redundant = -1;
my $median_redundant;
my $stddev_redundant;
my $sumsquared_redundant = 0;
$cutoff = "";
GetOptions("mash_exec|M=s"=>\$mash_exec, "out_prefix|o=s"=>\$out, "help|h|?"=>\$help, "cutoff|c=s"=>\$cutoff, "redundant|r=s"=>\$redundant, "max_reps|m=s"=>\$max_reps, "input_file|f=s"=>\$input_file, "type_strain|t=s"=>\$type_strain);
if ($help || !$input_file || !$mash_exec || !$type_strain) {
    if (!$help) {
	if (!$input_file) {
	    print STDERR "Must specify an input file of genome paths!\n\n";
	}
	if (!$mash_exec) {
	    print STDERR "Must specify a path to the MASH executable!\n\n";
	}
	if (!$type_strain) {
	    print STDERR "Must specify a path to the type strain genome MASH sketch!\n\n";
	}
    }
    print STDERR "Runs a type strain versus all MASH of a set of genome MASH sketches specified by the input_file(-f)\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome MASH sketch file paths, one per line (required)\n";
    print STDERR "	-o output file prefix (required if you do not want some generic silly name for your output files)\n";
    print STDERR "	-c cutoff value to use to ignore genomes below this ANI level to the type strain genome if desired\n";
    print STDERR "	-r cutoff value to use to ignore redundant genomes above this ANI level to the type strain genome if desired\n";
    print STDERR "	-m maximum number of genomes to keep as representative for the species (default 1000)\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-t genome MASH sketch file path to the type strain (required)\n";
    print STDERR "	-? or -h help\n";
		
    exit(1);
}

if (!$out) {
    $out = $cwd . "/temp_mash_out";
}

if ($max_reps ne "") { 
    if (!looks_like_number($max_reps)) {
	die ("ERROR: the maximum number of representative genomes $max_reps does not look like a number\n");
    }
    if (($max_reps < 100) || ($max_reps > 1000000)) {
	die ("ERROR: $max_reps is outside of the expected maximum number of representative genomes range of 100-1,000,000\n");
    }
} else {
    $max_reps = 1000;
    print STDERR "Maximum number of representative genomes set to default of $max_reps!n";
}

if ($cutoff ne "") { 
    if (!looks_like_number($cutoff)) {
	die ("ERROR: the ANI cutoff $cutoff does not look like a number\n");
    }
    if (($cutoff < 0) || ($cutoff > 100)) {
	die ("ERROR: $cutoff is outside of the expected ANI cutoff range of 0-100\n");
    }
}

if ($redundant ne "") { 
    if (!looks_like_number($redundant)) {
	die ("ERROR: the ANI redundant cutoff $redundant does not look like a number\n");
    }
    if (($redundant < 90) || ($redundant > 100)) {
	die ("ERROR: $redundant is outside of the expected ANI redundant cutoff range of 90-100\n");
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
my $genome_paths_size = `wc -l < $input_file`;
chomp($genome_paths_size);

my $input_fh;
unless (open ($input_fh, "<", $input_file) )  {
    die ("ERROR: Cannot open genome sketch paths input file $input_file!\n");
}

my $mash_out_file = $out . ".mash_dist_out";
`rm $mash_out_file`;
my $tmp_mash_out_file = $out . ".tmp_mash_dist_out";
my $stats_file = $out . ".stats";
my $out_file = $out . ".out";
my $filtered_file = $out . ".filtered";
my $redundant_file = $out . ".redundant";
my $index = 0;
my $mash_paths = "";
my $first = 1;
while (my $genome_path = <$input_fh>) {
    chomp ($genome_path); #remove newline
    push (@genome_paths, $genome_path);
    if (!(-e $genome_path)) {
	die ( "$genome_path the genome sketch path file does not exist!\n");
    }
    if (-d $genome_path) {
	die ( "$genome_path the genome sketch pathfile is a directory!\n");
    }
    $mash_paths .= " $genome_path";
    $index++;
    if ($index == 1000) {
	$index = 0;
	print STDERR "Executing command:\n$mash_exec dist -t $type_strain $mash_paths > $tmp_mash_out_file\n";
	`$mash_exec dist -t $type_strain $mash_paths > $tmp_mash_out_file`;
	if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
	    die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
	}
	if ($first) {
	    `cat $tmp_mash_out_file >> $mash_out_file`;
	    $first = 0;
	} else {
	    `tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
	}
	`rm $tmp_mash_out_file`;
	$mash_paths = "";
    }
}
if ($index != 1000) {
    $index = 0;
    print STDERR "Executing command:\n$mash_exec dist -t $type_strain $mash_paths > $tmp_mash_out_file\n";
    `$mash_exec dist -t $type_strain $mash_paths > $tmp_mash_out_file`;
    if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
	die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
    }
    if ($first) {
	`cat $tmp_mash_out_file >> $mash_out_file`;
    } else {
	`tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
    }
    `rm $tmp_mash_out_file`;
    $mash_paths = "";
}

my $dist_fh;
unless (open ($dist_fh, "<", $mash_out_file) )  {
    die ("ERROR: Cannot open mash distances file $mash_out_file!\n");
}
my $stats_fh;
unless (open ($stats_fh, ">", $stats_file) )  {
    die ("ERROR: Cannot open stats output file $stats_file!\n");
}
my $out_fh;
unless (open ($out_fh, ">", $out_file) )  {
    die ("ERROR: Cannot open output file $out_file!\n");
}
print $out_fh "#Type strain $type_strain_id\n";
my $filtered_fh;
if ($cutoff ne "") {
    unless (open ($filtered_fh, ">", $filtered_file) )  {
	die ("ERROR: Cannot open filtered output file $filtered_file!\n");
    }
    print $filtered_fh "#Genomes being filtered out for ANI below cutoff ($cutoff)\n";
}
my $redundant_fh;
if ($redundant ne "") {
    unless (open ($redundant_fh, ">", $redundant_file) )  {
	die ("ERROR: Cannot open filtered output file $redundant_file!\n");
    }
    print $redundant_fh "#Genomes being filtered out for ANI above redundant cutoff ($redundant)\n";
}
my $dist = <$dist_fh>; #process header line
my @fields = split(/\s/, $dist);
my $num_fields = @fields;
if ($num_fields != 2) {
    die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
}
$type_strain_id = $fields[1];
my $row_count = 0;
while ($dist = <$dist_fh>) { #process the MASH lines
    if ($row_count >= $genome_paths_size) {
	die ("ERROR: The number of distances in the tabular mash output exceeds the number of genome identifiers provided ($genome_paths_size).\n");
    }
    @fields = split(/\s/, $dist);
    $num_fields = @fields;
    if ($num_fields != 2) {
	die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
    }
    $genome_ids[$row_count] = $fields[0];
    my $ani_est = 100 * (1 - $fields[1]);
    if (($redundant ne "") && ($ani_est >= $redundant)) {
	$print_ids[$row_count] = 0;
	$num_redundant++;
	$total_redundant += $ani_est;
	$sumsquared_redundant += $ani_est * $ani_est;
	if ($ani_est < $min_redundant) {
	    $min_redundant = $ani_est;
	}
	if ($ani_est > $max_redundant) {
	    $max_redundant = $ani_est;
	}
    } elsif (($cutoff eq "") || ($ani_est >= $cutoff)) {
	$print_ids[$row_count] = 1;
	$num_kept++;
	$total_kept += $ani_est;
	$sumsquared_kept += $ani_est * $ani_est;
	if ($ani_est < $min_kept) {
	    $min_kept = $ani_est;
	}
	if ($ani_est > $max_kept) {
	    $max_kept = $ani_est;
	}
    } else {
	$print_ids[$row_count] = 0;
	$num_discard++;
	$total_discard += $ani_est;
	$sumsquared_discard += $ani_est * $ani_est;
	if ($ani_est < $min_discard) {
	    $min_discard = $ani_est;
	}
	if ($ani_est > $max_discard) {
	    $max_discard = $ani_est;
	}
    }
    $distances[$row_count] += $fields[1];
    $indices[$row_count] = $row_count;
    $num_all++;
    $total_all += $ani_est;
    $sumsquared_all += $ani_est * $ani_est;
    if ($ani_est < $min_all) {
	$min_all = $ani_est;
    }
    if ($ani_est > $max_all) {
	$max_all = $ani_est;
    }
    $row_count++;
}
close($dist_fh);
`rm $mash_out_file`;
@ordered_indices = sort { $distances[$a] <=> $distances[$b] || $a <=> $b } @indices; # sort indices from smallest to largest distance to type strain
for (my $i=0; $i < $num_redundant; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    print $redundant_fh "$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
}
my $red_plus_kept = $num_redundant + $num_kept;
for (my $i=$num_redundant; $i < $red_plus_kept; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    print $out_fh "$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
}
for (my $i=$red_plus_kept; $i < $num_all; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    print $filtered_fh "$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
}
close($out_fh);
if ($cutoff ne "") {
    close($filtered_fh);
}
if ($redundant ne "") {
    close($redundant_fh);
}
if ($num_all > 0) {
    $mean_all = $total_all / $num_all;
    $median_all = 100 * (1 - (($num_all % 2) ? $distances[$ordered_indices[($num_all / 2)]] : (($distances[$ordered_indices[(($num_all / 2) - 1)]] + $distances[$ordered_indices[($num_all / 2)]]) / 2)));
    $stddev_all = sqrt(($sumsquared_all - ($mean_all * $mean_all * $num_all)) / (($num_all > 1) ? ($num_all - 1) : 1));
} else {
    $mean_all = 0;
    $median_all = 0;
    $stddev_all = 0;
}
print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for all $genome_paths_size genomes: $mean_all, $median_all, $min_all, $max_all, $stddev_all\n";
if (($cutoff ne "") || ($redundant ne "")) {
    if ($cutoff ne "") {
	print $stats_fh "For cutoff ($cutoff) and type strain $type_strain_id: $num_discard discarded\n";
    }
    if ($redundant ne "") {
	print $stats_fh "For redundant cutoff ($redundant) and type strain $type_strain_id: $num_redundant redundant\n";
    }
    print $stats_fh "For type strain $type_strain_id: $num_kept kept\n";
}
if (($redundant ne "") && ($num_redundant > 0)) {
    $mean_redundant = $total_redundant / $num_redundant;
    $median_redundant = 100 * (1 - (($num_redundant % 2) ? $distances[$ordered_indices[($num_redundant / 2)]] : (($distances[$ordered_indices[(($num_redundant / 2) - 1)]] + $distances[$ordered_indices[($num_redundant / 2)]]) / 2)));
    $stddev_redundant = sqrt(($sumsquared_redundant - ($mean_redundant * $mean_redundant * $num_redundant)) / (($num_redundant > 1) ? ($num_redundant - 1) : 1));
    print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for just redundant genomes: $mean_redundant, $median_redundant, $min_redundant, $max_redundant, $stddev_redundant\n";
}
if (($cutoff ne "") && ($num_discard > 0)) {
    $mean_discard = $total_discard / $num_discard;
    $median_discard = 100 * (1 - (($num_discard % 2) ? $distances[$ordered_indices[$red_plus_kept + ($num_discard / 2)]] : (($distances[$ordered_indices[$red_plus_kept + (($num_discard / 2) - 1)]] + $distances[$ordered_indices[$red_plus_kept + ($num_discard / 2)]]) / 2)));
    $stddev_discard = sqrt(($sumsquared_discard - ($mean_discard * $mean_discard * $num_discard)) / (($num_discard > 1) ? ($num_discard - 1) : 1));
    print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for just discarded genomes: $mean_discard, $median_discard, $min_discard, $max_discard, $stddev_discard\n";
}
if (((($cutoff ne "") && ($num_discard > 0)) || (($redundant ne "") && ($num_redundant > 0))) && ($num_kept > 0)) {
    $mean_kept = $total_kept / $num_kept;
    $median_kept = 100 * (1 - (($num_kept % 2) ? $distances[$ordered_indices[($num_kept / 2)]] : (($distances[$ordered_indices[(($num_kept / 2) - 1)]] + $distances[$ordered_indices[($num_kept / 2)]]) / 2)));
    $stddev_kept = sqrt(($sumsquared_kept - ($mean_kept * $mean_kept * $num_kept)) / (($num_kept > 1) ? ($num_kept - 1) : 1));
    print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for just kept genomes: $mean_kept, $median_kept, $min_kept, $max_kept, $stddev_kept\n";
} else {
    $mean_kept = $mean_all;
    $median_kept = $median_all;
    $stddev_kept = $stddev_all;
}

if ($num_kept > 0) {
    print $stats_fh "Percentiles for ANI to the type strain for kept genomes:\n";
    my $slice = $num_kept / 100;
    for (my $i=1; $i <= 100; $i++) {
	$index = int(($i * $slice) + 0.499) - 1;
	if ($index < 0) {
	    $index = 0;
	} elsif ($index >= $num_kept) {
	    $index = $num_kept - 1;
	}
	my $percentile = 100 * (1 - $distances[$ordered_indices[$index]]);
	print $stats_fh "$i:$percentile\t$genome_ids[$ordered_indices[$index]]\n";
    }
    if ($num_kept > $max_reps) {
	my $reps_file = $out . ".reps";
	my $reps_sk_file = $out . ".reps_msh";
	my $reps_fh;
	unless (open ($reps_fh, ">", $reps_file) )  {
	    die ("ERROR: Cannot open stats output file $reps_file!\n");
	}
	my $reps_sk_fh;
	unless (open ($reps_sk_fh, ">", $reps_sk_file) )  {
	    die ("ERROR: Cannot open stats output file $reps_sk_file!\n");
	}
	my $tmp_reps = int($max_reps / 10);
	my $rep_slice = $num_kept / $tmp_reps;
	for (my $i=1; $i <= $tmp_reps; $i++) {
	    $index = int(($i * $slice) + 0.499) - 1;
	    if ($index < 0) {
		$index = 0;
	    } elsif ($index >= $num_kept) {
		$index = $num_kept - 1;
	    }
	    my $rep_ANI = 100 * (1 - $distances[$ordered_indices[$index]]);
	    print $reps_fh "$genome_ids[$ordered_indices[$index]]/t$rep_ANI\n";
	    print $reps_sk_fh "$genome_paths[$ordered_indices[$index]]\n";
	}
	close($reps_fh);
	close($reps_sk_fh);
    } else {
	print STDERR "Number of kept genomes $num_kept <= maximum number of genome representatives specified $max_reps - use kept genomes as representatives.\n";
    }
    for (my $i=1; $i < $num_kept; $i++) {
	$diffs[$i-1] = $distances[$ordered_indices[$i]] - $distances[$ordered_indices[$i-1]];
    }
    my $num_diffs = @diffs;
    if ($num_diffs > 0) {
	my $total_diff = 0;
	my $min_diff = 100;
	my $max_diff = 0;
	my $mean_diff;
	my $median_diff;
	foreach my $diff (@diffs) {
	    my $ani_est_diff = 100 * $diff;
	    $total_diff += $ani_est_diff;
	    if ($ani_est_diff < $min_diff) {
		$min_diff = $ani_est_diff;
	    }
	    if ($ani_est_diff > $max_diff) {
		$max_diff = $ani_est_diff;
	    }
	}
	@ordered_diffs = sort { $a <=> $b } @diffs; # sort diffs from smallest to largest
	$mean_diff = $total_diff / $num_diffs;
	$median_diff = 100 * ((($num_diffs % 2) ? $ordered_diffs[($num_diffs / 2)] : (($ordered_diffs[(($num_diffs / 2) - 1)] + $ordered_diffs[($num_diffs / 2)]) / 2)));
	print $stats_fh "Mean, median, min, max pairwise ANI differences for just kept genomes: $mean_diff, $median_diff, $min_diff, $max_diff\n";
	print $stats_fh "Percentiles for ANI differences for kept genomes:\n";
	my $slice = $num_diffs / 100;
	for (my $i=1; $i <= 100; $i++) {
	    $index = int(($i * $slice) + 0.499) - 1;
	    if ($index < 0) {
		$index = 0;
	    } elsif ($index >= $num_diffs) {
		$index = $num_diffs - 1;
	    }
	    my $percentile = 100 * $ordered_diffs[$index];
	    print $stats_fh "$i:$percentile\n";
	}
    }
}
close($stats_fh);
exit (0);
