#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use Cwd;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $dirname = dirname(__FILE__);
my $commandline = join (" ", @ARGV);
print STDERR "$commandline\n";

my ($mash_exec, $help, $cutoff, $redundant, $max_reps, $out, $input_file, $type_strain_path, $quit_ANI, $iterate, $increment, $stop_mean, $stop_min, $stop_reps, $fold_reduction);
my $redundant_dist;
my $cwd = getcwd;
my %genome_ids_paths; #hash for mapping genome ids to genome sketch paths
my @genome_ids; #array of genome identifiers used for labels
my @genome_paths; #array of genome sketch paths that are still active (not redundant, not filtered, not representative)
my @distances; #distances from type strain to other genomes used to determine median and other stats
my @indices; #indices of parallel arrays from 0 -> N-1
my @ordered_indices; #indices of parallel arrays from 0 -> N-1 ordered by distance from type strain smallest first
my @diffs; #differences between ordered ditances
my @ordered_diffs; #ordered differences between ordered ditances smallest first
my %hash_genome_paths; #hash to check that genome sketch paths are unique
my %hash_genome_ids; #hash to check that genome ids are unique
my $num_cur_genomes = 0;
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
my $num_total_redundant = 0;
my $total_redundant = 0;
my $min_redundant = 101;
my $max_redundant = -1;
my $median_redundant;
my $stddev_redundant;
my $sumsquared_redundant = 0;
my $num_new_reps = 1; #include the type strain as a representative genome from the beginning
my $num_total_reps = 1; #include the type strain as a representative genome from the beginning
my $max_increment_reps;
GetOptions("mash_exec|M=s"=>\$mash_exec, "out_prefix|o=s"=>\$out, "help|h|?"=>\$help, "cutoff|c=s"=>\$cutoff, "redundant|r=s"=>\$redundant, "max_reps|m=s"=>\$max_reps, "input_file|f=s"=>\$input_file, "type_strain|t=s"=>\$type_strain_path, "rep_dist|d=s"=>\$quit_ANI, "iterate|i"=>\$iterate, "increment|I=s"=>\$increment, "stop_mean|S=s"=>\$stop_mean, "stop_min|s=s"=>\$stop_min, "stop_reps|R=s"=>\$stop_reps, "fold_reduction|F=s"=>\$fold_reduction);
if ($help || !$input_file || !$mash_exec || !$type_strain_path) {
    if (!$help) {
	if (!$input_file) {
	    print STDERR "Must specify an input file of genome paths!\n\n";
	}
	if (!$mash_exec) {
	    print STDERR "Must specify a path to the MASH executable!\n\n";
	}
	if (!$type_strain_path) {
	    print STDERR "Must specify a path to the type strain genome MASH sketch!\n\n";
	}
    }
    print STDERR "Runs a type strain versus all MASH of a set of genome MASH sketches specified by the input_file(-f)\n\n";
    print STDERR "--------------------USAGE--------------------------------------\n";
    print STDERR "	-f input file of genome MASH sketch file paths, one per line (required)\n";
    print STDERR "	-o output file prefix (required if you do not want some generic silly name for your output files)\n";
    print STDERR "	-c cutoff value to use to ignore genomes below this ANI level to the type strain genome if desired\n";
    print STDERR "	-r cutoff value to use to ignore redundant genomes above this ANI level to the type strain genome or representative genomes if desired\n";
    print STDERR "	-d cutoff value to stop generating representative genomes above this ANI level to the representative genomes (default 100)\n";
    print STDERR "	-m maximum number of genomes to keep as representative for the species (default 1000)\n";
    print STDERR "	-M mash executable (required)\n";
    print STDERR "	-t genome MASH sketch file path to the type strain (required)\n";
    print STDERR "	-I maximum number of new representatives to be chosen per iteration after the first set (optional)\n";
    print STDERR "	-i iterate the generation of representative genomes\n";
    print STDERR "	-R stop iterations if the number of new representatives is less than or equal to this value\n";
    print STDERR "	-s stop iterations if the minimum ANI is greater than or equal to this value\n";
    print STDERR "	-S stop iterations if the mean ANI is greater than or equal to this value\n";
    print STDERR "	-F the desired fold reduction from total genomes to representative genomes\n";
    print STDERR "	-? or -h help\n";
		
    exit(1);
}

if (!$out) {
    $out = $cwd . "/temp_mash_out";
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
if (!(-e $type_strain_path)) {
    die ( "$type_strain_path the type strain genome MASH sketch path file does not exist!\n");
}
if (-d $type_strain_path) {
    die ( "$type_strain_path the type strain genome MASH sketch path file is a directory!\n");
}


###################################################################################################
sub process_group
# recursive processing of each group of ANI related genomes
{
    (my $group_ref, my $group_top, my $group_num, my $depth, my $reps_fh, my $redundant_fh)  = @_;
    my $group_size = @{ $group_ref };
    print STDERR "PG $group_num:$depth:$group_size\n";
    if ($group_size == 0) {
	die ("ERROR:empty group passed into process_group\n");
    }
    if ($group_size == 1) {
	print $reps_fh "$group_ref->[0]\t$group_num\n";
	print STDERR "Single $group_ref->[0]\n";
	return(1);
    }
    my $num_new_reps = 0;
    my $group_genome_median_index = int($group_size / 2);
    my $group_genome_median = $group_ref->[$group_genome_median_index];
    $group_size--; #this accounts for the median genome from the group being excluded
    $depth++;

    my $mash_out_file_prefix = $out . ".mash_dist_out_";
    my $mash_out_file = $mash_out_file_prefix . $group_top;
    my $tmp_mash_out_file = $mash_out_file . ".tmp";
    $group_num .= "_";
    my $num_redundant = 0;
    my $iteration = 0;
    my @genome_ids;
    my @indices;
    my @distances;
    my $num_kept = 0;
    my $mean_all;
    my $num_all = 0;
    my $total_all = 0;
    my $min_all = 101;
    my $max_all = -1;
    my $median_all;
    my $stddev_all;
    my $sumsquared_all = 0;
    print $reps_fh "$group_genome_median\t$group_num\n";
    $num_new_reps++;
    my $index = 0;
    my $mash_paths = "";
    my $first = 1;
    my $group_genome_median_path = $genome_paths[$genome_ids_paths{$group_genome_median}];
    foreach my $genome_id (@{ $group_ref }) {
	if ($genome_id eq $group_genome_median) {
	    next; #skip the median we are using as the query sequence
	}
	$mash_paths .= " $genome_paths[$genome_ids_paths{$genome_id}]";
	$index++;
	if ($index == 1000) {
	    print STDERR "Executing command:\n$mash_exec dist -t $group_genome_median_path mash_paths > $tmp_mash_out_file\n";
	    `$mash_exec dist -t $group_genome_median_path $mash_paths > $tmp_mash_out_file`;
	    if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
		die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
	    }
	    my $batch_size = `wc -l < $tmp_mash_out_file`;
	    $index++;
	    if ($batch_size != $index ) {
		die ( "ERROR: mash dist command failed: $tmp_mash_out_file has $batch_size lines - expecting $index lines!\n");
	    }
	    if ($first) {
		`cat $tmp_mash_out_file > $mash_out_file`;
		$first = 0;
	    } else {
		`tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
	    }
	    unlink $tmp_mash_out_file;
	    $mash_paths = "";
	    $index = 0;
	}
    }
    if ($index) {
	print STDERR "Executing command:\n$mash_exec dist -t $group_genome_median_path mash_paths > $tmp_mash_out_file\n";
	`$mash_exec dist -t $group_genome_median_path $mash_paths > $tmp_mash_out_file`;
	if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
	    die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
	}
	my $batch_size = `wc -l < $tmp_mash_out_file`;
	$index++;
	if ($batch_size != $index ) {
	    die ( "ERROR: mash dist command failed: $tmp_mash_out_file has $batch_size lines - expecting $index lines!\n");
	}
	if ($first) {
	    `cat $tmp_mash_out_file > $mash_out_file`;
	} else {
	    `tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
	}
	unlink $tmp_mash_out_file;
	$mash_paths = "";
	$index = 0;
    }
	
    my $dist_fh;
    unless (open ($dist_fh, "<", $mash_out_file) )  {
	die ("ERROR: Cannot open mash distances file $mash_out_file!\n");
    }
    my $dist = <$dist_fh>; #process header line
    my @fields = split(/\s/, $dist);
    my $num_fields = @fields;
    @genome_ids = ();
    @indices = ();
    @distances = ();
    if ($num_fields != 2) {
	die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
    }
    if ($fields[1] ne $group_genome_median) {
	die ("ERROR: query name in tablular mash output is $fields[1] expecting $group_genome_median\n");
    }
    my $row_count = 0;
    $num_redundant = 0;
    $num_kept = 0;
    $mean_all;
    $num_all = 0;
    $total_all = 0;
    $min_all = 101;
    $max_all = -1;
    $median_all;
    $stddev_all;
    $sumsquared_all = 0;
    while ($dist = <$dist_fh>) { #process the MASH lines
	if ($row_count >= $group_size) {
	    die ("ERROR: The number of distances in the tabular mash output exceeds the number of genome identifiers provided ($group_size).\n");
	}
	@fields = split(/\s/, $dist);
	$num_fields = @fields;
	if ($num_fields != 2) {
	    die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
	}
	$genome_ids[$row_count] = $fields[0];
	my $ani_est = 100 * (1 - $fields[1]);
	if ($ani_est >= $redundant) {
	    print $redundant_fh "$group_genome_median\t$genome_ids[$row_count]\t$ani_est\n";
	    $num_redundant++;
	} else {
	    $num_kept++;
	}
	$num_all++;
	$total_all += $ani_est;
	$sumsquared_all += $ani_est * $ani_est;
	if ($ani_est < $min_all) {
	    $min_all = $ani_est;
	}
	if ($ani_est > $max_all) {
	    $max_all = $ani_est;
	}
	$distances[$row_count] = $fields[1];
	$indices[$row_count] = $row_count;
	$row_count++;
    }
    if ($row_count != $group_size) {
	die ("ERROR: The number of distances in the tabular mash output ($row_count) is not the same as the number of genome identifiers provided ($group_size).\n");
    }
    close($dist_fh);
    unlink $mash_out_file;
    $num_total_redundant += $num_redundant;
    print STDERR "#reds for $group_genome_median $num_redundant:$num_total_redundant:$group_num\n";
    my @ordered_indices = sort { $distances[$a] <=> $distances[$b] || $a <=> $b } @indices; # sort indices from smallest to largest distance to genome median
    if ($num_all > 0) {
	$mean_all = $total_all / $num_all;
	$median_all = 100 * (1 - (($num_all % 2) ? $distances[$ordered_indices[($num_all / 2)]] : (($distances[$ordered_indices[(($num_all / 2) - 1)]] + $distances[$ordered_indices[($num_all / 2)]]) / 2)));
	$stddev_all = sqrt(($sumsquared_all - ($mean_all * $mean_all * $num_all)) / (($num_all > 1) ? ($num_all - 1) : 1));
    } else {
	$mean_all = 0;
	$median_all = 0;
	$stddev_all = 0;
    }
    print STDERR "Mean, median, min, max, std_dev pairwise ANI for $num_all genomes: $mean_all, $median_all, $min_all, $max_all, $stddev_all\n";
    if ($num_kept > 0) {
	my @group; #array of genome ids which might be redundant to each other
	my $prev_ordered_distance = $distances[$ordered_indices[$num_redundant]];
	my $begin_ordered_distance = $distances[$ordered_indices[$num_redundant]];
	my $group_num_suffix = 1;
	push(@group, $genome_ids[$ordered_indices[$num_redundant]]);
	for (my $i=$num_redundant + 1; $i < $group_size; $i++) {
	    my $ordered_distance = $distances[$ordered_indices[$i]];
	    my $diff_distance_begin = $ordered_distance - $begin_ordered_distance;
	    my $diff_distance_prev = $ordered_distance - $prev_ordered_distance;
	    if (($diff_distance_begin <= (2 * $redundant_dist)) && ($diff_distance_prev <= $redundant_dist)) {
		#	    if (($diff_distance_begin <= $redundant_dist) && ($diff_distance_prev <= ($redundant_dist / 2))) {
		push(@group, $genome_ids[$ordered_indices[$i]]);
	    } else {
		my $new_group_num = $group_num . $group_num_suffix;
		print STDERR "DR:$begin_ordered_distance:$prev_ordered_distance:$ordered_distance:$diff_distance_prev:$diff_distance_begin\n";
		$num_new_reps += &process_group(\@group, $group_top, $new_group_num, $depth, $reps_fh, $redundant_fh);
		print STDERR "#reps $num_new_reps:$new_group_num\n";
		$begin_ordered_distance = $ordered_distance;
		@group = ();
		$group_num_suffix++;
		push(@group, $genome_ids[$ordered_indices[$i]]);
	    }
	    $prev_ordered_distance = $ordered_distance;
	}
	if (@group > 0) {
	    my $new_group_num = $group_num . $group_num_suffix;
	    $num_new_reps += &process_group(\@group, $group_top, $new_group_num, $depth, $reps_fh, $redundant_fh);
	    print STDERR "#reps $num_new_reps:$new_group_num\n";
	}
    }

    return($num_new_reps);
}

#################################################################################
#main process

my $input_fh;
unless (open ($input_fh, "<", $input_file) )  {
    die ("ERROR: Cannot open genome sketch paths input file $input_file!\n");
}

my $mash_out_file_prefix = $out . ".mash_dist_out_";
my $mash_out_file = $mash_out_file_prefix . "0";
my $tmp_mash_out_file = $mash_out_file . ".tmp";
my $stats_file = $out . ".stats";
my $out_file = $out . ".out";
my $filtered_file = $out . ".filtered";
my $redundant_file = $out . ".redundant";
my $reps_file = $out . ".reps";
my $reps_fh;
my $index = 0;
my $mash_paths = "";
my $first = 1;
while (my $genome_path = <$input_fh>) {
    chomp ($genome_path); #remove newline
    if (defined $hash_genome_paths{$genome_path}) {
	die ("ERROR: $genome_path occurs more than once in $input_file\n");
    } else {
	$hash_genome_paths{$genome_path} = 1;
    }
    push (@genome_paths, $genome_path);
    $num_cur_genomes++;
    if (!(-e $genome_path)) {
	die ( "$genome_path the genome sketch path file does not exist!\n");
    }
    if (-d $genome_path) {
	die ( "$genome_path the genome sketch pathfile is a directory!\n");
    }
    $mash_paths .= " $genome_path";
    $index++;
    if ($index == 1000) {
	print STDERR "Executing command:\n$mash_exec dist -t $type_strain_path mash_paths > $tmp_mash_out_file\n";
	`$mash_exec dist -t $type_strain_path $mash_paths > $tmp_mash_out_file`;
	if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
	    die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
	}
	my $batch_size = `wc -l < $tmp_mash_out_file`;
	$index++;
	if ($batch_size != $index ) {
	    die ( "ERROR: mash dist command failed: $tmp_mash_out_file has $batch_size lines - expecting $index lines!\n");
	}
	if ($first) {
	    `cat $tmp_mash_out_file > $mash_out_file`;
	    $first = 0;
	} else {
	    `tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
	}
	unlink $tmp_mash_out_file;
	$mash_paths = "";
	$index = 0;
    }
}
if ($num_cur_genomes == 0) {
    die ("ERROR: There were no valid genome sketch paths in $input_file\n");
}
%hash_genome_paths = (); #free memory
close($input_fh);
if ($index) {
    print STDERR "Executing command:\n$mash_exec dist -t $type_strain_path mash_paths > $tmp_mash_out_file\n";
    `$mash_exec dist -t $type_strain_path $mash_paths > $tmp_mash_out_file`;
    if (!(-e $tmp_mash_out_file) || !(-s $tmp_mash_out_file)) {
	die ( "ERROR: mash dist command failed: $tmp_mash_out_file does not exist or is zero size!\n");
    }
    my $batch_size = `wc -l < $tmp_mash_out_file`;
    $index++;
    if ($batch_size != $index ) {
	die ( "ERROR: mash dist command failed: $tmp_mash_out_file has $batch_size lines - expecting $index lines!\n");
    }
    if ($first) {
	`cat $tmp_mash_out_file > $mash_out_file`;
    } else {
	`tail -n +2 $tmp_mash_out_file >> $mash_out_file`;
    }
    unlink $tmp_mash_out_file;
    $mash_paths = "";
    $index = 0;
}

if ($max_reps ne "") { 
    if (!looks_like_number($max_reps)) {
	die ("ERROR: the maximum number of representative genomes ($max_reps) does not look like a number\n");
    }
    if (($max_reps < 50) || ($max_reps > 1000000)) {
	die ("ERROR: $max_reps is outside of the expected maximum number of representative genomes range of 50-1,000,000\n");
    }
} else {
    $max_reps = 1000;
    if ($max_reps > $num_cur_genomes) {
	$max_reps = int(($num_cur_genomes + 1) / 2);
    }
    print STDERR "Maximum number of representative genomes set to default of $max_reps!n";
}

if ($stop_reps ne "") { 
    if (!looks_like_number($stop_reps)) {
	die ("ERROR: the minimum number of new representative genomes to keep iterating ($stop_reps) does not look like a number\n");
    }
} else {
    $stop_reps = 0;
    print STDERR "Minimum number of new representative genomes to keep iterating set to default of $stop_reps!n";
}

if ($increment ne "") { 
    if (!looks_like_number($increment)) {
	die ("ERROR: the maximum number of representative genomes to be chosen each iteration after the first ($increment) does not look like a number\n");
    }
    if (($increment < 1) || ($increment > $max_reps)) {
	die ("ERROR: $increment is outside of the expected maximum number of representative genomes per iteration range of 1-$max_reps\n");
    }
    $increment = int($increment);
}

if ($fold_reduction ne "") { 
    if (!looks_like_number($fold_reduction)) {
	die ("ERROR: the desired fold reduction of genomes ($fold_reduction) does not look like a number\n");
    }
    if (($fold_reduction < 1) || ($fold_reduction > $max_reps)) {
	die ("ERROR: $fold_reduction is outside of the expected maximum number of representative genomes per iteration range of 1-$max_reps\n");
    }
    $fold_reduction = int($fold_reduction);
} else {
    $fold_reduction = 1 + int($num_cur_genomes / $max_reps);
}

if ($cutoff ne "") { 
    if (!looks_like_number($cutoff)) {
	die ("ERROR: the ANI cutoff ($cutoff) does not look like a number\n");
    }
    if (($cutoff < 80) || ($cutoff > 100)) {
	die ("ERROR: $cutoff is outside of the expected ANI cutoff range to exclude genomesof 80-100\n");
    }
} else {
    $cutoff = 95;
    print STDERR "ANI cutoff to exclude genomes set to default of $cutoff!n";
}

if ($redundant ne "") { 
    if (!looks_like_number($redundant)) {
	die ("ERROR: the ANI redundant cutoff ($redundant) does not look like a number\n");
    }
    if (($redundant < 99) || ($redundant > 100)) {
	die ("ERROR: $redundant is outside of the expected ANI redundant cutoff range of 99-100\n");
    }
} else {
    $redundant = 99.99;
    print STDERR "ANI redundant cutoff set to default $redundant!\n";
}
$redundant_dist = (100 - $redundant) / 100;

if ($quit_ANI ne "") { 
    if (!looks_like_number($quit_ANI)) {
	die ("ERROR: the ANI representative minimum distance cutoff ($quit_ANI) does not look like a number\n");
    }
    if (($quit_ANI < 99) || ($quit_ANI > 100)) {
	die ("ERROR: $quit_ANI is outside of the expected ANI representative minimum distance cutoff range of 99-100\n");
    }
} else {
    $quit_ANI = $redundant;
    print STDERR "ANI representative minimum distance cutoff set to ANI redundant cutoff $quit_ANI!\n";
}

if ($stop_mean ne "") { 
    if (!looks_like_number($stop_mean)) {
	die ("ERROR: the ANI mean cutoff ($stop_mean) does not look like a number\n");
    }
    if (($stop_mean < $cutoff) || ($stop_mean > $redundant)) {
	die ("ERROR: $stop_mean is outside of the expected ANI mean range of $cutoff-$redundant\n");
    }
} else {
    $stop_mean = $redundant;
    print STDERR "ANI stopping mean cutoff set to ANI redundant cutoff $stop_mean!\n";
}

if ($stop_min ne "") { 
    if (!looks_like_number($stop_min)) {
	die ("ERROR: the ANI minimum cutoff ($stop_min) does not look like a number\n");
    }
    if (($stop_min < $cutoff) || ($stop_min > $redundant)) {
	die ("ERROR: $stop_min is outside of the expected ANI minimum range of $cutoff-$redundant\n");
    }
} else {
    $stop_min = $redundant;
    print STDERR "ANI stopping minimum cutoff set to ANI redundant cutoff $stop_min!\n";
}

my $dist_fh;
unless (open ($dist_fh, "<", $mash_out_file) )  {
    die ("ERROR: Cannot open mash distances file $mash_out_file!\n");
}
my $stats_fh;
unless (open ($stats_fh, ">", $stats_file) )  {
    die ("ERROR: Cannot open stats output file $stats_file!\n");
}
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
my $out_fh;
unless (open ($out_fh, ">", $out_file) )  {
    die ("ERROR: Cannot open output file $out_file!\n");
}
print $out_fh "#Type strain $type_strain_id\n";
print $stats_fh "Number of genome sketches input: $num_cur_genomes\n";
my $row_count = 0;
my $type_strain_present = 0;
while ($dist = <$dist_fh>) { #process the MASH lines
    if ($row_count >= $num_cur_genomes) {
	die ("ERROR: The number of distances in the tabular mash output exceeds the number of genome identifiers provided ($num_cur_genomes).\n");
    }
    @fields = split(/\s/, $dist);
    $num_fields = @fields;
    if ($num_fields != 2) {
	die ("ERROR: Unexpected number of fields ($num_fields) in tablular mash output - expecting 2.\n$dist\n");
    }
    $genome_ids[$row_count] = $fields[0];
    if (defined $hash_genome_ids{$genome_ids[$row_count]}) {
	die ("ERROR: $genome_ids[$row_count] occurs more than once in $mash_out_file\n");
    } else {
	$hash_genome_ids{$genome_ids[$row_count]} = 1;
    }
    $genome_ids_paths{$genome_ids[$row_count]} = $row_count;
    if ($type_strain_id eq $genome_ids[$row_count]) {
	print $stats_fh "Type strain $type_strain_id found in genome sketches input - skipping\n";
	$type_strain_present = 1;
	$distances[$row_count] = $fields[1];
	if ($distances[$row_count] > 0) {
	    die ("ERROR: mash distance for type strain to itself > 0!\n");
	}
	$indices[$row_count] = $row_count;
	$row_count++;
	next; #skip type strain if present
    }
    my $ani_est = 100 * (1 - $fields[1]);
    if ($ani_est >= $redundant) {
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
    $distances[$row_count] = $fields[1];
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
if ($row_count != $num_cur_genomes) {
    die ("ERROR: The number of distances in the tabular mash output ($row_count) is not the same as the number of genome identifiers provided ($num_cur_genomes).\n");
}
%hash_genome_ids = (); #free memory
$num_total_redundant += $num_redundant;
close($dist_fh);
unlink($mash_out_file);
@ordered_indices = sort { $distances[$a] <=> $distances[$b] || $a <=> $b } @indices; # sort indices from smallest to largest distance to type strain
my $type_plus_red = $type_strain_present + $num_redundant;
for (my $i=0; $i < $type_plus_red; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    if ($ani_est < $redundant) {
	die ("ERROR: Sorted distance/ANI values to type strain does not agree with number redundant count!\n");
    }
    if ($type_strain_id ne $genome_ids[$ordered_indices[$i]]) {
	print $redundant_fh "$type_strain_id\t$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
    }
}
my $red_plus_kept = $type_plus_red + $num_kept;
for (my $i=$type_plus_red; $i < $red_plus_kept; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    if ($ani_est >= $redundant) {
	die ("ERROR: Sorted distance/ANI values to type strain does not agree with number redundant count!\n");
    }
    print $out_fh "$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
}
my $type_plus_all = $type_strain_present + $num_all;
if ($type_plus_all != $num_cur_genomes) {
    die ("ERROR: The number of genomes in the input including the type strain ($type_plus_all) is not equal to the number of input genome sketches ($num_cur_genomes)!\n");
}
for (my $i=$red_plus_kept; $i < $type_plus_all; $i++) {
    my $ani_est = 100 * (1 - $distances[$ordered_indices[$i]]);
    print $filtered_fh "$genome_ids[$ordered_indices[$i]]\t$ani_est\n";
}
close($out_fh);
if ($cutoff ne "") {
    close($filtered_fh);
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
print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for all $num_all genomes: $mean_all, $median_all, $min_all, $max_all, $stddev_all\n";
if ($cutoff ne "") {
    print $stats_fh "For cutoff ($cutoff) and type strain $type_strain_id: $num_discard discarded\n";
}
print $stats_fh "For redundant cutoff ($redundant) and type strain $type_strain_id: $num_redundant redundant\n";
print $stats_fh "For type strain $type_strain_id: $num_kept kept\n";
if ($num_redundant > 0) {
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
if (((($cutoff ne "") && ($num_discard > 0)) || ($num_redundant > 0)) && ($num_kept > 0)) {
    $mean_kept = $total_kept / $num_kept;
    $median_kept = 100 * (1 - (($num_kept % 2) ? $distances[$ordered_indices[($num_kept / 2)]] : (($distances[$ordered_indices[(($num_kept / 2) - 1)]] + $distances[$ordered_indices[($num_kept / 2)]]) / 2)));
    $stddev_kept = sqrt(($sumsquared_kept - ($mean_kept * $mean_kept * $num_kept)) / (($num_kept > 1) ? ($num_kept - 1) : 1));
    print $stats_fh "Mean, median, min, max, std_dev pairwise ANI for just kept genomes: $mean_kept, $median_kept, $min_kept, $max_kept, $stddev_kept\n";
} else {
    $mean_kept = $mean_all;
    $median_kept = $median_all;
    $stddev_kept = $stddev_all;
}
print $stats_fh "Number new representatives $num_new_reps ($num_total_reps)\nNumber new redundant $num_redundant ($num_total_redundant)\n";

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
	print $stats_fh "$i\t$percentile\t$genome_ids[$ordered_indices[$index]]\n";
    }
    $index = 0;
    for (my $i=$type_plus_red + $fold_reduction; $i < $red_plus_kept; $i++) {
	$diffs[$index] = $distances[$ordered_indices[$i]] - $distances[$ordered_indices[$i-$fold_reduction]];
	$index++;
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
	    print $stats_fh "$i\t$percentile\n";
	}
    }
    if ($num_kept > $max_reps) {
	unless (open ($reps_fh, ">", $reps_file) )  {
	    die ("ERROR: Cannot open reprentative genomes output file $reps_file!\n");
	}
	my $tmp_reps = int($max_reps / 20);
	if ($tmp_reps >= ($num_kept / 20)) {
	    $tmp_reps = int($num_kept / 20);
	}
	if ($tmp_reps < 1) {
	    $tmp_reps = 1;
	}
	if ($increment eq "") { 
	    $max_increment_reps = $tmp_reps;
	} else {
	    $max_increment_reps = $increment;
	}
	print $reps_fh "$type_strain_id\t0\n";

	my @group; #array of genome ids which might be redundant to each other
	my $prev_ordered_distance = $distances[$ordered_indices[$type_plus_red]];
	my $begin_ordered_distance = $distances[$ordered_indices[$type_plus_red]];
	my $group_num = 1;
	push(@group, $genome_ids[$ordered_indices[$type_plus_red]]);
	for (my $i=$type_plus_red + 1; $i < $red_plus_kept; $i++) {
	    my $ordered_distance = $distances[$ordered_indices[$i]];
	    my $diff_distance_begin = $ordered_distance - $begin_ordered_distance;
	    my $diff_distance_prev = $ordered_distance - $prev_ordered_distance;
	    if (($diff_distance_begin <= (2 * $redundant_dist)) && ($diff_distance_prev <= $redundant_dist)) {
#	    if (($diff_distance_begin <= $redundant_dist) && ($diff_distance_prev <= ($redundant_dist / 2))) {
		push(@group, $genome_ids[$ordered_indices[$i]]);
	    } else {
		print STDERR "DR:$begin_ordered_distance:$prev_ordered_distance:$ordered_distance:$diff_distance_prev:$diff_distance_begin\n";
		$num_total_reps += &process_group(\@group, $group_num, $group_num, 1, $reps_fh, $redundant_fh);
		print STDERR "#reps $num_total_reps:$group_num\n";
		$begin_ordered_distance = $ordered_distance;
		@group = ();
		$group_num++;
		push(@group, $genome_ids[$ordered_indices[$i]]);
	    }
	    $prev_ordered_distance = $ordered_distance;
	}
	if (@group > 0) {
	    $num_total_reps += &process_group(\@group, $group_num, $group_num, 1, $reps_fh, $redundant_fh);
	    print STDERR "#reps $num_total_reps:$group_num\n";
	}
	print $stats_fh "Number representative genomes: $num_total_reps\nNumber redundant genomes: $num_total_redundant\n";
#	print $stats_fh "Number current active genomes: $num_cur_genomes\n";
    } else {
	print STDERR "Number of kept genomes $num_kept <= maximum number of genome representatives specified $max_reps - use kept genomes as representatives.\n";
    }
}
close($redundant_fh);
close($reps_fh);

exit (0);
