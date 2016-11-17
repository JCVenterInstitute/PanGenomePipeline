#!/usr/bin/env perl
#Copyright (C) 2014-2015  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Granger Sutton, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use feature 'state';
getopts ('PDhW:M:m:p:l:C:I:i:A:a:g:t:T:L:V');# M is _adjacency_vector.txt from PanOCT
our ($opt_a, $opt_g, $opt_h, $opt_i, $opt_l, $opt_m, $opt_p, $opt_t, $opt_A, $opt_C, $opt_D, $opt_I, $opt_L, $opt_M, $opt_P, $opt_T, $opt_V, $opt_W);

my $high_adjvec_file_name;
my $low_adjvec_file_name;
my $cluster_file_name;
my $centroids_file_name;
my $attribute_file_name;
my $aux_attribute_file_name;
my $genome_tags_file_name;
my $genome_labels_file_name;
my $insert_file_name;
my $input_insert_file_name;
my $specified = 0;
my $genes_file_name;
my %adj_vec = (); # adjacency vector read in from the input
my %high_adj_mat = (); # adjacency matrix read in from the input
my %high_adj_best = (); # best (maximum weight) edge
my %high_adj_mutual = (); # edges which are mutually best
my %high_chain_length = (); # Key = cluster number, Value = length of chain starting at this cluster - only defined for start points
my %high_chain_gene_length = (); # Key = cluster number, Value = length of chain starting at this cluster - only defined for start points
my @high_clus_start = ();
my @high_clus_present = ();
my %low_adj_mat = (); # adjacency matrix read in from the input
my %low_adj_best = (); # best (maximum weight) edge
my %low_adj_mutual = (); # edges which are mutually best
my %low_chain_length = (); # Key = cluster number, Value = length of chain starting at this cluster - only defined for start points
my %low_chain_gene_length = (); # Key = cluster number, Value = length of chain starting at this cluster - only defined for start points
my %max_chain_node = ();
my %max_chain_edge = ();
my @low_clus_start = ();
my @low_clus_present = ();
my @centroid_length = ();
my @centroid_anno = ();
my @context_count = ();
my @cluster_size = ();
my @assembly_num_clus = ();
my @cluster_strand = ();
my @genome_tags = ();
my @groups = ();
my %tag_group = ();
my %global_group_counts = ();
my %end_coord = ();
my $output_centroids;
my $high_global_largest_edge = 0;
my $low_global_largest_edge = 0;
my $DEBUG = 0;
my $comma_tags;
my $cluster_prefix;
my $num_clusters;
my $context_length;
my $scale_protein;
my $combine_high_low;
my $global_assembly_num = 1;
my $num_genomes;
my $version = "1.3";

if ($opt_D) {
  $DEBUG = 1;
} else {
  $DEBUG = 0;
} # Debug mode is off as default.
if ($opt_T) {
  $comma_tags = $opt_T;
} else {
  $comma_tags = 0;
} # stye 0 is default.
if ($opt_P) {
  $scale_protein = 3;
} else {
  $scale_protein = 1;
} # Scale protein lengths to nucleotide in attribute file
if ($opt_h) { 
  &option_help;
} # quit with help menu
if ($opt_V) {die "$prog version $version\n";}
if ($opt_W) {
  if (-s "$opt_W") {
    $cluster_file_name = $opt_W;
  } else {
    print STDERR "Error $opt_W is not an existing file\n";
    &option_help;
  }
} else {
  print STDERR "Error must specify -W\n";
  &option_help;
}
if ($opt_M) {
  if (-s "$opt_M") {
    $high_adjvec_file_name = $opt_M;
  } else {
    print STDERR "Error $opt_M is not an existing file\n";
    &option_help;
  }
} else {
  print STDERR "Error must specify -M\n";
  &option_help;
}
if ($opt_m) {
  if (-s "$opt_m") {
    $low_adjvec_file_name = $opt_m;
  } else {
    print STDERR "Error $opt_m is not an existing file\n";
    &option_help;
  }
  $combine_high_low = 1;
  if ($opt_a) {
    $aux_attribute_file_name = $opt_a;
  } else {
    print STDERR "Error must specify -a with -m\n";
    &option_help;
  }
  if ($opt_t) {
    if (-s "$opt_t") {
      $genome_tags_file_name = $opt_t;
    } else {
      print STDERR "Error $opt_t is not an existing file\n";
      &option_help;
    }
  } else {
    print STDERR "Error must specify -t with -m\n";
    &option_help;
  }
  if ($opt_L) {
    if (-s "$opt_L") {
      $genome_labels_file_name = $opt_L;
    } else {
      print STDERR "Error $opt_L is not an existing file\n";
      &option_help;
    }
  } else {
    if ($comma_tags eq "3") {
      print STDERR "Error must specify -L with -m if -T 3 specified\n";
      &option_help;
    }
  }
  if ($opt_i) {
    if (-s "$opt_i") {
      $input_insert_file_name = $opt_i;
      $specified = 1;
    } else {
      print STDERR "Error $opt_i is not an existing file\n";
      &option_help;
    }
  }
  if ($opt_I) {
    $insert_file_name = $opt_I;
  } else {
    print STDERR "Error must specify -I with -m\n";
    &option_help;
  }
} else {
  $combine_high_low = 0;
}
if ((!$opt_m) && $opt_a) {
  print STDERR "Error cannot have -a without -m\n";
  &option_help;
}
if ((!$opt_m) && $opt_I) {
  print STDERR "Error cannot have -I without -m\n";
  &option_help;
}
if ((!$opt_m) && $opt_i) {
  print STDERR "Error cannot have -i without -m\n";
  &option_help;
}
if (($opt_C) && (-s "$opt_C")) {
  $centroids_file_name = $opt_C;
} else {
  print STDERR "Error must specify -C\n";
  &option_help;
}
if ($opt_A) {
  $attribute_file_name = $opt_A;
} else {
  print STDERR "Error must specify -A\n";
  &option_help;
}
if ($opt_g) {
  $genes_file_name = $opt_g;
  $output_centroids = 1;
} else {
  $output_centroids = 0;
}
if ($opt_p) {
  $cluster_prefix = $opt_p;
} else {
  $cluster_prefix = "CL";
}
if ($opt_l) {
  if (!(looks_like_number($opt_l))) {
    die ("ERROR: $opt_l is not a number for length of context (-l)!\n");
  }
  if ($opt_l <= 0) {
    die ("ERROR: $opt_l is not > 0 for length of context (-l)!\n");
  }
  $context_length = int($opt_l);
} else {
  $context_length = 5;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog:

           Takes an adjacency_vector.txt file from PanOCT and determines a consensus gene order file. 
           The gene order file is sent to standard out as specified prefix concatenated onto cluster 
           number separated by a tab from the strand +/-. A blank line indicates a break in the gene 
           order as for multiple contigs. A genes primary context is taken from it's first appearance 
           in the file but can show up later in the file for other genes' context. A PanOCT centroids 
           file must be input using -C and a gene attribute file name must be specified using -A and 
           will be output using the gene order and centroid information.

Copyright (C) 2014-2015  The J. Craig Venter Institute (JCVI).  All rights reserved

License:   This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

Citation:  Agnes P. Chan, Granger Sutton, Jessica DePew, Radha Krishnakumar, Yongwook Choi, Xiao-Zhe Huang, 
           Erin Beck, Derek M. Harkins, Maria Kim, Emil P. Lesho, Mikeljon P. Nikolich and Derrick E. Fouts 
           (2015) "A novel method of consensus pan-chromosome assembly and large-scale comparative analysis 
           reveal the highly flexible pan-genome of Acinetobacter baumannii" Genome Biol. 16(143):1-28.

  Usage: $prog <options>
Example: $prog -M 0_core_adjacency_vector.txt -W cluster_weights.txt -p L1B1C -l 5 -A L1B1.att -g L1B1.pep > gene_order.txt # for hierarchical/iterative pan-genome runs
Example: $prog -P -W cluster_weights.txt -M 75_core_adjacency_vector.txt -m 0_core_adjacency_vector.txt -t ../example_tags.txt -l 5 -A Core.att -a fGI.att -I fGI_report.txt -C centroids.fasta > gene_order.txt # fGI analysis
Version: $version
Options:
     -h: print this help page
     -W: cluster weights file from PanOCT output, tab delimited, column 1 is cluster number, column 2 is cluster size [REQUIRED]
     -M: num_core_adjacency_vector.txt file from PanOCT, num should be 75 for consensus/fGI analysis or 0 for hierarchical PanOCT runs [REQUIRED]
     -m: num_core_adjacency_vector.txt file from PanOCT, num should be 0 for consensus/fGI analysis, optional (not used for hierarchical PanOCT runs)
     -t: genome tags file name used as input for PanOCT, required for -m and ignored otherwise
     -L: genome groups label file name: first column is genome tag second column is group label used with -m and -T 3
     -T: controls how genome vectors are output: 0 vectors of 0s and 1s tab delimited, 1 genome tags comma delimited
         2 also include for edges in consensus file, 3 use group labels from -L file
     -C: centroids.fasta file from PanOCT [REQUIRED]
     -I: output file name for a report on fGI insertions into the Core pan-genome, optional (not used for hierarchical PanOCT runs)
     -A: output file name for a PanOCT gene attribute file for either the Core pan-genome or entire genome for hierarchical PanOCT runs [REQUIRED]
     -a: output file name for a PanOCT gene attribute file for flexible Genomic Islands (fGI), optional (not used for hierarchical PanOCT runs)
     -g: output file name for a PanOCT genes/proteins file [for iterative pan-genome only]
     -p: prefix for clusters in output, optional (default CL)
     -l: context length for output, optional (default 5)
     -P: no argument - specifies protein lengths should be scaled to nucleotide values in attribute files
     -D: DEBUG MODE (default = off)
 Output: consensus gene order for clusters sent to STDOUT, cluster type P is primary, type S is for context only
         a PanOCT gene attribute file for the consensus gene order (-A and -a)
	 a report on fGI insertions (-I)
	 a multifasta peptide file of centroids for hierarchical PanOCT runs (-g)
     -V: print just (V)ersion information
 Author: Granger Sutton, Ph.D.
 Date: 09/01/2015
_EOB_
    exit;
}

sub get_cluster_sizes {

    my $number = "";
    my $size = "";
    my $cluster_file;
    $num_genomes = 0;

    unless (open ($cluster_file, $cluster_file_name) )  {
	die ("ERROR: can not open file $cluster_file_name.\n");
    }
    $num_clusters = 0;
    while (<$cluster_file>) {
	my @cluster_line = ();
	chomp;
	@cluster_line = split(/\t/, $_);  # split the scalar $cluster_line on tab
	$number = $cluster_line[0];
	if ($number eq "") {
	    die ("ERROR: cluster number must not be empty/null in cluster file\n");
	}
	if (!(looks_like_number($number))) {
	    die ("ERROR: $number is not a number for cluster number in cluster file!\n");
	}
	if ($number <= 0) {
	    die ("ERROR: $number is not > 0 for cluster number in cluster file!\n");
	}
	$size = $cluster_line[1];
	if (!(looks_like_number($size))) {
	    die ("ERROR: $size is not a number for cluster size in cluster file!\n");
	}
	if ($size <= 0) {
	    die ("ERROR: $size is not > 0 for cluster size in cluster file!\n");
	}
	if ($size eq "") {
	    die ("ERROR: cluster size must not be empty/null in cluster file\n");
	}
	$cluster_size[$number] = $size;
	if ($size > $num_genomes) {
	    $num_genomes = $size;
	}
	$num_clusters++;
	if ($number != $num_clusters) {
	    die ("ERROR: cluster $num_clusters was not present in $cluster_file_name - skipped to $number\n");
	}
    }
    for my $clus (1 .. $num_clusters) {
	$assembly_num_clus[$clus] = 0; #if cluster is not present then have assembly number for cluster be 0
    }
    close ($cluster_file);
    return;
}

sub get_genome_tags {

    my $genome_tags_file;
    my %unique_genome_tags = ();

    unless (open ($genome_tags_file, $genome_tags_file_name) )  {
	die ("ERROR: can not open file $genome_tags_file_name.\n");
    }
    while (<$genome_tags_file>) {
	chomp;
	(my $genome_tag, my $group) =  split(/\t/, $_);  # split on tab
	push(@genome_tags, $genome_tag);
	$unique_genome_tags{$genome_tag} = $genome_tag;
    }
    close ($genome_tags_file);
    my $genome_tags_size = @genome_tags;
    if ($genome_tags_size != $num_genomes) {
	die "ERROR: mismatch between number of genome tags ($genome_tags_size) and number of genomes in clusters ($num_genomes) from file $cluster_file_name\n";
    }

    if ($comma_tags eq "3") {
	my $genome_labels_file;
	
	unless (open ($genome_labels_file, $genome_labels_file_name) )  {
	    die ("ERROR: can not open file $genome_labels_file_name.\n");
	}
	while (<$genome_labels_file>) {
	    chomp;
	    (my $genome_tag, my $group) =  split(/\t/, $_);  # split on tab
	    if (!defined $unique_genome_tags{$genome_tag}) {
		die ("ERROR: genome tag ( $genome_tag ) was present in the group labels file $genome_labels_file_name but not in the genome tags file $genome_tags_file_name\n");
	    } elsif ($unique_genome_tags{$genome_tag} eq $genome_tag) {
		$unique_genome_tags{$genome_tag} = "";
	    } else {
		die ("ERROR: genome tag ( $genome_tag ) was present in the group labels file $genome_labels_file_name more than once\n");
	    }
	    if (defined $group) {
		$tag_group{$genome_tag} = $group;
		if (!defined $global_group_counts{$group}) {
		    $global_group_counts{$group} = 1;
		    push(@groups, $group);
		} else {
		    $global_group_counts{$group}++;
		}
	    } else {
		$tag_group{$genome_tag} = '$*OTHER*$';
		if (!defined $global_group_counts{'$*OTHER*$'}) {
		    $global_group_counts{'$*OTHER*$'} = 1;
		    push(@groups, '$*OTHER*$');
		} else {
		    $global_group_counts{'$*OTHER*$'}++;
		}
	    }
	}
	close ($genome_labels_file);
	my $fatal_error = 0;
	for my $genome_tag (@genome_tags) {
	    if ($unique_genome_tags{$genome_tag} eq $genome_tag) {
		$fatal_error = 1;
		print STDERR "$genome_tag appeared in the genome tags file $genome_tags_file_name but not in the group labels file $genome_labels_file_name\n";
	    }
	}
	if ($fatal_error) {
	    die;
	}
    }
    
    return;
}

sub read_adj_mat {

    my ($adjvec_file_name, $clus_present, $adj_mat, $screen_nodes, $screen_present, $store_adj_vec) = @_;

    my $end1 = "";
    my $end2 = "";
    my $weight;
    my $adjvecfile;
    my $global_largest_edge = 0;
    my $clus;
    my $end_type;

    for my $clus (1 .. $num_clusters) {
	$clus_present->[$clus] = 0;
    }
    unless (open ($adjvecfile, "<", $adjvec_file_name) )  {
	die ("ERROR: can not open file $adjvec_file_name.\n");
    }
    while (<$adjvecfile>) {
	my @adjvec_line = ();

	chomp;
	@adjvec_line = split(/\t/, $_);  # split the scalar $adjvec_line on tab
	my $edge = shift(@adjvec_line);
	if ($edge =~ /\((\d+_[35]),(\d+_[35])\)/) {
	    ($end1, $end2) = ($edge =~ /\((\d+_[35]),(\d+_[35])\)/);
	} else {
	    die ("ERROR: input lines must be tab delimited with first field of the form:\n(#_[3 or 5],#_[3 or 5])\nNOT:\n$_\n");
	}
	if ($store_adj_vec) {
	    $adj_vec{$end1}{$end2} = \@adjvec_line;
	}
	($clus, $end_type) = split('_', $end1);
	if (($clus < 1) || ($clus > $num_clusters) || (($end_type ne "5") && ($end_type ne "3"))) {
	    die ("ERROR: cluster end ($end1) must be between 1 and $num_clusters followed by _5 or _3\n");
	}
	if ($screen_nodes && $screen_present->[$clus]) { # do not allow edges from clusters in the screen list
	    next; # skip edges that involve clusters in screen_present
	}
	$clus_present->[$clus] = 1;
	($clus, $end_type) = split('_', $end2);
	if (($clus < 1) || ($clus > $num_clusters) || (($end_type ne "5") && ($end_type ne "3"))) {
	    die ("ERROR: cluster end ($end2) must be between 1 and $num_clusters followed by _5 or _3\n");
	}
	if ($screen_nodes && $screen_present->[$clus]) { #allow edges to clusters in the screen list
	   # next;
	} else { # but do not count the clusters as present
	    $clus_present->[$clus] = 1;
	}
	$weight = 0;
	foreach my $gen_pres (@adjvec_line) {
	    $weight += $gen_pres;
	}
	if (defined($adj_mat->{$end1}{$end2})) {
	    die ("ERROR: $end1 , $end2 was previously defined in file: $adjvec_file_name\n");
	}
	$adj_mat->{$end1}{$end2} = $weight;
	if ($weight > $global_largest_edge) {
	    $global_largest_edge = $weight;
	}
    }
    close ($adjvecfile);
    for my $clus (1 .. $num_clusters) {
	if (!$clus_present->[$clus] && (!$screen_nodes || !$screen_present->[$clus])) {
	    print STDERR "cluster $clus was not present in the input file $adjvec_file_name\n" if $DEBUG;
	}
    }
    for $end1 (keys %{ $adj_mat }) {
	($clus, $end_type) = split('_', $end1);
	if ($screen_nodes && $screen_present->[$clus]) {
	    die ("ERROR: cluster end ($end1) should have been screened\n");
	}
	for $end2 (keys %{ $adj_mat->{$end1} }) {
	    ($clus, $end_type) = split('_', $end2);
	    if ($screen_nodes && $screen_present->[$clus]) {
		next; # skip edges that involve clusters in screen_present
	    }
	    if ($adj_mat->{$end1}{$end2} ne $adj_mat->{$end2}{$end1}) {
		die ("ERROR: adjacency matrix not symmetric for $end1 , $end2 in file: $adjvec_file_name\n");
	    }
	}
    }
    return($global_largest_edge);
}

sub read_centroids {

    my @line = ();
    my @centroid_id;
    my $centroid_num;
    my $centroid_name;
    my $size;
    my $id;
    my $title = "";
    my $sequence = "";
    my $centroidfile;
    my $genesfile;
    
    unless (open ($centroidfile, "<", $centroids_file_name) )  {
	die ("Can't open file $centroids_file_name.\n");
    }
    if ($output_centroids) {
	unless (open ($genesfile, ">", $genes_file_name) )  {
	    die ("Can't open file $genes_file_name.\n");
	}
    }
    my ($save_input_separator) = $/;
    $/="\n>";
    while (<$centroidfile>) {
	($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool) also removes the leading > and new line from the header
	@line = split(/\s+/, $title);  # split the scalar $title on space or tab (to separate the identifier from the header and store in array @line
	$id = $line[0]; #centroid_number expected here
	@centroid_id = split('_', $id);  # split the scalar $id at the _ character
	if ($centroid_id[0] ne "centroid") {
	    die ("ERROR: centroid fasta header line id: $id - not in expected format of >centroid_\n");
	}
	if ($centroid_id[1] !~ /[0-9]+/) {
	    die ("ERROR: centroid fasta header line id: $id - not in expected format of >centroid_(an integer)\n");
	}
	$centroid_num = $centroid_id[1];
	$centroid_name = $cluster_prefix . "_" . $centroid_num;
	if ($output_centroids) {
	    print $genesfile ">$centroid_name\n";
	    print $genesfile $sequence;
	}
	$sequence =~ s/\n//g;
	$centroid_length[$centroid_num] = length($sequence);
	$centroid_anno[$centroid_num] = join(' ', @line[2 .. $#line]);
	$title = ""; # clear the title for the next round.
	$sequence = ""; #clear out the sequence for the next round.
    }
    $/ = $save_input_separator; # restore the input separator
    close ($centroidfile);
    if ($output_centroids) {
	close ($genesfile);
    }

    return;
}

sub process_adj_mat_cross {

    my ($high_adj_mat, $low_adj_mat, $clus_present, $clus_start, $high_adj_best, $low_adj_best, $high_adj_mutual, $low_adj_mutual) = @_;

    for my $end1 (keys %{ $low_adj_mat }) { # determine best edges for each cluster end
	(my $clus, my $end) = split('_', $end1);
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_start->[$clus] != 1) {
	    next; # we are only adding best edges between low chain ends and high chain clusters for context
	}
	my @sort_end2 = (sort {
	    if ($low_adj_mat->{$end1}{$b} <=> $low_adj_mat->{$end1}{$a}) {
		return ($low_adj_mat->{$end1}{$b} <=> $low_adj_mat->{$end1}{$a});
	    } else {
		(my $clus_a, my $end_a) = split('_', $a);
		(my $clus_b, my $end_b) = split('_', $b);
		if ($cluster_size[$clus_b] <=> $cluster_size[$clus_a]) {
		    return ($cluster_size[$clus_b] <=> $cluster_size[$clus_a]);
		} else {
		    return ($clus_a <=> $clus_b);
		}
	    }
			 } keys %{ $low_adj_mat->{$end1} });
	$low_adj_best->{$end1} = $sort_end2[0];
	print STDERR "best edge $end1 $sort_end2[0] $low_adj_mat->{$end1}{$sort_end2[0]}\n" if ($DEBUG);
    }
    for my $end1 (keys %{ $high_adj_mat }) { # add edge weights from high to low
	for my $end2 (keys %{ $high_adj_mat->{$end1} }) {
	    $low_adj_mat->{$end1}{$end2} = $high_adj_mat->{$end1}{$end2};
	}
    }
    for my $end (keys %{ $high_adj_best }) { # add high best edges to low best edges
	$low_adj_best->{$end} = $high_adj_best->{$end};
	print STDERR "best edge added from high to low $end $high_adj_best->{$end}\n" if ($DEBUG);
    }
    for my $end (keys %{ $high_adj_mutual }) { # add high mutual best edges to low mutual best edges
	$low_adj_mutual->{$end} = $high_adj_mutual->{$end};
	print STDERR "mutual best edge added from high to low $end $high_adj_mutual->{$end}\n" if ($DEBUG);
    }

    return;
}

sub process_not_present { # give minimal definition to clusters without edges so we can output them

    if ($combine_high_low) {
	for my $clus (1 .. $num_clusters) {
	    if ($high_clus_present[$clus] || $low_clus_present[$clus]) { # cluster was already processed
		next;
	    }
	    $low_clus_present[$clus] = 1;
	    $low_clus_start[$clus] = 1;
	    $low_chain_length{$clus} = 1;
	    $low_chain_gene_length{$clus} = $centroid_length[$clus] * $scale_protein;
	    $max_chain_node{$clus} = $cluster_size[$clus];
	    $max_chain_edge{$clus} = 0;
	    # these clusters have no edges so set them up minimally
	    print STDERR "$clus:$low_clus_start[$clus]:$low_chain_length{$clus}:$max_chain_node{$clus}\n" if ($DEBUG);
	}
    } else {
	for my $clus (1 .. $num_clusters) {
	    if ($high_clus_present[$clus]) { # cluster was already processed
		next;
	    }
	    $high_clus_present[$clus] = 1;
	    $high_clus_start[$clus] = 1;
	    $high_chain_length{$clus} = 1;
	    $high_chain_gene_length{$clus} = $centroid_length[$clus] * $scale_protein;
	    $max_chain_node{$clus} = $cluster_size[$clus];
	    $max_chain_edge{$clus} = 0;
	    # these clusters have no edges so set them up minimally
	    print STDERR "$clus:$high_clus_start[$clus]:$high_chain_length{$clus}:$max_chain_node{$clus}\n" if ($DEBUG);
	}
    }

    return;
}

sub process_adj_mat {

    my ($global_largest_edge, $clus_present, $adj_mat, $clus_start, $adj_best, $adj_mutual, $chain_length, $chain_gene_length, $restrict_mutual, $is_core) = @_;

    my @clus_visited = ();

    for my $end1 (keys %{ $adj_mat }) { # determine best edges for each cluster end
	my @sort_end2 = (sort {
	    if ($adj_mat->{$end1}{$b} <=> $adj_mat->{$end1}{$a}) {
		return ($adj_mat->{$end1}{$b} <=> $adj_mat->{$end1}{$a});
	    } else {
		(my $clus_a, my $end_a) = split('_', $a);
		(my $clus_b, my $end_b) = split('_', $b);
		if ($cluster_size[$clus_b] <=> $cluster_size[$clus_a]) {
		    return ($cluster_size[$clus_b] <=> $cluster_size[$clus_a]);
		} else {
		    return ($clus_a <=> $clus_b);
		}
	    }
			 } keys %{ $adj_mat->{$end1} });
	$adj_best->{$end1} = $sort_end2[0];
	print STDERR "best edge $end1 $sort_end2[0] $adj_mat->{$end1}{$sort_end2[0]}\n" if ($DEBUG);
    }
    for my $end1 (keys %{ $adj_mat }) { # determine mutually best edges if they exist for each cluster end
	my $end2 = $adj_best->{$end1};
	if (!defined($end2) || !defined($adj_best->{$end2})) {
	    next;
	}
	if ($end1 eq $adj_best->{$end2}) {
	    $adj_mutual->{$end1} = $end2;
	    print STDERR "mutual best edge $end1 $end2\n" if ($DEBUG);
	}
    }
    for my $end1 (keys %{ $adj_mat }) { # check mutually best edges are symmetric
	my $end2 = $adj_mutual->{$end1};
	if (!defined($end2)) {
	    next;
	    }
	if (!defined($adj_mutual->{$end2}) || ($end1 ne $adj_mutual->{$end2})) {
	    $adj_mutual->{$end2} = $end1;
	    print STDERR "WARNING: adj_mutual hash should have been symmetric but was not for edge $end1 $end2\n";
	}
    }
    for my $clus (1 .. $num_clusters) { # define chain start clusters to have no mutual best edge on at least one end
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	$clus_visited[$clus] = 0;
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	if (!defined ($adj_mutual->{$clus_5}) || !defined ($adj_mutual->{$clus_3})) {
	    $clus_start->[$clus] = 1;
	} else {
	    $clus_start->[$clus] = 0;
	}
	print STDERR "start1 $clus:$clus_start->[$clus]\n" if ($DEBUG);
    }
    for my $clus (1 .. $num_clusters) { # attempt to repair small bubbles in the graph by adjusting mutual best edges
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_start->[$clus] == 0) { # we haven't found any cyles yet so this is just chains
	    next;
	}
	my $next_end;
	my $cur_end;
	my $prev_end;
	my $start_end;
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	my $min_edge1;
	my $min_edge2;
	my $max_edge1;
	my $max_edge2;
	my $path1_len = 0;
	my $path2_len = 0;
	print STDERR "bubble search cluster $clus\n" if ($DEBUG);
	# a bubble looks like  ----------------\/
	#                           /\------------------
	# where the ---s are mutually best edges and the \/ or /\ are the end of a chain whose
	# best edge point to the other chain; we want to collapse this into a single chain bey
	# choosing the "stronger" overlapping segment
	if (defined ($adj_mutual->{$clus_5})) {
	    $cur_end = $clus_3;
	} elsif (defined ($adj_mutual->{$clus_3})) {
	    $cur_end = $clus_5;
	} else {
	    print STDERR "skipping singleton\n" if ($DEBUG);
	    next; # do not try to repair singletons
	}
	print STDERR "end $cur_end " if ($DEBUG);
	$next_end = $adj_best->{$cur_end}; # this is the end of the chain that does not have a mutually best edge
	if (!defined $next_end) {
	    print STDERR "has no best edge\n" if ($DEBUG);
	    next; # no best edge to try to upgrade to mutually best
	} else {
	    print STDERR "best edge to $next_end\n" if ($DEBUG);
	    # initialize the strengths of the two competing paths
	    $max_edge1 = $min_edge1 = $adj_mat->{$cur_end}{$next_end};
	    $max_edge2 = 0;
	    $min_edge2 = $global_largest_edge;
	    # need to push edge to other side of cluster for following loop to work
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    if ($end eq "_5") {
		$next_end = $next_clus . "_3";
	    } else {
		$next_end = $next_clus . "_5";
	    }
	}
	$start_end = $next_end;
	while (defined($next_end)) { # follow path2 back to the end of the potential overlap
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    if ($end eq "_5") {
		$prev_end = $next_clus . "_3";
	    } else {
		$prev_end = $next_clus . "_5";
	    }
	    $next_end = $adj_mutual->{$prev_end};
	    if (defined $next_end) {
		$path2_len++;
		if ($adj_mat->{$prev_end}{$next_end} > $max_edge2) {
		    $max_edge2 = $adj_mat->{$prev_end}{$next_end};
		}
		if ($adj_mat->{$prev_end}{$next_end} < $min_edge2) {
		    $min_edge2 = $adj_mat->{$prev_end}{$next_end};
		}
		print STDERR "following mutual best edge to $next_end min2:$min_edge2 max2:$max_edge2\n" if ($DEBUG);
		if ($next_end eq $start_end) {  # we are in a cycle so quit
		    print STDERR "In a cycle so quitting\n" if ($DEBUG);
		    last;
		}
	    }
	} # $prev_end has the last end of path2
	if ((defined $prev_end) && (defined $adj_best->{$prev_end})) {
	    if (defined $adj_mutual->{$prev_end}) {
		print STDERR "In a cycle so quitting again\n" if ($DEBUG);
		next; # we are in a cycle so quit
	    }
	    if ($prev_end eq $cur_end) {
		print STDERR "Looped back on the same chain so quitting\n" if ($DEBUG);
		next; # we have looped back on the same chain so quit
	    }
	    # need to see where $prev_end best edge goes
	    my $alt_end = $prev_end; # store end of path2 in $alt_end
	    $next_end = $adj_best->{$prev_end};
	    print STDERR "alt_end = $alt_end : cur_end = $cur_end : next_end = $next_end\n" if ($DEBUG);
	    if ($adj_mat->{$prev_end}{$next_end} > $max_edge2) {
		$max_edge2 = $adj_mat->{$prev_end}{$next_end};
	    }
	    if ($adj_mat->{$prev_end}{$next_end} < $min_edge2) {
		$min_edge2 = $adj_mat->{$prev_end}{$next_end};
	    }
	    print STDERR "min2:$min_edge2 max2:$max_edge2\n" if ($DEBUG);
	    # need to push edge to other side of cluster for following loop to work
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    if ($end eq "_5") {
		$next_end = $next_clus . "_3";
	    } else {
		$next_end = $next_clus . "_5";
	    }
	    $start_end = $next_end;
	    while (defined($next_end)) { # follow path1 back to the end of the potential overlap
		(my $next_clus, my $end) = split('_', $next_end);
		$end = '_' . $end;
		if ($end eq "_5") {
		    $prev_end = $next_clus . "_3";
		} else {
		    $prev_end = $next_clus . "_5";
		}
		$next_end = $adj_mutual->{$prev_end};
		if (defined $next_end) {
		    $path1_len++;
		    if ($adj_mat->{$prev_end}{$next_end} > $max_edge1) {
			$max_edge1 = $adj_mat->{$prev_end}{$next_end};
		    }
		    if ($adj_mat->{$prev_end}{$next_end} < $min_edge1) {
			$min_edge1 = $adj_mat->{$prev_end}{$next_end};
			}
		    print STDERR "following mutual best edge to $next_end min1:$min_edge1 max1:$max_edge1\n" if ($DEBUG);
		    if ($next_end eq $start_end) {  # we are in a cycle so quit
			print STDERR "In a cycle so quitting\n" if ($DEBUG);
			last;
		    }
		}
	    } # $prev_end has the end of path1
	    if ((defined $prev_end) && ($prev_end eq $cur_end)) { # we have confirmed the bubble with path1 overlapping path2
		if (defined $adj_mutual->{$prev_end}) { # this should probably not happen
		    print STDERR "In a cycle so quitting again\n" if ($DEBUG);
		    next; # we are in a cycle so quit
		}
		# repair bubble by using strongest insertion variant
		my $tot_edge1 = $min_edge1 + $max_edge1;
		my $tot_edge2 = $min_edge2 + $max_edge2;
		print STDERR "path1($path1_len):$cur_end:$min_edge1:$max_edge1:$tot_edge1\npath2($path2_len):$alt_end:$min_edge2:$max_edge2:$tot_edge2\n" if ($DEBUG);
		if ($tot_edge1 < $tot_edge2) { # choose which path to use in the joined chain
		    $cur_end = $alt_end; # chose path 2
		}
		(my $next_clus, my $end) = split('_', $cur_end);
		my $opp_end;
		$end = '_' . $end;
		if ($end eq "_5") {
		    $opp_end = $next_clus . "_3";
		} else {
		    $opp_end = $next_clus . "_5";
		}
		if (defined ($adj_mutual->{$opp_end})) { # need to make sure we did not add a singleton before resetting start which only applies for choosing path 2
		    $clus_start->[$next_clus] = 0; # since we are joing the two chains we have to remove the chain start mark
		}
		$next_end = $adj_best->{$cur_end};
		$adj_mutual->{$cur_end} = $next_end; # join the two chains by adding a mutually best edge - add the other later (see below)
		$prev_end = $next_end;
		if (defined $adj_mutual->{$next_end}) {
		    $next_end = $adj_mutual->{$next_end}; # this is where the chain used to go
		} else {
		    $next_end = $adj_best->{$next_end}; # in case the edge wasn't mutually best if path 2 or path 1 was zero length
		}
		$adj_mutual->{$prev_end} = $cur_end; # finish the join to maintain symmetry of mutually best edges (see above)
		$adj_best->{$prev_end} = $cur_end; # need to keep best consistent with mutual
		print STDERR "repair bubble $cur_end:$prev_end $next_clus $clus_start->[$next_clus] $next_end" if ($DEBUG);
		if ((defined ($adj_mutual->{$next_end})) && ($adj_mutual->{$next_end} eq $prev_end)) { # the path not chosen was not zero length
		    delete $adj_mutual->{$next_end}; # we reassigned the other end of this so need to delete this one to keep symmetry and break this path from the chain
		    ($next_clus, $end) = split('_', $next_end);
		    $clus_start->[$next_clus] = 1; # need to mark this as the start of the broken chain
		    print STDERR " delete:$next_end $next_clus $clus_start->[$next_clus]\n" if ($DEBUG);
		} else { # the path not chosen was 0 length
		    (my $prev_clus, my $end) = split('_', $prev_end);
		    my $other_end;
		    $end = '_' . $end;
		    if ($end eq "_5") {
			$other_end = $prev_clus . "_3";
		    } else {
			$other_end = $prev_clus . "_5";
		    }
		    if (defined ($adj_mutual->{$other_end})) { # need to make sure we did not add a singleton before resetting start which only applies for choosing path 1
			$clus_start->[$prev_clus] = 0;
		    }
		    print STDERR " $prev_clus $clus_start->[$prev_clus]\n" if ($DEBUG);
		}
	    }
	}
    }
    for my $clus (1 .. $num_clusters) { # mark chains as visited so we can find the cycles
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_visited[$clus] > 0){
	    next;
	}
	if ($clus_start->[$clus] == 0) {
	    next;
	}
	$clus_visited[$clus] = 1;
	my $next_end;
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	if (defined ($adj_mutual->{$clus_5})) {
	    $next_end = $adj_mutual->{$clus_5};
	} elsif (defined ($adj_mutual->{$clus_3})) {
	    $next_end = $adj_mutual->{$clus_3};
	}
	while (defined($next_end)) {
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    $clus_visited[$next_clus] = 1;
	    if ($end eq "_5") {
		$next_end = $adj_mutual->{$next_clus . "_3"};
	    } else {
		$next_end = $adj_mutual->{$next_clus . "_5"};
	    }
	}
	print STDERR "end of start $clus\n" if ($DEBUG);
    }
    for my $clus (1 .. $num_clusters) { # find cycles break at weakest edge (or smallest cluster for core)
	# set start to 2 for cycles
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_visited[$clus] > 0){
	    next;
	}
	$clus_visited[$clus] = 1;
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	if (!defined ($adj_mutual->{$clus_5})) {
	    die("ERROR: nonvisited clusters should have defined adj_mutual edges!\n");
	}
	my $cur_end = $clus_5;
	my $prev_end = $clus_3;
	my $next_end = $adj_mutual->{$cur_end};
	my $min_clus = $clus;
	my $min_clus_num = $clus;
	my $min_weight = $adj_mat->{$cur_end}{$next_end};
	print STDERR "cycle start $prev_end $cur_end $next_end $min_clus $min_weight $min_clus_num\n" if ($DEBUG);
	while ($next_end ne $prev_end) {
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    $clus_visited[$next_clus] = 1;
	    if ($end eq "_5") {
		$cur_end = $next_clus . "_3";
		$next_end = $adj_mutual->{$cur_end};
	    } else {
		$cur_end = $next_clus . "_5";
		$next_end = $adj_mutual->{$cur_end};
	    }
	    if (!defined($next_end)) {
		die("ERROR: at this stage adj_mutual edges should form a cycle!\n");
	    }
	    if ($adj_mat->{$cur_end}{$next_end} < $min_weight) {
		$min_weight = $adj_mat->{$cur_end}{$next_end};
		$min_clus = $next_clus;
	    }
	    if ($next_clus < $min_clus_num) {
		$min_clus_num = $next_clus;
	    }
	    print STDERR "cycle $cur_end $next_end $min_clus $min_weight $min_clus_num\n" if ($DEBUG);
	}
	if ($is_core) {
	    $clus_start->[$min_clus_num] = 2;
	} else {
	    $clus_start->[$min_clus] = 2;
	}
    }
    for my $clus (1 .. $num_clusters) {
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_visited[$clus] == 0){
	    die("ERROR: all clusters should have been visited during input: $clus missing!\n");
	}
    }
    if ($restrict_mutual) { # remove mutual best edges which are less than half the weight of previous mutual edge
	print STDERR "Restricting mutual best edges\n" if ($DEBUG);
	for my $clus (1 .. $num_clusters) {
	    if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
		next;
	    }
	    $clus_visited[$clus] = 0;
	}
	for my $clus (1 .. $num_clusters) {
	    if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
		next;
	    }
	    if ($clus_start->[$clus] == 0) {
		next;
	    }
	    if ($clus_visited[$clus]) {
		next; # we've already been here
	    }
	    $clus_visited[$clus] = 1;
	    my $clus_5 = $clus . "_5";
	    my $clus_3 = $clus . "_3";
	    my $prev_weight;
	    my $next_weight;
	    my $prev_end;
	    my $next_end;
	    my $cycle_clus = 0; # cluster numbering starts at 1 so this means there is no cycle
	    if (defined ($adj_mutual->{$clus_5})) {
		if (defined ($adj_mutual->{$clus_3})) {
		    # we are at a cycle start
		    $cycle_clus = $clus;
		    my $tmp1_end = $adj_mutual->{$clus_3};
		    my $tmp2_end = $adj_mutual->{$clus_5};
		    if ($adj_mat->{$clus_5}{$tmp2_end} < $adj_mat->{$clus_3}{$tmp1_end}) { # for cycles start at weakest edge
			$prev_end = $clus_3;
		    } else {
			$prev_end = $clus_5;
		    }
		} else {
		    # we are not at a cycle start
		    $prev_end = $clus_5;
		}
	    } elsif (defined ($adj_mutual->{$clus_3})) {
		# we are not at a cycle start
		$prev_end = $clus_3;
	    } else {
		# singleton cluster
		next;
	    }
	    $next_end = $adj_mutual->{$prev_end};
	    if (defined($next_end)) {
		$prev_weight = $adj_mat->{$prev_end}{$next_end};
		print STDERR "starting search for weak edges $prev_end $next_end prev weight $prev_weight\n" if ($DEBUG);
		while (defined($next_end)) {
		    $next_weight = $adj_mat->{$prev_end}{$next_end};
		    (my $next_clus, my $end) = split('_', $next_end);
		    if (((4 * $next_weight) < $prev_weight) || ((4 * $prev_weight) < $next_weight)){ # remove mutual edges
			(my $prev_clus, my $p_end) = split('_', $prev_end);
			$clus_start->[$prev_clus] = 1;
			$clus_start->[$next_clus] = 1;
			$clus_start->[$cycle_clus] = 0; # if we had a cycle we are removing the cycle start, if not we set a meaningless value
			delete $adj_mutual->{$prev_end};
			delete $adj_mutual->{$next_end};
			print STDERR "deleting mutual best edges $prev_end $next_end prev weight $prev_weight next weight $next_weight\n" if ($DEBUG);
		    }
		    $prev_weight = $next_weight;
		    if ($clus_visited[$next_clus]) {
			last;
		    }
		    $clus_visited[$next_clus] = 1;
		    $end = '_' . $end;
		    if ($end eq "_5") {
			$prev_end = $next_clus . "_3";
		    } else {
			$prev_end = $next_clus . "_5";
		    }
		    $next_end = $adj_mutual->{$prev_end};
		    if (defined $next_end) {
			print STDERR "following mutual best edge to $next_end next weight $next_weight\n" if ($DEBUG);
		    }
		}
	    }
	}
    }
    for my $clus (1 .. $num_clusters) { #compute chain and cycle lengths and maximum cluster size
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_start->[$clus] == 0) {
	    next;
	}
	my $length = 1;
	my $gene_length = $centroid_length[$clus] * $scale_protein;
	my $next_end;
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	my $next_clus = $clus;
	my $end;
	my $max_size = $cluster_size[$clus];
	my $max_edge = 0;
	my $prev_end;
	if (!defined ($adj_mutual->{$clus_3})) {
	    $prev_end = $clus_5;
	    $next_end = $adj_mutual->{$clus_5};
	    if (defined $adj_best->{$clus_3}) {
		$max_edge = $adj_mat->{$clus_3}{$adj_best->{$clus_3}};
	    } else {
		$max_edge = 0;
	    }
	} else {
	    $prev_end = $clus_3;
	    $next_end = $adj_mutual->{$clus_3};
	    if (defined $adj_best->{$clus_5}) {
		$max_edge = $adj_mat->{$clus_5}{$adj_best->{$clus_5}};
	    } else {
		$max_edge = 0;
	    }
	}
	while (defined($next_end)) {
	    if ($adj_mat->{$prev_end}{$next_end} > $max_edge) {
		$max_edge = $adj_mat->{$prev_end}{$next_end};
	    }
	    ($next_clus, $end) = split('_', $next_end);
	    if ($clus_start->[$next_clus] == 2) {
		last; # we have completed the cycle
	    }
	    if ($cluster_size[$next_clus] > $max_size) {
		$max_size = $cluster_size[$next_clus];
	    }
	    $end = '_' . $end;
	    if ($end eq "_5") {
		$prev_end = $next_clus . "_3";
	    } else {
		$prev_end = $next_clus . "_5";
	    }
	    $next_end = $adj_mutual->{$prev_end};
	    $length++;
	    $gene_length += $centroid_length[$next_clus] * $scale_protein;
	}
	if ((defined $prev_end) && (defined $adj_best->{$prev_end})) {
	    my $last_edge = $adj_mat->{$prev_end}{$adj_best->{$prev_end}};
	    if ($last_edge > $max_edge) {
		$max_edge = $last_edge;
	    }
	}
	$chain_length->{$clus} = $length;
	$chain_length->{$next_clus} = $length;
	$chain_gene_length->{$clus} = $gene_length;
	$chain_gene_length->{$next_clus} = $gene_length;
	$max_chain_node{$clus} = $max_size;
	$max_chain_node{$next_clus} = $max_size;
	$max_chain_edge{$clus} = $max_edge;
	$max_chain_edge{$next_clus} = $max_edge;
	print STDERR "end of chain length $clus $next_clus$length\n" if ($DEBUG);
    }
    for my $end1 (keys %{ $adj_mat }) { # check mutually best edges are symmetric
	my $end2 = $adj_mutual->{$end1};
	if (!defined($end2)) {
	    next;
	    }
	if (!defined($adj_mutual->{$end2}) || ($end1 ne $adj_mutual->{$end2})) {
	    $adj_mutual->{$end2} = $end1;
	    print STDERR "WARNING: adj_mutual hash should have been symmetric but was corrupted for edge $end1 $end2\n";
	}
    }

    return;
}

sub print_node  {

    state $base_coordinate = 1;
    state $last_assembly = 0;
    my ($file, $clus, $prefix, $type, $strand, $assembly_num, $beforeORafter, $edge_weight, $prev_end, $next_end) = @_;
    my $label = $prefix . "_" . $clus;
    my $type_label = $type ? "P" : "S";
    my $coord_5;
    my $coord_3;

    if ($last_assembly != $assembly_num) {
	$last_assembly = $assembly_num;
	$base_coordinate = 1;
    }
    if (($edge_weight > 0) && ($beforeORafter eq "before")) {
	my $genomes = "";
	if ($comma_tags eq "2") {
	    $genomes = "\t" . &genome_tags_string($adj_vec{$prev_end}{$next_end});
	} elsif ($comma_tags eq "3") {
	    $genomes = "\t" . &genome_groups_string($adj_vec{$prev_end}{$next_end});
	}
	print STDOUT "edge:$edge_weight$genomes\n";
    }
    print STDOUT "$label\t$strand\t$type_label\t$cluster_size[$clus]\t$centroid_anno[$clus]\n";
    if (($edge_weight > 0) && ($beforeORafter eq "after")) {
	my $genomes = "";
	if ($comma_tags eq "2") {
	    $genomes = "\t" . &genome_tags_string($adj_vec{$prev_end}{$next_end});
	} elsif ($comma_tags eq "3") {
	    $genomes = "\t" . &genome_groups_string($adj_vec{$prev_end}{$next_end});
	}
	print STDOUT "edge:$edge_weight$genomes\n";
    }
    if (!$type) {
	$label = "CONTEXT" . $context_count[$clus] . ":" . $label;
	$context_count[$clus]++;
    }
    $coord_5 = $base_coordinate;
    $base_coordinate += $centroid_length[$clus] * $scale_protein;
    $coord_3 = $base_coordinate - 1;
    if ($strand eq "-") {
	my $swap = $coord_5;
	$coord_5 = $coord_3;
	$coord_3 = $swap;
    }
    if ($type) {
	$cluster_strand[$clus] = $strand;
	$end_coord{$clus . "_5"} = $coord_5;
	$end_coord{$clus . "_3"} = $coord_3;
	print STDERR "$assembly_num $clus $cluster_strand[$clus] $end_coord{$clus . '_5'} $end_coord{$clus . '_3'}\n" if $DEBUG;
    }
    print $file "$assembly_num\t$label\t$coord_5\t$coord_3\t$centroid_anno[$clus]\t$prefix\t$centroid_length[$clus]\n";

    return;
}

sub print_up_context  {

    my ($file, $stop_end, $prev_end, $next_end, $number, $prefix, $assembly_num, $beforeORafter, $edge_weight, $adj_mat, $adj_mutual) = @_;
    my @tmp_stack = ();

    while ((defined($next_end)) && $number-- && ($next_end ne $stop_end)) {
	(my $next_clus, my $end) = split('_', $next_end);
	$end = '_' . $end;
	if ($end eq "_5") {
	    push (@tmp_stack, [$next_clus, $prefix, 0, "-", $edge_weight, $prev_end, $next_end]);
	    $prev_end = $next_clus . "_3";
	} else {
	    push (@tmp_stack, [$next_clus, $prefix, 0, "+", $edge_weight, $prev_end, $next_end]);
	    $prev_end = $next_clus . "_5";
	}
	$next_end = $adj_mutual->{$prev_end};
	if (defined($next_end)) {
	    $edge_weight = $adj_mat->{$prev_end}{$next_end};
	}
    }
    while (defined (my $next_params = pop (@tmp_stack))) {
	&print_node ($file, $next_params->[0], $next_params->[1], $next_params->[2], $next_params->[3], $assembly_num, $beforeORafter, $next_params->[4], $next_params->[5], $next_params->[6]);
    }

    return;
}

sub print_down_context  {

    my ($file, $stop_end, $prev_end, $next_end, $number, $prefix, $assembly_num, $beforeORafter, $edge_weight, $adj_mat, $adj_mutual) = @_;

    while ((defined($next_end)) && $number-- && ($next_end ne $stop_end)) {
	(my $next_clus, my $end) = split('_', $next_end);
	$end = '_' . $end;
	if ($end eq "_5") {
	    &print_node($file, $next_clus, $prefix, 0, "+", $assembly_num, $beforeORafter, $edge_weight, $prev_end, $next_end);
	    $prev_end = $next_clus . "_3";
	} else {
	    &print_node($file, $next_clus, $prefix, 0, "-", $assembly_num, $beforeORafter, $edge_weight, $prev_end, $next_end);
	    $prev_end = $next_clus . "_5";
	}
	$next_end = $adj_mutual->{$prev_end};
	if (defined($next_end)) {
	    $edge_weight = $adj_mat->{$prev_end}{$next_end};
	}
    }
    return;
}

sub print_consensus {

    my ($attribute_file_name, $clus_present, $adj_mat, $clus_start, $adj_best, $adj_mutual, $chain_length, $chain_gene_length, $suffix, $assembly_num, $is_core) = @_;

    my $attributefile;
    my @clus_visited = ();
    my @sorted_by_length = (sort {
	if ($clus_present->[$b] <=> $clus_present->[$a]) {
	    return ($clus_present->[$b] <=> $clus_present->[$a]);
	} elsif (!$clus_present->[$b]) {
	    return (0);
	} elsif ($clus_start->[$b] <=> $clus_start->[$a]) {
	    return ($clus_start->[$b] <=> $clus_start->[$a]);
	} elsif (!$clus_start->[$b]) {
	    return (0);
	} elsif ($chain_length->{$b} <=> $chain_length->{$a}) {
	    return ($chain_length->{$b} <=> $chain_length->{$a});
	} elsif ($max_chain_edge{$b} <=> $max_chain_edge{$a}) {
	    return ($max_chain_edge{$b} <=> $max_chain_edge{$a});
	} else {
	    return ($max_chain_node{$b} <=> $max_chain_node{$a});
	}
			    } (1 .. $num_clusters));

    unless (open ($attributefile, ">", $attribute_file_name) )  {
	die ("Can't open file $attribute_file_name.\n");
    }
    for my $clus (1 .. $num_clusters) {
	$context_count[$clus] = 1;
	$clus_visited[$clus] = 0;
    }
    for my $clus (@sorted_by_length) {
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_visited[$clus] > 0){
	    next;
	}
	if ($clus_start->[$clus] == 0) {
	    next;
	}
	my $chain_type = ($clus_start->[$clus] == 1) ? "chain" : "cycle";
	my $prot_type = ($scale_protein == 3) ? "nt" : "aa";
	print STDOUT "#Assembly$suffix\t$assembly_num\t$chain_type\t$chain_length->{$clus}\t$chain_gene_length->{$clus} $prot_type\n";
	print STDERR "$clus:$clus_start->[$clus]\n" if ($DEBUG);
	$clus_visited[$clus] = 1;
	$assembly_num_clus[$clus] = $assembly_num;
	my $local_context_length;
	if (($is_core) && ($clus_start->[$clus] == 2)) { # do not output initial context for core cycles
	    $local_context_length = 0;
	} else {
	    $local_context_length = $context_length;
	}
	my $next_end;
	my $prev_end;
	my $cur_end; # needed for cycles in the graph
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	if (defined ($adj_mutual->{$clus_5})) {
	    $prev_end = $clus_5;
	    $cur_end = $clus_3;
	    $next_end = $adj_mutual->{$clus_5};
	    if (defined ($adj_mutual->{$clus_3})) {
		my $tmp_end = $adj_mutual->{$clus_3};
		(my $next_clus_5, my $end1) = split('_', $next_end);
		(my $next_clus_3, my $end2) = split('_', $tmp_end);
		if (($is_core && ($next_clus_3 < $next_clus_5)) || ((!$is_core) && ($adj_mat->{$clus_5}{$next_end} < $adj_mat->{$clus_3}{$tmp_end}))) { #for cycles start at weakest edge or smallest cluster for core
		    print STDERR "branch1\n" if ($DEBUG);
		    $next_end = $tmp_end;
		    $prev_end = $clus_3;
		    $cur_end = $clus_5;
		    if (defined ($adj_best->{$clus_5})) {
			&print_up_context($attributefile, $clus_3, $clus_5, $adj_best->{$clus_5}, $local_context_length, $cluster_prefix, $assembly_num, "after", $adj_mat->{$clus_5}{$adj_best->{$clus_5}}, $adj_mat, $adj_mutual);
		    }
		    &print_node($attributefile, $clus, $cluster_prefix, 1, "+", $assembly_num, "", 0, "", "");
		} else {
		    print STDERR "branch2\n" if ($DEBUG);
		    if (defined ($adj_best->{$clus_3})) {
			&print_up_context($attributefile, $clus_5, $clus_3, $adj_best->{$clus_3}, $local_context_length, $cluster_prefix, $assembly_num, "after", $adj_mat->{$clus_3}{$adj_best->{$clus_3}}, $adj_mat, $adj_mutual);
		    }
		    &print_node($attributefile, $clus, $cluster_prefix, 1, "-", $assembly_num, "", 0, "", "");
		}
	    } else {
		    print STDERR "branch3\n" if ($DEBUG);
		    if (defined ($adj_best->{$clus_3})) {
			&print_up_context($attributefile, $clus_5, $clus_3, $adj_best->{$clus_3}, $local_context_length, $cluster_prefix, $assembly_num, "after", $adj_mat->{$clus_3}{$adj_best->{$clus_3}}, $adj_mat, $adj_mutual);
		    }
		    &print_node($attributefile, $clus, $cluster_prefix, 1, "-", $assembly_num, "", 0, "", "");
	    }
	} elsif (defined ($adj_mutual->{$clus_3})) {
	    print STDERR "branch4\n" if ($DEBUG);
	    $next_end = $adj_mutual->{$clus_3};
	    $prev_end = $clus_3;
	    $cur_end = $clus_5;
	    if (defined ($adj_best->{$clus_5})) {
		&print_up_context($attributefile, $clus_3, $clus_5, $adj_best->{$clus_5}, $local_context_length, $cluster_prefix, $assembly_num, "after", $adj_mat->{$clus_5}{$adj_best->{$clus_5}}, $adj_mat, $adj_mutual);
	    }
	    &print_node($attributefile, $clus, $cluster_prefix, 1, "+", $assembly_num, "", 0, "", "");
	} else { #singleton cluster
	    print STDERR "branch5\n" if ($DEBUG);
	    if (defined ($adj_best->{$clus_3})) {
		&print_up_context($attributefile, $clus_5, $clus_3, $adj_best->{$clus_3}, $local_context_length, $cluster_prefix, $assembly_num, "after", $adj_mat->{$clus_3}{$adj_best->{$clus_3}}, $adj_mat, $adj_mutual);
	    }
	    &print_node($attributefile, $clus, $cluster_prefix, 1, "+", $assembly_num, "", 0, "", "");
	    if (defined ($adj_best->{$clus_5})) {
		&print_down_context($attributefile, $clus_3, $clus_5, $adj_best->{$clus_5}, $local_context_length, $cluster_prefix, $assembly_num, "before", $adj_mat->{$clus_5}{$adj_best->{$clus_5}}, $adj_mat, $adj_mutual);
	    }
	    $assembly_num++;
	    print STDOUT "\n\n";
	    next;
	}
	my $cycle_end;
	while ((defined($next_end)) && ($next_end ne $cur_end)) {
	    $cycle_end = $next_end;
	    print STDERR "next $next_end\n" if ($DEBUG);
	    (my $next_clus, my $end) = split('_', $next_end);
	    $end = '_' . $end;
	    $clus_visited[$next_clus] = 1;
	    $assembly_num_clus[$next_clus] = $assembly_num;
	    if ($end eq "_5") {
		&print_node($attributefile, $next_clus, $cluster_prefix, 1, "+", $assembly_num, "before", $adj_mat->{$prev_end}{$next_end}, $prev_end, $next_end);
		$prev_end = $next_clus . "_3";
	    } else {
		&print_node($attributefile, $next_clus, $cluster_prefix, 1, "-", $assembly_num, "before", $adj_mat->{$prev_end}{$next_end}, $prev_end, $next_end);
		$prev_end = $next_clus . "_5";
	    }
	    $next_end = $adj_mutual->{$prev_end};
	}
	print STDERR "prev $prev_end\n" if ($DEBUG);
	if (($is_core) && ($clus_start->[$clus] == 2)) { # only output one context for core cycles
	    $local_context_length = 1;
	} else {
	    $local_context_length = $context_length;
	}
	if (defined ($adj_best->{$prev_end})) {
	    &print_down_context($attributefile, $cycle_end, $prev_end, $adj_best->{$prev_end}, $local_context_length, $cluster_prefix, $assembly_num, "before", $adj_mat->{$prev_end}{$adj_best->{$prev_end}}, $adj_mat, $adj_mutual);
	}
	$assembly_num++;
	print STDOUT "\n\n";
    }
    for my $clus (1 .. $num_clusters) {
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_visited[$clus] == 0){
	    die("ERROR: all clusters should have been visited during output: $clus missing!\n");
	}
    }
    close ($attributefile);

    return ($assembly_num);
}

sub follow_best {

    my ($next_end, $stop_end, $adj_mat, $adj_best, $clus_present, $length, $gene_length, $prevalence, $prev_assembly) = @_;
    my $next_clus;
    my $end;
    my @tmp_array = ();
    my %tmp_hash_assembly = ();
    my %tmp_hash_end = ();
    my $next_assembly;
    my $stop_end2 = $stop_end;
    my $total_length = 0;
    my $total_gene_length = 0;
    my $max_prevalence = 0;

    print STDERR "FB:$next_end:$stop_end:$length:$gene_length:$prevalence:$prev_assembly\n" if ($DEBUG);
    $tmp_hash_assembly{$prev_assembly} = 1;

    while (defined($next_end)) {
	my $prev_end;
	($next_clus, $end) = split('_', $next_end);
	$next_assembly = $assembly_num_clus[$next_clus];
	print STDERR "$next_end:$next_clus:$next_assembly:$clus_present->[$next_clus]\n" if ($DEBUG);
	if (($next_end eq $stop_end) || ($next_end eq $stop_end2)) {
	    print STDERR "in a cycle $next_end:$stop_end:$stop_end2\n" if ($DEBUG);
	    $next_end = undef; # somehow we have gotten into a cycle and need to stop
	    last;
	}
	if (defined $tmp_hash_end{$next_end}) {
	    print STDERR "in a cycle $next_end visited already\n" if ($DEBUG);
	    $next_end = undef; # somehow we have gotten into a cycle and need to stop
	    last;
	} else {
	    $tmp_hash_end{$next_end} = 1;
	}
	if ($prev_assembly ne $next_assembly) {
	    if (defined $tmp_hash_assembly{$next_assembly}) {
		$next_end = undef; # somehow we have gotten into a cycle and need to stop
		print STDERR "in a cycle $next_assembly\n" if ($DEBUG);
		last;
	    } else {
		$tmp_hash_assembly{$next_assembly} = 1;
	    }
	    my @tmp2_array = ($prevalence, $length, $gene_length, $prev_assembly);
	    push(@tmp_array, \@tmp2_array);
	    $total_length += $length;
	    $total_gene_length += $gene_length;
	    if ($max_prevalence < $prevalence) {
		$max_prevalence = $prevalence;
	    }
	    $stop_end2 = $next_end;
	    print STDERR "$stop_end2:$length:$gene_length:$prevalence:$total_length:$total_gene_length:$max_prevalence\n" if ($DEBUG);
	    $length = 0;
	    $gene_length = 0;
	    $prevalence = 0;
	    $prev_assembly = $next_assembly;
	}
	if (!$clus_present->[$next_clus]) {
	    last; # we have found a core gene
	}
	$length++;
	$gene_length += $centroid_length[$next_clus] * $scale_protein;
	$end = '_' . $end;
	if ($end eq "_5") {
	    $prev_end = $next_clus . "_3";
	} else {
	    $prev_end = $next_clus . "_5";
	}
	$next_end = $adj_best->{$prev_end};
	if (defined $next_end) {
	    if ($prevalence < $adj_mat->{$prev_end}{$next_end}) {
		$prevalence = $adj_mat->{$prev_end}{$next_end};
	    }
	}
    }

    return ($next_end, $total_length, $total_gene_length, $max_prevalence, \@tmp_array);
}

sub print_insertions {

    my ($attribute_file_name, $clus_present, $clus_start, $adj_mat, $adj_best, $adj_mutual, $chain_length, $chain_gene_length) = @_;
    my $attributefile;
    my %insertions;
    my %islands;
    my $size_threshold = $num_genomes / 10;

    for my $clus (1 .. $num_clusters) { # store insertion points
	if (!$clus_present->[$clus]) { # cluster was not in the inputted edge file so skip
	    next;
	}
	if ($clus_start->[$clus] != 1) {
	    next; # only interested in chains for fGI's and where they insert
	}
	if (($chain_length->{$clus} < 3) && ($max_chain_edge{$clus} < $size_threshold)) {
	    print STDERR "fGI $clus:$assembly_num_clus[$clus] too short $chain_length->{$clus} and too rare $max_chain_edge{$clus}\n" if ($DEBUG);
	    next; # not counting short insertion events which are rare
	}
	if ($max_chain_edge{$clus} < 3) {
	    print STDERR "fGI $clus:$assembly_num_clus[$clus] too rare $max_chain_edge{$clus}\n" if ($DEBUG);
	    next; # not counting very rare insertion events
	}
	print STDERR "fGI $clus:$assembly_num_clus[$clus]\n" if ($DEBUG);
	my $clus_5 = $clus . "_5";
	my $clus_3 = $clus . "_3";
	my $next_end;
	my $assemblies;
	my $total_length = 0;
	my $total_gene_length = 0;
	my $max_prevalence = 0;
	if (defined ($adj_mutual->{$clus_5})) {
	    ($next_end, $total_length, $total_gene_length, $max_prevalence, $assemblies) = &follow_best($adj_best->{$clus_3}, $clus_5, $adj_mat, $adj_best, $clus_present, $chain_length->{$clus}, $chain_gene_length->{$clus}, $max_chain_edge{$clus}, $assembly_num_clus[$clus]);
	    print STDERR "insert1:$total_length:$total_gene_length:$max_prevalence\n" if ($DEBUG);
	    if (defined($next_end)) {
		$insertions{$next_end}{$clus_3}{'size'} = $max_prevalence;
		$insertions{$next_end}{$clus_3}{'length'} = $total_length;
		$insertions{$next_end}{$clus_3}{'gene'} = $total_gene_length;
		$insertions{$next_end}{$clus_3}{'ass'} = $assemblies;
	    }
	} elsif (defined ($adj_mutual->{$clus_3})) {
	    ($next_end, $total_length, $total_gene_length, $max_prevalence, $assemblies) = &follow_best($adj_best->{$clus_5}, $clus_3, $adj_mat, $adj_best, $clus_present, $chain_length->{$clus}, $chain_gene_length->{$clus}, $max_chain_edge{$clus}, $assembly_num_clus[$clus]);
	    print STDERR "insert2:$total_length:$total_gene_length:$max_prevalence\n" if ($DEBUG);
	    if (defined($next_end)) {
		$insertions{$next_end}{$clus_5}{'size'} = $max_prevalence;
		$insertions{$next_end}{$clus_5}{'length'} = $total_length;
		$insertions{$next_end}{$clus_5}{'gene'} = $total_gene_length;
		$insertions{$next_end}{$clus_5}{'ass'} = $assemblies;
	    }
	} else { # this is a singleton chain/insertion which we are currently not counting but put this here in case we decide to later
	    ($next_end, $total_length, $total_gene_length, $max_prevalence, $assemblies) = &follow_best($adj_best->{$clus_3}, $clus_5, $adj_mat, $adj_best, $clus_present, $chain_length->{$clus}, $chain_gene_length->{$clus}, $max_chain_edge{$clus}, $assembly_num_clus[$clus]);
	    print STDERR "insert3:$total_length:$total_gene_length:$max_prevalence\n" if ($DEBUG);
	    if (defined($next_end)) {
		$insertions{$next_end}{$clus_3}{'size'} = $max_prevalence;
		$insertions{$next_end}{$clus_3}{'length'} = $total_length;
		$insertions{$next_end}{$clus_3}{'gene'} = $total_gene_length;
		$insertions{$next_end}{$clus_3}{'ass'} = $assemblies;
	    }
	    ($next_end, $total_length, $total_gene_length, $max_prevalence, $assemblies) = &follow_best($adj_best->{$clus_5}, $clus_3, $adj_mat, $adj_best, $clus_present, $chain_length->{$clus}, $chain_gene_length->{$clus}, $max_chain_edge{$clus}, $assembly_num_clus[$clus]);
	    print STDERR "insert4:$total_length:$total_gene_length:$max_prevalence\n" if ($DEBUG);
	    if (defined($next_end)) {
		$insertions{$next_end}{$clus_5}{'size'} = $max_prevalence;
		$insertions{$next_end}{$clus_5}{'length'} = $total_length;
		$insertions{$next_end}{$clus_5}{'gene'} = $total_gene_length;
		$insertions{$next_end}{$clus_5}{'ass'} = $assemblies;
	    }
	}
    }

    if ($DEBUG) {
	print STDERR "Done storing insertion points\n";
	for my $clus (1 .. $num_clusters) {
	    print STDERR "$assembly_num_clus[$clus] $clus $cluster_strand[$clus] $end_coord{$clus . '_5'} $end_coord{$clus . '_3'}\n";
	}
    }
    unless (open ($attributefile, "<", $attribute_file_name) )  {
	die ("Can't open file $attribute_file_name.\n");
    }
    my $fGI_file_name = $attribute_file_name . "fGI";
    my $fGIfile;
    unless (open ($fGIfile, ">", $fGI_file_name) )  {
	die ("Can't open file $fGI_file_name.\n");
    }
    my $insertfile;
    unless (open ($insertfile, ">", $insert_file_name) )  {
	die ("Can't open file $insert_file_name.\n");
    }
    my $fGI_detail_file_name = $insert_file_name . ".details";
    my $fGI_detail_file;
    unless (open ($fGI_detail_file, ">", $fGI_detail_file_name) )  {
	die ("Can't open file $fGI_detail_file_name.\n");
    }
    my $prev_gene_len = 0;
    my $prev_core_end;
    my $prev_size = 0;
    my $prev_count = 0;
    my $insert_num = 0;
    my $prev_assembly_num = 0;
    for my $core_end (sort {
	(my $clus_a, my $end_a) = split('_', $a);
	(my $clus_b, my $end_b) = split('_', $b);
	if ($assembly_num_clus[$clus_a] <=> $assembly_num_clus[$clus_b]) {
	    return ($assembly_num_clus[$clus_a] <=> $assembly_num_clus[$clus_b]);
	} else {
	    return ($end_coord{$a} <=> $end_coord{$b});
	}
			    } (keys %insertions)) {
	my $count = 0;
	my $max_len = 0;
	my $max_gene_len = 0;
	my $max_size = 0;
	for my $fGI_end (keys %{ $insertions{$core_end} }) {
	    if ($insertions{$core_end}{$fGI_end}{'length'} > $max_len) {
		$max_len = $insertions{$core_end}{$fGI_end}{'length'};
	    }
	    if ($insertions{$core_end}{$fGI_end}{'gene'} > $max_gene_len) {
		$max_gene_len = $insertions{$core_end}{$fGI_end}{'gene'};
	    }
	    if ($insertions{$core_end}{$fGI_end}{'size'} > $max_size) {
		$max_size = $insertions{$core_end}{$fGI_end}{'size'};
	    }
	    $count++;
	}
	(my $clus, my $end) = split('_', $core_end);
	if ($assembly_num_clus[$clus] ne $prev_assembly_num) {
	    undef $prev_core_end;
	}
	if ((defined $prev_core_end) && (defined $high_adj_mutual{$prev_core_end}) && ($core_end eq $high_adj_mutual{$prev_core_end})) {
	    if ($max_gene_len > $prev_gene_len) {
		$islands{$prev_core_end}{'gene'} = $max_gene_len;
	    }
	    if ($max_size > $prev_size) {
		$islands{$prev_core_end}{'size'} = $max_size;
	    }
	    if ($count > $prev_count) {
		$islands{$prev_core_end}{'count'} = $count;
	    }
	    $islands{$prev_core_end}{'end'} = $core_end;
	} else {
	    $islands{$core_end}{'gene'} = $max_gene_len;
	    $islands{$core_end}{'size'} = $max_size;
	    $islands{$core_end}{'count'} = $count;
	    $islands{$core_end}{'num'} = ++$insert_num;
	    $islands{$core_end}{'end'} = $high_adj_mutual{$core_end};
	}
	$prev_gene_len = $max_gene_len;
	$prev_core_end = $core_end;
	$prev_assembly_num = $assembly_num_clus[$clus];
	$prev_size = $max_size;
	print $insertfile "#Core\t$core_end\t$assembly_num_clus[$clus]\t$cluster_strand[$clus]\t$end_coord{$core_end}\tmaxlen:$max_len,$max_gene_len\thas $count fGI insertions\n";
	my $insert_label = $cluster_prefix . "_" . "INS_" . $insert_num;
	print $insertfile "#$insert_label\n";
	for my $fGI_end (sort {
	    if ($insertions{$core_end}{$b}{'size'} <=> $insertions{$core_end}{$a}{'size'}) {
		return ($insertions{$core_end}{$b}{'size'} <=> $insertions{$core_end}{$a}{'size'});
	    } elsif ($insertions{$core_end}{$b}{'length'} <=> $insertions{$core_end}{$a}{'length'}) {
		return ($insertions{$core_end}{$b}{'length'} <=> $insertions{$core_end}{$a}{'length'});
	    } else {
		return ($insertions{$core_end}{$b}{'gene'} <=> $insertions{$core_end}{$a}{'gene'});
	    }
			 } (keys %{ $insertions{$core_end} })) {
	    print $insertfile "\tprev:$insertions{$core_end}{$fGI_end}{'size'}\tlen:$insertions{$core_end}{$fGI_end}{'length'}\tgene_len:$insertions{$core_end}{$fGI_end}{'gene'}";
	    foreach my $ass (@ { $insertions{$core_end}{$fGI_end}{'ass'} }) {
		my $print_string = join(',', (@ { $ass }));
		print $insertfile "\t$print_string";
	    }
	    print $insertfile "\n";
	}
    }
    undef $prev_assembly_num;
    my $offset = 0;
    my $max_start_offset = 0;
    my $assembly_num;
    my $label;
    my $coord_5;
    my $coord_3;
    my $anno;
    my $glabel;
    my $prot_len;
    my $cur_clus;
    my $coord_beg;
    my $coord_end;
    my $before_gene_end;
    my $after_gene_end;
    my $tmp_prefix = $cluster_prefix . "_"; # add underscore to prefix
    my @cur_vec;
    
    for my $index (0 .. ($num_genomes - 1)) {
	$cur_vec[$index] = 1;
    }
    for my $core_end (sort {
	(my $clus_a, my $end_a) = split('_', $a);
	(my $clus_b, my $end_b) = split('_', $b);
	if ($assembly_num_clus[$clus_a] <=> $assembly_num_clus[$clus_b]) {
	    return ($assembly_num_clus[$clus_a] <=> $assembly_num_clus[$clus_b]);
	} else {
	    return ($end_coord{$a} <=> $end_coord{$b});
	}
			    } (keys %islands)) {
	my $stop_clus;
	my $start_clus;
	my $other_end = $islands{$core_end}{'end'};
	my $insert_label = $tmp_prefix . "INS_" . $islands{$core_end}{'num'};
	print $fGI_detail_file "$insert_label\t$core_end";
	if (!defined $other_end) {
	    print $fGI_detail_file "\n";
	    $stop_clus = 0; # set to impossible value
	} else {
	    print $fGI_detail_file "\t$other_end\n";
	    (my $clus, my $end) = split('_', $other_end);
	    $stop_clus = $clus;
	}
	(my $clus, my $end) = split('_', $core_end);
	$start_clus = $clus;
	my $output_counts;
	if (($output_counts = &dfs_fgi_ins ($fGI_detail_file, $start_clus, $stop_clus, $core_end, \@cur_vec, $num_genomes, "", 1, 0, 0, 0)) != $num_genomes) {
	    die "ERROR: all genomes should be accounted for in depth first search of the inserts at $core_end instead of $output_counts of $num_genomes\n";
	}
	if (defined $other_end) {
	    if (($output_counts = &dfs_fgi_ins ($fGI_detail_file, $stop_clus, $start_clus, $other_end, \@cur_vec, $num_genomes, "", 0, 0, 0, 0)) != $num_genomes) {
		die "ERROR: all genomes should be accounted for in depth first search of the inserts at $other_end instead of $output_counts of $num_genomes\n";
	    }
	}
	my $fGI_type = "";
	if (((100 * $islands{$core_end}{'size'}) / $num_genomes) >= 75) {
	    $fGI_type = "fGI_75";
	} elsif (((100 * $islands{$core_end}{'size'}) / $num_genomes) >= 50) {
	    $fGI_type = "fGI_50";
	} elsif (((100 * $islands{$core_end}{'size'}) / $num_genomes) >= 25) {
	    $fGI_type = "fGI_25";
	} elsif (((100 * $islands{$core_end}{'size'}) / $num_genomes) >= 10) {
	    $fGI_type = "fGI_10";
	} else {
	    $fGI_type = "fGI_1";
	}
	if ((defined $after_gene_end) && ($core_end eq $after_gene_end)) {
	    my $ins_coord_5 = $coord_end + 1;
	    my $ins_coord_3 = $coord_end + $islands{$core_end}{'gene'};
	    my $ins_len = int ($islands{$core_end}{'gene'} / 3);
	    $offset += $islands{$core_end}{'gene'};
	    print $fGIfile "$assembly_num\t$insert_label\t$ins_coord_5\t$ins_coord_3\t$fGI_type;$islands{$core_end}{'count'}\tfGI\t$ins_len\n";
	    next; #go to next core_end
	}
	while (<$attributefile>) {
	    chomp;
	    ($assembly_num, $label, $coord_5, $coord_3, $anno, $glabel, $prot_len) = split(/\t/, $_);  # split the scalar $attribute_line on tab
	    if ((defined $prev_assembly_num) && ($prev_assembly_num ne $assembly_num)) {#use ne in case we don't always use numbers
		$offset = 0;
		$max_start_offset = 0;
	    }
	    $prev_assembly_num = $assembly_num;
	    if ($label =~ /CONTEXT(\d+):$tmp_prefix/) {
		if ($coord_5 > $coord_3) {
		    $max_start_offset = $coord_5;
		} else {
		    $max_start_offset = $coord_3;
		}
		next;#skip CONTEXT lines
	    }
	    if ($label =~ /$tmp_prefix(\d+)/) {
		($cur_clus) = ($label =~ /$tmp_prefix(\d+)/);
		$offset -= $max_start_offset;
		$max_start_offset = 0;
	    } else {
		die ("ERROR: $label is not in expected format in $attribute_file_name should start with $tmp_prefix\n");
	    }
	    $coord_5 += $offset;
	    $coord_3 += $offset;
	    if ($coord_5 > $coord_3) {
		$before_gene_end = $cur_clus . "_3";
		$after_gene_end = $cur_clus . "_5";
		$coord_beg = $coord_3;
		$coord_end = $coord_5;
	    } else {
		$before_gene_end = $cur_clus . "_5";
		$after_gene_end = $cur_clus . "_3";
		$coord_beg = $coord_5;
		$coord_end = $coord_3;
	    }
	    if ($core_end eq $before_gene_end) {
		my $ins_coord_5 = $coord_beg;
		my $ins_coord_3 = $coord_beg + $islands{$core_end}{'gene'} - 1;
		my $ins_len = int ($islands{$core_end}{'gene'} / 3);
		$offset += $islands{$core_end}{'gene'};
		print $fGIfile "$assembly_num\t$insert_label\t$ins_coord_5\t$ins_coord_3\t$fGI_type;$islands{$core_end}{'count'}\tfGI\t$ins_len\n";
		$coord_5 += $islands{$core_end}{'gene'};
		$coord_3 += $islands{$core_end}{'gene'};
		$coord_beg += $islands{$core_end}{'gene'};
		$coord_end += $islands{$core_end}{'gene'};
		print $fGIfile "$assembly_num\t$label\t$coord_5\t$coord_3\t$anno\t$glabel\t$prot_len\n";
		last; #go to next core_end
	    }
	    if ($core_end eq $after_gene_end) {
		my $ins_coord_5 = $coord_end + 1;
		my $ins_coord_3 = $coord_end + $islands{$core_end}{'gene'};
		my $ins_len = int ($islands{$core_end}{'gene'} / 3);
		print $fGIfile "$assembly_num\t$label\t$coord_5\t$coord_3\t$anno\t$glabel\t$prot_len\n";
		$offset += $islands{$core_end}{'gene'};
		print $fGIfile "$assembly_num\t$insert_label\t$ins_coord_5\t$ins_coord_3\t$fGI_type;$islands{$core_end}{'count'}\tfGI\t$ins_len\n";
		last; #go to next core_end
	    }
	    print $fGIfile "$assembly_num\t$label\t$coord_5\t$coord_3\t$anno\t$glabel\t$prot_len\n";
	}
    }
    while (<$attributefile>) {
	chomp;
	($assembly_num, $label, $coord_5, $coord_3, $anno, $glabel, $prot_len) = split(/\t/, $_);  # split the scalar $attribute_line on tab
	if ((defined $prev_assembly_num) && ($prev_assembly_num ne $assembly_num)) {#use ne in case we don't always use numbers
	    $offset = 0;
	    $max_start_offset = 0;
	}
	$prev_assembly_num = $assembly_num;
	if ($label =~ /CONTEXT(\d+):$tmp_prefix/) {
	    if ($coord_5 > $coord_3) {
		$max_start_offset = $coord_5;
	    } else {
		$max_start_offset = $coord_3;
	    }
	    next;#skip CONTEXT lines
	}
	if ($label =~ /$tmp_prefix(\d+)/) {
	    ($cur_clus) = ($label =~ /$tmp_prefix(\d+)/);
	} else {
	    die ("ERROR: $label is not in expected format in $attribute_file_name should start with $tmp_prefix\n");
	}
	$coord_5 += $offset;
	$coord_3 += $offset;
	print $fGIfile "$assembly_num\t$label\t$coord_5\t$coord_3\t$anno\t$glabel\t$prot_len\n";
    }
    close($insertfile);
    close($attributefile);
    close($fGIfile);
    close($fGI_detail_file);

    if (!$specified) {
	return;
    }
    my $inputinsertfile;
    unless (open ($inputinsertfile, "<", $input_insert_file_name) )  {
	die ("Can't open file $input_insert_file_name.\n");
    }
    my $fGI_specified_file_name = $insert_file_name . ".specified";
    my $fGI_specified_file;
    unless (open ($fGI_specified_file, ">", $fGI_specified_file_name) )  {
	die ("Can't open file $fGI_specified_file_name.\n");
    }
    for my $index (0 .. ($num_genomes - 1)) {
	$cur_vec[$index] = 1;
    }
    while (<$inputinsertfile>) {
	my $stop_clus;
	my $start_clus;
	chomp;
        (my $start_end, my $stop_end, my $depth_limit) = split(/\t/, $_);  # split on tab
	if ($start_end =~ /(\d+_[35])/) {
	    (my $clus, my $end) = split('_', $start_end);
	    $start_clus = $clus;
	} else {
	    die ("ERROR: input lines must be tab delimited of the form:\n#_[3 or 5] tab #_[3 or 5] tab #\nNOT:\n$_\n");
	}
	if ($stop_end =~ /(\d+_[35])/) {
	    (my $clus, my $end) = split('_', $stop_end);
	    $stop_clus = $clus;
	} else {
	    die ("ERROR: input lines must be tab delimited of the form:\n#_[3 or 5] tab #_[3 or 5] tab #\nNOT:\n$_\n");
	}
	if ((defined $depth_limit) && (($depth_limit >= 1) && ($depth_limit <= 100))) {
	    $depth_limit = int($depth_limit);
	} else {
	    $depth_limit = 100;
	}
	print $fGI_specified_file "Inserts between $start_end and $stop_end max length $depth_limit\n";
	my $output_counts;
	if (($output_counts = &dfs_fgi_ins ($fGI_specified_file, $start_clus, $stop_clus, $start_end, \@cur_vec, $num_genomes, "", 1, 1, $depth_limit, 0)) != $num_genomes) {
	    die "ERROR: all genomes should be accounted for in depth first search of the inserts at $start_end instead of $output_counts of $num_genomes\n";
	}
	if (($output_counts = &dfs_fgi_ins ($fGI_specified_file, $stop_clus, $start_clus, $stop_end, \@cur_vec, $num_genomes, "", 0, 1, $depth_limit, 0)) != $num_genomes) {
	    die "ERROR: all genomes should be accounted for in depth first search of the inserts at $stop_end instead of $output_counts of $num_genomes\n";
	}
    }
    close($inputinsertfile);
    close($fGI_specified_file);
    return;
}

sub genome_tags_string {

    my ($genome_vec) = @_;
    my $out_string = "";

    for my $index (0 .. ($num_genomes - 1)) {
	if ($genome_vec->[$index]) {
	    $out_string .= $genome_tags[$index] . ",";
	}
    }
    chop($out_string);
    
    return $out_string;
}

sub genome_groups_string {

    my ($genome_vec) = @_;
    my $out_string = "";
    my %local_group_counts = ();

    for my $group (@groups) {
	$local_group_counts{$group} = 0;
    }

    for my $index (0 .. ($num_genomes - 1)) {
	if ($genome_vec->[$index]) {
	    $local_group_counts{$tag_group{$genome_tags[$index]}}++;
	}
    }

    for my $group (@groups) {
	if ($local_group_counts{$group} > 0) {
	    my $rounded = sprintf "%.1f", ((100 * $local_group_counts{$group}) / $global_group_counts{$group});
	    $out_string .= $group . '(' . $rounded . '%),';
	}
    }
    chop($out_string);
    
    return $out_string;
}

sub dfs_fgi_ins {

    no warnings 'recursion';
    
    my ($out_file, $start_clus, $stop_clus, $cur_end, $cur_vec, $number, $clusters, $always_print, $limited, $depth, $length) = @_;
    $length++;
    my $output_counts = 0;
    (my $clus, my $end_type) = split('_', $cur_end);

    print STDERR "dfs_fgi_ins:$start_clus:$stop_clus:$cur_end\t$number\t$clusters\t$always_print\t@{ $cur_vec }\n" if ($DEBUG);
    if ($clus eq $stop_clus) {
	if ($high_clus_present[$clus]) {
	    $clusters = $clusters . ":STOP_CORE" . $clus;
	} else {
	    $clusters = $clusters . ":STOP_" . $clus;
	}
	if ($end_type eq "3") {
	    $clusters = $clusters . "+";
	} else {
	    $clusters = $clusters . "-";
	}
	if ($always_print) {
	    my $genomes;
	    if (($comma_tags eq "1") || ($comma_tags eq "2")) {
		$genomes = &genome_tags_string($cur_vec);
	    } elsif ($comma_tags eq "3") {
		$genomes = &genome_groups_string($cur_vec);
	    } else {
		$genomes = join("\t", @{ $cur_vec });
	    }
	    print $out_file "$clusters\t$number\t$genomes\tlength:$length\n";
	}
	print STDERR "return1\n" if ($DEBUG);
	return $number;
    }

    if ($limited) { # when using a depth limit we do not stop at unexpected core nodes
	if ($depth <= 0) {
	    if ($high_clus_present[$clus]) {
		$clusters = $clusters . ":DEPTH_CORE" . $clus;
	    } else {
		$clusters = $clusters . ":DEPTH_" . $clus;
	    }
	    if ($end_type eq "3") {
		$clusters = $clusters . "+";
	    } else {
		$clusters = $clusters . "-";
	    }
	    my $genomes;
	    if (($comma_tags eq "1") || ($comma_tags eq "2")) {
		$genomes = &genome_tags_string($cur_vec);
	    } elsif ($comma_tags eq "3") {
		$genomes = &genome_groups_string($cur_vec);
	    } else {
		$genomes = join("\t", @{ $cur_vec });
	    }
	    print $out_file "$clusters\t$number\t$genomes\tlength:$length\n";
	    print STDERR "return2a\n" if ($DEBUG);
	    return $number;
	}
    } elsif ($high_clus_present[$clus] && ($clus ne $start_clus)) { # we have reached an unexpected "core" node indicating a rearrangement so stop
	$clusters = $clusters . ":U_CORE" . $clus;
	if ($end_type eq "3") {
	    $clusters = $clusters . "+";
	} else {
	    $clusters = $clusters . "-";
	}
	my $genomes;
	if (($comma_tags eq "1") || ($comma_tags eq "2")) {
	    $genomes = &genome_tags_string($cur_vec);
	} elsif ($comma_tags eq "3") {
	    $genomes = &genome_groups_string($cur_vec);
	} else {
	    $genomes = join("\t", @{ $cur_vec });
	}
	print $out_file "$clusters\t$number\t$genomes\tlength:$length\n";
	print STDERR "return2b\n" if ($DEBUG);
	return $number;
    }
    
    if (!defined $adj_vec{$cur_end}) { # edges petered out before reaching $stop_clus so we are done
	$clusters = $clusters . ":" . $clus;
	if ($end_type eq "3") {
	    $clusters = $clusters . "+";
	} else {
	    $clusters = $clusters . "-";
	}
	my $genomes;
	if (($comma_tags eq "1") || ($comma_tags eq "2")) {
	    $genomes = &genome_tags_string($cur_vec);
	} elsif ($comma_tags eq "3") {
	    $genomes = &genome_groups_string($cur_vec);
	} else {
	    $genomes = join("\t", @{ $cur_vec });
	}
	print $out_file "$clusters\t$number\t$genomes\tlength:$length\n";
	print STDERR "return3\n" if ($DEBUG);
	return $number;
    }

    my $sum_counts = 0;
    my @next_vec = [];
    my @leftover_vec = (@{ $cur_vec });
    for my $next_end ( keys %{ $adj_vec{$cur_end} } )  { # go through all edges from $cur_end
	my $count = 0;
	for my $index (0 .. ($num_genomes - 1)) {
	    $count += ($next_vec[$index] = $cur_vec->[$index] && $adj_vec{$cur_end}{$next_end}->[$index]);
	    $leftover_vec[$index] = $leftover_vec[$index] && (!($adj_vec{$cur_end}{$next_end}->[$index]));
	}
	print STDERR "@{ $adj_vec{$cur_end}{$next_end} }\n" if ($DEBUG);
	print STDERR "$next_end\t$count\t@next_vec\n" if ($DEBUG);
	if ($count) {
	    #($clus, $end_type) = split('_', $cur_end); this is done at the start of the subroutine
	    my $next_clusters;
	    if ($clus eq $start_clus) { #this should be the start of the fGI
		$next_clusters = $clusters . "START_";
	    } else {
		$next_clusters = $clusters . ":";
	    }
	    if ($high_clus_present[$clus]) {
		$next_clusters = $next_clusters . "CORE" . $clus;
	    } else {
		$next_clusters = $next_clusters . $clus;
	    }
	    if ($end_type eq "3") {
		$next_clusters = $next_clusters . "+";
	    } else {
		$next_clusters = $next_clusters . "-";
	    }
	    (my $next_clus, my $next_end_type) = split('_', $next_end);
	    if ($next_end_type eq "3") {
		$next_end = $next_clus . "_5";
	    } else {
		$next_end = $next_clus . "_3";
	    }
	    $output_counts += &dfs_fgi_ins($out_file, $start_clus, $stop_clus, $next_end, \@next_vec, $count, $next_clusters, $always_print, $limited, ($depth - 1), $length);
	}
	$sum_counts += $count;
    }
    if ($sum_counts < $number) { # edges petered out before reaching $end_finish so output this
	my $next_clusters;
	if ($clus eq $start_clus) { #this should be the start of the fGI
	    $next_clusters = $clusters . "START_";
	} else {
	    $next_clusters = $clusters . ":";
	}
	if ($high_clus_present[$clus]) {
	    $next_clusters = $next_clusters . "CORE" . $clus;
	} else {
	    $next_clusters = $next_clusters . $clus;
	}
	if ($end_type eq "3") {
	    $next_clusters = $next_clusters . "+";
	} else {
	    $next_clusters = $next_clusters . "-";
	}
	$number -= $sum_counts;
	my $genomes;
	if (($comma_tags eq "1") || ($comma_tags eq "2")) {
	    $genomes = &genome_tags_string(\@leftover_vec);
	} elsif ($comma_tags eq "3") {
	    $genomes = &genome_groups_string(\@leftover_vec);
	} else {
	    $genomes = join("\t", @leftover_vec);
	}
	print $out_file "$next_clusters\t$number\t$genomes\tlength:$length\n";
	$output_counts += $number;
    }

    print STDERR "return4\n" if ($DEBUG);
    return $output_counts;
}

sub print_index {

    my $fGI_index_file_name = $insert_file_name . ".index";
    my $fGI_index_file;
    unless (open ($fGI_index_file, ">", $fGI_index_file_name) )  {
	die ("Can't open file $fGI_index_file_name.\n");
    }
    for my $clus (1 .. $num_clusters) {
	print $fGI_index_file "$clus\t$assembly_num_clus[$clus]\n";
    }
}
########################################  M A I N  #####################################################
print STDERR "Reading cluster sizes\n";
&get_cluster_sizes;
print STDERR "Reading centroids\n";
&read_centroids;
print STDERR "Reading high adjacency matrix $high_adjvec_file_name\n";
$high_global_largest_edge = &read_adj_mat ($high_adjvec_file_name, \@high_clus_present, \%high_adj_mat, 0, \@high_clus_present, 1);
if (!$combine_high_low) {
    print STDERR "Processing adjacency matrix\n";
    &process_adj_mat ($high_global_largest_edge, \@high_clus_present, \%high_adj_mat, \@high_clus_start, \%high_adj_best, \%high_adj_mutual, \%high_chain_length, \%high_chain_gene_length, 1, 1);
    print STDERR "Processing clusters with no edges in the input\n";
    &process_not_present;
    print STDERR "Printing high consensus order with context\n";
    $global_assembly_num = &print_consensus ($attribute_file_name, \@high_clus_present, \%high_adj_mat, \@high_clus_start, \%high_adj_best, \%high_adj_mutual, \%high_chain_length, \%high_chain_gene_length, "", $global_assembly_num, 0);
} else {
    print STDERR "Reading genome tags\n";
    &get_genome_tags;
    print STDERR "Processing high adjacency matrix\n";
    &process_adj_mat ($high_global_largest_edge, \@high_clus_present, \%high_adj_mat, \@high_clus_start, \%high_adj_best, \%high_adj_mutual, \%high_chain_length, \%high_chain_gene_length, 0, 1);
    print STDERR "Printing high consensus order with context\n";
    $global_assembly_num = &print_consensus ($attribute_file_name, \@high_clus_present, \%high_adj_mat, \@high_clus_start, \%high_adj_best, \%high_adj_mutual, \%high_chain_length, \%high_chain_gene_length, "_Core", $global_assembly_num, 1);
    print STDERR "Reading low adjacency matrix $low_adjvec_file_name\n";
    %adj_vec = (); # reset adj_vec for the low adjvec file
    $low_global_largest_edge = &read_adj_mat ($low_adjvec_file_name, \@low_clus_present, \%low_adj_mat, 1, \@high_clus_present, 1);
    print STDERR "Processing low adjacency matrix\n";
    &process_adj_mat ($low_global_largest_edge, \@low_clus_present, \%low_adj_mat, \@low_clus_start, \%low_adj_best, \%low_adj_mutual, \%low_chain_length, \%low_chain_gene_length, 1, 0);
    print STDERR "Processing cross adjacency matrix\n";
    &process_adj_mat_cross (\%high_adj_mat, \%low_adj_mat, \@low_clus_present, \@low_clus_start, \%high_adj_best, \%low_adj_best, \%high_adj_mutual, \%low_adj_mutual);
    print STDERR "Processing clusters with no edges in the input\n";
    &process_not_present;
    print STDERR "Printing low consensus order with context\n";
    $global_assembly_num = &print_consensus ($aux_attribute_file_name, \@low_clus_present, \%low_adj_mat, \@low_clus_start, \%low_adj_best, \%low_adj_mutual, \%low_chain_length, \%low_chain_gene_length, "_fGI", $global_assembly_num, 0);
    &print_insertions ($attribute_file_name, \@low_clus_present, \@low_clus_start, \%low_adj_mat, \%low_adj_best, \%low_adj_mutual, \%low_chain_length, \%low_chain_gene_length);
    &print_index;
}
exit(0);
