#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.


use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $vector_path = '', my $cluster_weights;
my $core, my $clusters,  my $debug, my $help_flag;

my %core_clusters = ();
my %cluster_neighbors = ();

my @edges = ();
my @neighbors = ();

my $help= "The core neighbor finder program goes through the Pan-Genome Graph (PGG)
which is the 0_core_adjacency_vector.txt
file that is outputted by PanOCT to determine how often two clusters are
neighbors on either end. The program will also use a specified core cutoff (0-100)
value or a list of known singleton core clusters to determine which core cluster each
cluster is closest to on either end and how often it is closest to that
core cluster.

INPUT FLAGS:
-vector: The PGG (0_core_adjacency_matrix.txt file that PanOCT automatically generates) 
-weights: The cluster_weights.txt file that PanOCT automatically generates
(just needs the cluster size in column two and cluster number in column 1 of a tab delimited file)
-core: A value between 1-100 indicating the percentage of genomes a given
cluster needs to be in to be considered core (also needs the -weights flag).
-clusters: A file containing a newline separated list of cluster numbers that
are believed to be singleton core clusters.";

GetOptions('core=i' => \$core,                   # The percentage of genomes a cluster needs to be in to be considered core
	'clusters=s' => \ $clusters,          # A list of cluster numbers to use
	'vector=s' => \$vector_path,          # The location of the vector file
	'weights=s' => \$cluster_weights,     # The location of the cluster weights file
	'debug' => \$debug,                   # Output debug text
	'help' => \$help_flag);               # Print help text

if($help_flag){
    die("$help\n");
}
if(($core and $clusters) or ($core and !$cluster_weights) or ($clusters and $cluster_weights) or (!$core and $cluster_weights) or (!$core and !$cluster_weights and !$clusters)){
    die("You must provide only a core percentage and a cluster sizes/weights file, or a list of singleton core clusters\n\n$help\n");
}
if($core && $cluster_weights){
    print("!!!FINDING CORE CLUSTERS using $core threshold!!!\n---------------------------------------------------------------------------------\n");
		# Determine which clusters meet the user's specified core cluster requirement
		# by dividing the number of genomes the cluster is in by the total number of
		# genomes supplied.
		open(VECTOR, "<", "$vector_path");          # Open a file that keeps track of which edges occur in which genomes
		my $count = split(/\t/, <VECTOR>) - 1;    # Read a single line and tab separate to get the number of genomes from the initial run
		close(VECTOR);
		
		open(WEIGHTS, "<", "$cluster_weights");     # Open the weights folder
		while(my $line = <WEIGHTS>){              # Read through each line of the weights file
			$line =~ /^(\d+?)\t(\d+?)\t/;      # Get the first two tab separated files (first # = cluster id, second # = number of genomes it's in)
			my $ratio = (100 * $2) / $count;         # Get percentage of genomes the cluster is in and add to list if it meets core requirement
			if($debug){print("$1 - $2 | $core - $ratio\n");}
			if($ratio >= $core){
				if($debug){print("$1 -> $2 -> $core <= $ratio\n");}
				$core_clusters{$1} = $ratio;
			}
		}
		close(WEIGHTS);
} elsif($clusters){# Go through each value in the list of clusters and add each one to the internal list of clusters
		# open and read through each line of the list of clusters
		open(CLUSTERS, "$clusters");
		while(my $line = <CLUSTERS>){
			$line =~ /(\w+)/; # Get just the cluster id without new line separators (\r, \n, etc.)
			$core_clusters{$1} = 1;
		}
		close(CLUSTERS);
} else {
    die("You must provide only a core percentage and a cluster sizes/weights file, or a list of singleton core clusters\n\n$help\n");
}

my $num_genomes;

# Open the adjacency vector file and go through it line-by-line
open(VECTOR, "$vector_path");
while(my $line = <VECTOR>){
        chomp $line;
	my @values = split(/\t/, $line);
        my $edge_id = shift @values;
	$edge_id =~ /^\s*\((\d+?)_(\d+),(\d+?)_(\d+)\)/; #Look at the first column to get the start cluster ($1), it's prime end ($2), the end cluster ($3), and it's prime end ($4)
	my $cluster_start = $1;
	my $start_side = $2;
	my $start = "$1" . "_" . "$2";
	my $cluster_end = $3;
	my $end_side = $4;
	my $end = "$3" . "_" . "$4";
	$num_genomes = @values;
	my$prime_end_index = ($start_side eq "3" ? 0 : 1); #Determine which index the edge information will be applied to (3' = 0, 5' = 1;)
	
	# Go through each genome
	for(my $i = 0; $i <= $#values; $i++){
		# If the edge exists for genome #i
		if($values[$i] == 1){
			$edges[$i]{$start} = $end;
			if(!defined($neighbors[$cluster_start][$prime_end_index]{$end})){
				$neighbors[$cluster_start][$prime_end_index]{$end} = 1;
				#print("GENOME $i: CREATING EDGE BETWEEN $cluster_start $start_side' AND $cluster_end $end_side' -> 1\n");
			}else{
				$neighbors[$cluster_start][$prime_end_index]{$end}++;
				#print("GENOME $i: COUNTING EDGE BETWEEN $cluster_start $start_side AND $cluster_end $end_side -> $neighbors[$cluster_start][$prime_end_index]{$end}\n");
			}
		}
	}
}
close(VECTOR);

#Print debug information about the edges
if($debug){
	print("\n\n!!!EDGES TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i <= $#edges; $i++){
		print("\nGenome $i:\n");
		if(defined($edges[$i])) {
		foreach my $start_node (sort keys(%{$edges[$i]})){
			print("$start_node -> $edges[$i]{$start_node}\n");
		}
                }
	}
	
	print("\n\n!!!CLOSEST NEIGHBORS TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i <= $#neighbors; $i++){
		for(my $j = 0; $j <= 1; $j++){
			my $start_end = ($j == 0 ? "3" : "5");
			foreach my $end_node (sort keys(%{$neighbors[$i][$j]})){
				my $count = $neighbors[$i][$j]{$end_node};
				print("($i" . "_" . "$start_end,$end_node) -> $count\n");
			}
		}
	}
        print("\n\n");
}

# Calculate the nearest core node for each node
foreach my $core_cluster (keys(%core_clusters)){
	#print("\nPROCESSING CORE: $core_cluster\n");
	
	for(my $i = 0; $i < $num_genomes; $i++){
		#print("Genome $i:\n");
		
		my $start_node = $core_cluster . "_3";
		my $core_node = $start_node;
		my $end_node = $edges[$i]{$start_node};
                my $dist = 1;

                if (defined $edges[$i]{$start_node}) {
		while($end_node =~ /(\d+?)_(\d+?)/){
                        my $cluster = $1;
                        my $which_end = $2;
			my $index = ($which_end == 3 ? 2 : 3);
			my $end = ($which_end == 3 ? 5 : 3);
			
			if(!defined($neighbors[$cluster][$index]{$core_node})){
				$neighbors[$cluster][$index]{$core_node} = 1;
				$neighbors[$cluster][$index]{$core_node . "_d"} = $dist;
			}else{
				$neighbors[$cluster][$index]{$core_node}++;
				$neighbors[$cluster][$index]{$core_node . "_d"} += $dist;
			}
			#print("$cluster $which_end' -> $core_node\n");
			if ($core_clusters{$cluster}) {
                                last;
                        }
			$start_node = $cluster . "_" . $end;
			$end_node = $edges[$i]{$start_node};
                        $dist++;
			#print("$start_node connects to: $end_node\n");
		}
		}
		#print("Done with 3'\n\n");
		
		$start_node = $core_cluster . "_5";
		$core_node = $start_node;
		$end_node = $edges[$i]{$start_node};
                $dist = 1;
		
                if (defined $edges[$i]{$start_node}) {
		while($end_node =~ /(\d+?)_(\d+?)/){
                        my $cluster = $1;
                        my $which_end = $2;
			my $index = ($which_end == 3 ? 2 : 3);
			my $end = ($which_end == 3 ? 5 : 3);
			
			if(!defined($neighbors[$cluster][$index]{$core_node})){
				$neighbors[$cluster][$index]{$core_node} = 1;
				$neighbors[$cluster][$index]{$core_node . "_d"} = $dist;
			}else{
				$neighbors[$cluster][$index]{$core_node}++;
				$neighbors[$cluster][$index]{$core_node . "_d"} += $dist;
			}
			#print("$cluster $which_end' -> $core_node\n");
			if ($core_clusters{$cluster}) {
                                last;
                        }
			$start_node = $cluster . "_" . $end;
			$end_node = $edges[$i]{$start_node};
                        $dist++;
			#print("$start_node connects to: $end_node\n");
		}
		}
		#print("Done with 5'\n\n");
	}
	#print("\n\/\/\n");
}

#Print debug information about the edges
if($debug){
	print("\n\n!!!CLOSEST CORES TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i <= $#neighbors; $i += 1){
		for(my $j = 2; $j <= 3; $j++){
			my $start_prime_end = ($j == 2 ? "3" : "5");
			foreach my $core_neighbor(sort keys(%{$neighbors[$i][$j]})){
				my $count = $neighbors[$i][$j]{$core_neighbor};
				print("($i" . "_" . "$start_prime_end,$core_neighbor) -> $neighbors[$i][$j]{$core_neighbor}\n");
			}
		}
	}
}

open(OUT, ">", "core_neighbors");

print OUT "#Cluster	3' neighbors and count	5' neighbors and count	3' core neighbors and count	5' core neighbors and count\n";

for(my $cluster = 1; $cluster <= $#neighbors; $cluster++){
    my $three_prime_output = "";
    my $first = 1;
    
    if(defined $neighbors[$cluster][0]){
	my %three_prime_neighbors = %{$neighbors[$cluster][0]};
	foreach my $neighbor (sort{$three_prime_neighbors{$b} <=> $three_prime_neighbors{$a}} keys(%three_prime_neighbors)){
	    if(!($first)){
		$three_prime_output .= ",";
	    }else{
		$first = 0;
	    }
	    $three_prime_output .= "$neighbor " . $three_prime_neighbors{$neighbor};
	}
    }else{
	$three_prime_output = "NONE"
    }
    
    my $five_prime_output = "";
    $first = 1;
    
    if(defined $neighbors[$cluster][1]){
	my %five_prime_neighbors = %{$neighbors[$cluster][1]};
	foreach my $neighbor (sort{$five_prime_neighbors{$b} <=> $five_prime_neighbors{$a}} keys(%five_prime_neighbors)){
	    if(!$first){
		$five_prime_output .= ",";
	    }else{
		$first = 0;
	    }
	    $five_prime_output .= "$neighbor " . $five_prime_neighbors{$neighbor};
	}
    }else{
	$five_prime_output = "NONE"
    }
    
    my $three_prime_core_output = "";
    $first = 1;
    
    if(defined $neighbors[$cluster][2]){
	my %three_prime_core_neighbors = ();
	foreach my $poss (keys %{ $neighbors[$cluster][2] }) { # remove the distance entries
	    if ($poss !~ /.*_d$/) {
		$three_prime_core_neighbors{$poss} = $neighbors[$cluster][2]{$poss};
	    }
	}
	foreach my $neighbor (sort{$three_prime_core_neighbors{$b} <=> $three_prime_core_neighbors{$a}} keys(%three_prime_core_neighbors)){
	    if(!$first){
		$three_prime_core_output .= ",";
	    }else{
		$first = 0;
	    }
	    $three_prime_core_output .= "$neighbor " . $three_prime_core_neighbors{$neighbor} . " " . sprintf("%.0f", ($neighbors[$cluster][2]{$neighbor . "_d"} / $three_prime_core_neighbors{$neighbor}));
	}
    }else{
	$three_prime_core_output = "NONE"
    }
	
    my $five_prime_core_output = "";
    $first = 1;
    if($neighbors[$cluster][3]){
	my %five_prime_core_neighbors = ();
	foreach my $poss (keys %{ $neighbors[$cluster][3] }) { # remove the distance entries
	    if ($poss !~ /.*_d$/) {
		$five_prime_core_neighbors{$poss} = $neighbors[$cluster][3]{$poss};
	    }
	}
	foreach my $neighbor (sort{$five_prime_core_neighbors{$b} <=> $five_prime_core_neighbors{$a}} keys(%five_prime_core_neighbors)){
	    if(!$first){
		$five_prime_core_output .= ",";
	    }else{
		$first = 0;
	    }
	    $five_prime_core_output .= "$neighbor " . $five_prime_core_neighbors{$neighbor} . " " . sprintf("%.0f", ($neighbors[$cluster][3]{$neighbor . "_d"} / $five_prime_core_neighbors{$neighbor}));
	}
    }else{
	$five_prime_core_output = "NONE"
    }
	
    print OUT "$cluster\t$three_prime_output\t$five_prime_output\t$three_prime_core_output\t$five_prime_core_output\n";
}

close(OUT);
