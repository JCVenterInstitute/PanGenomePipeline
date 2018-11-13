#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.


use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $vector_path = '', my $cluster_weights;
my $directory, my $core, my $clusters, my $centroids, my $debug, my $help_flag;

my %core_clusters = ();
my %cluster_neighbors = ();

my @edges = ();
my @neighbors = ();

my $help= "The core neighbor finder program goes through the 0_core_adjacency_vector.txt
file that is outputted by PanOCT to determine how often two clusters are
neighbors on either prime end. The programwill also use a specified core
value or a list of known core clusters to determine which core cluster each
cluster is closest to on either prime end and how often it is closest to that
core cluster.

PANOCT INPUT FLAGS:
-vector: The 0_core_adjacency_matrix.txt file that PanOCT automatically 
generates (REQUIRED if -dir is not used)
-weights: The cluster_weights.txt file that PanOCT automatically generates
(REQUIRED only if -dir is not used but -core is)
-dir: The directory to the folder PanOCT creates (REQUIRED if neither -vector
 nor -weights (if needed) is given)

CORE CLUSTERS INPUT FLAGS (use only one of these):
-core: A value between 1-100 indicating the percentage of genomes a given
cluster needs to be in to be considered (also needs either the -weight or
the -dir flags to use).
-clusters: A file containing a newline separated list of cluster numbers that
are believed to be core.";

GetOptions('dir=s' => \$directory,        # The directory where the PanOCT file are located
	'core=i' => \$core,                   # The percentage of genomes a cluster needs to be in to be considered core
	'clusters=s' => \ $clusters,          # A list of cluster numbers to use
	'vector=s' => \$vector_path,          # The location of the vector file
	'weights=s' => \$cluster_weights,     # The location of the cluster weights file
	'debug' => \$debug,                   # Output debug text
	'help' => \$help_flag);               # Print help text
	

if($help_flag){
	die("$help\n");
}

=pod
if(!($directory) and (!($vector_path) or !($cluster_weights))){
	die("A directory to the outputted pan-genome files is required\n\n$help\n");
}

elsif($directory){
	#$vector_path = "$directory/0_core_adjacency_vector_new_alt.txt";
	$vector_path = "$directory/0_core_adjacency_vector.txt";
	$cluster_weights = "$directory/cluster_weights.txt";
}

if(!($core or $clusters or $centroids)){
	die("A core percentage, a list of clusters, or a list of centroids must be provided\n\n$help\n");
}
=cut
if(($core and $clusters) or ($core and $centroids) or ($clusters and $centroids)){
#if((!!$core + !!$clusters + !!$centroids) != 1){
	die("You must provide only a core percentage, a list of clusters, or a list of centroids\n\n$help\n");
}

else{
	if($core){print("!!!FINDING CORE CLUSTERS!!!\n---------------------------------------------------------------------------------\n");
	
		# Determine which clusters meet the user's specified core cluster requirement
		# by dividing the number of genomes the cluster is in by the total number of
		# genomes supplied.
		open(VECTOR, "$vector_path");          # Open a file that keeps track of which edges occur in which genomes
		my $count = split(/\t/, <VECTOR>) - 1;    # Read a single line and tab separate to get the number of genomes from the initial run
		close(VECTOR);
		
		open(WEIGHTS, "$cluster_weights");     # Open the weights folder
		
		while(my $line = <WEIGHTS>){              # Read through each line of the weights file
			$line =~ /^(\d+?)\t(\d+?)\t/;      # Get the first two tab separated files (first # = cluster id, second # = number of genomes it's in)
			my $ratio = 100 * (int($2) / $count);         # Get percentage of genomes the cluster is in and add to list if it meets core requirement
			if($debug){print("$1 - $2 | $core - $ratio\n");}
			if($ratio >= $core){
				if($debug){print("$1 -> $2 -> $core <= $ratio\n");}
				$core_clusters{$1} = $ratio;
			}
		}
		
		close(WEIGHTS);
	}
	
	# Go through each value in the list of clusters and add each one to the internal list of clusters
	elsif($clusters){
		# open and read through each line of the list of clusters
		open(CLUSTERS, "$clusters");
		while(my $line = <CLUSTERS>){
			$line =~ /(\w+)/; # Get just the cluster id without new line separators (\r, \n, etc.)
			$core_clusters{$1} = 1;
		}
		
		close(CLUSTERS);
	}
=pod
	# Go through the list of known centroids and find and track the clusters for them
	else{
		
		# open and read through each line of the list of clusters
		my %centroids = ();
		open(CENTROIDS, "$centroids");
		while(my $line = <CENTROIDS>){
			$line =~ /(\w+)/; # Get just the centroid name without new line separators (\r, \n, etc.)
			$centroids{$1} = 1;
		}
		close(CENTROIDS);
		
		# Go through each line of the centroids.fasta file to find the cluster id's for each given centroid
		open(FASTA, "$directory/centroids.fasta");
		while(my $line = <FASTA>){
			if($line =~ /^>centroid_(\d+?) (.+?) .*\n/){        #For each header line get the cluster id and the centroid representing the cluster id
				if($centroids{$2}){                             #If the centroid is one of the ones the user supplied, track it's corresponding cluster id
					$core_clusters{$1} = 1;
				}
			}
		}
	}
=cut
} 

my $count;

# Open the adjacency vector file and go through it line-by-line
open(VECTOR, "$vector_path");
while(my $line = <VECTOR>){
	my @values = split(/\t/, $line);
	$values[0] =~ /^\s*\((\d+?)_(\d+),(\d+?)_(\d+)\)/; #Look at the first column to get the start cluster ($1), it's prime end ($2), the end cluster ($3), and it's prime end ($4)
	my $cluster_start = $1;
	my $start_side = $2;
	my $start = "$1" . "_" . "$2";
	my $cluster_end = $3;
	my $end_side = $4;
	my $end = "$3" . "_" . "$4";
	$count = @values;
	
	
	my$prime_end_index = ($start_side eq "3" ? 0 : 1); #Determine which index the edge information will be applied to (3' = 0, 5' = 1;)
	
	# Go through each genome
	for(my $i = 1; $i < @values; $i++){
		# If the edge exists for genome #i
		if($values[$i] == 1){
			$edges[$i - 1]{$start} = $end;
			
			if(!($neighbors[$cluster_start][$prime_end_index]{$end})){
				$neighbors[$cluster_start][$prime_end_index]{$end} = 1;
				#print("GENOME $i: CREATING EDGE BETWEEN $cluster_start $start_side' AND $cluster_end $end_side' -> 1\n");
			}
			
			else{
				$neighbors[$cluster_start][$prime_end_index]{$end} += 1;
				#$edge_count = $neighbors[$cluster_start][$prime_end_index]{$end};
				#print("GENOME $i: COUNTING EDGE BETWEEN $cluster_start $start_side' AND $cluster_end $end_side' -> $edge_count\n");
			}
		}
	}
	
	#print("\n\/\/\n\n");
}

#Print debug information about the edges
if($debug){
	print("\n\n!!!EDGES TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i < @edges; $i++){
		print("\nGenome $i:\n");
		my %genome_hash = ();
		if(defined($edges[$i])) {%genome_hash = %{$edges[$i]}};
		foreach my $start_node(sort keys(%genome_hash)){
			my $end_node = $genome_hash{$start_node};
			print("$start_node -> $end_node\n");
		}
	}
	
	print("\n\n!!!CLOSEST NEIGHBORS TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i < @neighbors; $i += 1){
		for(my $j = 0; $j < 2; $j += 1){
			my $start_end = ($j == 0 ? "3" : "5");
			foreach my $end_node(sort keys(%{$neighbors[$i][$j]})){
				my $count = $neighbors[$i][$j]{$end_node};
				print("($i" . "_" . "$start_end,$end_node) -> $count\n");
			}
		}
	}
}

if($debug){print("\n\n");}

# Calculate the nearest core node for each node
foreach my $core_cluster(sort keys(%core_clusters)){
	#print("\nPROCESSING CORE: $core_cluster\n");
	
	for(my $i = 0; $i < $count; $i++){
		#print("Genome $i:\n");
		
		my $start_node = $core_cluster . "_3";
		my $core_node = $start_node;
		my $end_node = $edges[$i]{$start_node};
                my $dist = 1;
		
		while($end_node =~ /(\d+?)_(\d+?)/){
			my $index = ($2 == 3 ? 2 : 3);
			my $end = ($2 == 3 ? 5 : 3);
			
			if(!($neighbors[$1][$index]{$core_node})){
				$neighbors[$1][$index]{$core_node} = 1;
				$neighbors[$1][$index]{$core_node . "_d"} = $dist;
			}
			
			else{
				$neighbors[$1][$index]{$core_node} += 1;
				$neighbors[$1][$index]{$core_node . "_d"} += $dist;
			}
			
			#print("$1 $2' -> $core_node\n");

			if ($core_clusters{$1}) {
                                last;
                        }
			$start_node = $1 . "_" . $end;
			$end_node = $edges[$i]{$start_node};
                        $dist++;
			
			#print("$start_node connects to: $end_node\n");
			
		}
		
		#print("Done with 3'\n\n");
		
		$start_node = $core_cluster . "_5";
		$core_node = $start_node;
		$end_node = $edges[$i]{$start_node};
                $dist = 1;
		
		while($end_node =~ /(\d+?)_(\d+?)/){
			my $index = ($2 == 3 ? 2 : 3);
			my $end = ($2 == 3 ? 5 : 3);
			
			if(!($neighbors[$1][$index]{$core_node})){
				$neighbors[$1][$index]{$core_node} = 1;
				$neighbors[$1][$index]{$core_node . "_d"} = $dist;
			}
			
			else{
				$neighbors[$1][$index]{$core_node} += 1;
				$neighbors[$1][$index]{$core_node . "_d"} += $dist;
			}
			
			#print("$1 $2' -> $core_node\n");

			if ($core_clusters{$1}) {
                                last;
                        }
			$start_node = $1 . "_" . $end;
			$end_node = $edges[$i]{$start_node};
                        $dist++;
			
			#print("$start_node connects to: $end_node\n");
			
		}
		
		#print("Done with 5'\n\n");
	}
	
	#print("\n\/\/\n");
}

#Print debug information about the edges
if($debug){
	print("\n\n!!!CLOSEST CORES TEST!!!\n---------------------------------------------------------------------------------\n");
	for(my $i = 0; $i < @neighbors; $i += 1){
		for(my $j = 2; $j < 4; $j++){
			my $start_prime_end = ($j == 2 ? "3" : "5");
			foreach my $core_neighbor(sort keys(%{$neighbors[$i][$j]})){
				my $count = $neighbors[$i][$j]{$core_neighbor};
				print("($i" . "_" . "$start_prime_end,$core_neighbor) -> $count\n");
			}
		}
	}
}

open(OUT, ">./core_neighbors");

print OUT "#Cluster	3' neighbors and count	5' neighbors and count	3' core neighbors and count	5' core neighbors and count\n";

for(my $cluster = 1; $cluster < @neighbors; $cluster++){
	my $three_prime_output = "", my $first = 1;
	
	if($neighbors[$cluster][0]){
		my %three_prime_neighbors = %{$neighbors[$cluster][0]};
		foreach my $neighbor(sort{$three_prime_neighbors{$b} <=> $three_prime_neighbors{$a}} 
		keys(%three_prime_neighbors)){
			if(!($first)){
				$three_prime_output .= ",";
			}
			
			else{
				$first = 0;
			}
			
			$three_prime_output .= "$neighbor " . $three_prime_neighbors{$neighbor};
		}
		
		%three_prime_neighbors = ();
	}
	
	if(!($three_prime_output)){
		$three_prime_output = "NONE"
	}
	
	my $five_prime_output = "";
        $first = 1;
	
	if($neighbors[$cluster][1]){
		my %five_prime_neighbors = %{$neighbors[$cluster][1]};
		foreach my $neighbor(sort{$five_prime_neighbors{$b} <=> $five_prime_neighbors{$a}} 
		keys(%five_prime_neighbors)){
			if(!$first){
				$five_prime_output .= ",";
			}
			
			
			else{
				$first = 0;
			}
			
			$five_prime_output .= "$neighbor " . $five_prime_neighbors{$neighbor};
		}
		
		%five_prime_neighbors = ();
	}
	
	if(!($five_prime_output)){
		$five_prime_output = "NONE"
	}
	
	my $three_prime_core_output = "";
        $first = 1;
	if($neighbors[$cluster][2]){
	        my %three_prime_core_neighbors = ();
                foreach my $poss (keys %{ $neighbors[$cluster][2] }) {
                        if ($poss !~ /.*_d$/) {
                                $three_prime_core_neighbors{$poss} = $neighbors[$cluster][2]{$poss};
                        }
                }
	        foreach my $neighbor(sort{$three_prime_core_neighbors{$b} <=> $three_prime_core_neighbors{$a}} 
		keys(%three_prime_core_neighbors)){
			if(!$first){
				$three_prime_core_output .= ",";
			}
			
			else{
				$first = 0;
			}
			
			$three_prime_core_output .= "$neighbor " . $three_prime_core_neighbors{$neighbor} . " " . sprintf("%.0f", ($neighbors[$cluster][2]{$neighbor . "_d"} / $three_prime_core_neighbors{$neighbor}));
		}
		
		%three_prime_core_neighbors = ();
	}
	
	if(!($three_prime_core_output)){
		$three_prime_core_output = "NONE"
	}
	
	
	my $five_prime_core_output = "";
        $first = 1;
	if($neighbors[$cluster][3]){
	        my %five_prime_core_neighbors = ();
                foreach my $poss (keys %{ $neighbors[$cluster][3] }) {
                        if ($poss !~ /.*_d$/) {
                                $five_prime_core_neighbors{$poss} = $neighbors[$cluster][3]{$poss};
                        }
                }
	        foreach my $neighbor(sort{$five_prime_core_neighbors{$b} <=> $five_prime_core_neighbors{$a}}
	        keys(%five_prime_core_neighbors)){
		        if(!$first){
		   	         $five_prime_core_output .= ",";
		        }
		
		
		        else{
			        $first = 0;
		        }
		
			$five_prime_core_output .= "$neighbor " . $five_prime_core_neighbors{$neighbor} . " " . sprintf("%.0f", ($neighbors[$cluster][3]{$neighbor . "_d"} / $five_prime_core_neighbors{$neighbor}));
		}
		
		%five_prime_core_neighbors = ();
	}
	
	
	if(!($five_prime_core_output)){
		$five_prime_core_output = "NONE"
	}
	
	print OUT "$cluster\t$three_prime_output\t$five_prime_output\t$three_prime_core_output\t$five_prime_core_output\n";
}

close(OUT);
