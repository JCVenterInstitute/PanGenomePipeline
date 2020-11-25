#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in a 3 column tab delimitted file from the distances option file. The 3 columns are the first three from mash dist runs: name1, name2, distance. This program outputs to files in the current directory single-linkage clusters using the threshold. The files are named after the founding genome and have suffix .cluster.

Input Flags:
-threshold - used for column 3 to see if two genomes should be clustered <= (default 0.02)
-distance - file name of distances data
-help - outputs this help text";

GetOptions('threshold=s' => \my $threshold,
	   'distances=s' => \my $distance_file,
	   'help' => \my $help);
	
if(!$threshold) {
    $threshold = 0.02;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my %cluster_name = ();         # key = genome name, value is the name of the founding cluster - with adjustments when clusters are merged.
my %dist_sum = ();             # key = genome name, value is the sum of distances within the cluster this genome is in
my %medoid_dist_sum = ();      # key = cluster head name, value is the minimum sum of distances
my %medoid = ();               # key = cluster head name, value is the medoid genome name
my %cluster_size = ();         # key = cluster head name, value is the size of the cluster

#####################################################################################################
sub read_input    # Read in the 3 column input (name1, name2, distance)
{
    my $distance_handle;
    unless (open ($distance_handle, "<", $distance_file) )  {
	die ("cannot open file $distance_file!\n");
    }
    my $line;
    while ($line = <$distance_handle>) {
	(my $name1, my $name2, my $distance) = split(/\t/, $line, 3);
	if ($name1 eq $name2) {
	    $cluster_name{$name1} = $name1;
	    next;
	}
	if ($distance <= $threshold) {
	    if (!defined($cluster_name{$name1})) {
		if (!defined($cluster_name{$name2})) {
		    $cluster_name{$name1} = $cluster_name{$name2} = $name1;
		} else {
		    $cluster_name{$name1} = $cluster_name{$name2};
		}
	    } else {
		if (!defined($cluster_name{$name2})) {
		    $cluster_name{$name2} = $cluster_name{$name1};
		} else {
		    if ($cluster_name{$name1} ne $cluster_name{$name2}) {
			my $head1 = $cluster_name{$name1};
			while ($head1 ne $cluster_name{$head1}) {
			    $head1 = $cluster_name{$head1};
			}
			my $head2 = $cluster_name{$name2};
			while ($head2 ne $cluster_name{$head2}) {
			    $head2 = $cluster_name{$head2};
			}
			$cluster_name{$head2} = $head1;
		    }
		}
	    }
	}
    }
    close($distance_handle);
	
    return;
}

#####################################################################################################
sub write_cluster_files    # Write out the single-linkage cluster files with names of the head genome and suffix .cluster
{
    foreach my $name (keys %cluster_name)  { # go through all the nodes of the cluster trees
	my $head = $cluster_name{$name};
	while ($head ne $cluster_name{$head}) {
	    $head = $cluster_name{$head};
	}
	$cluster_name{$name} = $head;
	my $cluster_file = $head . ".cluster";
	my $cluster_handle;
	unless (open ($cluster_handle, ">>", $cluster_file) )  {
	    die ("cannot open file $cluster_file!\n");
	}
	print $cluster_handle "$name\n";
	close($cluster_handle);
    }
    my $medoid_file = $distance_file . ".medoids";
    my $medoid_handle;
    unless (open ($medoid_handle, ">", $medoid_file) )  {
	die ("cannot open file $medoid_file!\n");
    }
    my $distance_handle;
    unless (open ($distance_handle, "<", $distance_file) )  {
	die ("cannot reopen file $distance_file!\n");
    }
    my $line;
    while ($line = <$distance_handle>) {
	(my $name1, my $name2, my $distance) = split(/\t/, $line, 3);
	if ($cluster_name{$name1} eq $cluster_name{$name2}) { # distance for same genome to same genome is 0 so ok to include
	    if (!defined($dist_sum{$name1})) {
		$dist_sum{$name1} = $distance;
	    } else {
		$dist_sum{$name1} += $distance;
	    }
	    if (!defined($dist_sum{$name2})) {
		$dist_sum{$name2} = $distance;
	    } else {
		$dist_sum{$name2} += $distance;
	    }
	}
    }
    close($distance_handle);
    foreach my $name (keys %cluster_name)  { # go through all the nodes of the cluster trees
	print STDERR "$name\t$cluster_name{$name}\t$dist_sum{$name}\n";
	if ((!defined($medoid_dist_sum{$cluster_name{$name}})) || ($dist_sum{$name} < $medoid_dist_sum{$cluster_name{$name}})) {
	    $medoid_dist_sum{$cluster_name{$name}} = $dist_sum{$name};
	    $medoid{$cluster_name{$name}} = $name;
	}
	if (!defined($cluster_size{$cluster_name{$name}})) {
	    $cluster_size{$cluster_name{$name}} = 1;
	} else {
	    $cluster_size{$cluster_name{$name}}++;
	}
	print STDERR "$name\t$cluster_name{$name}\t$dist_sum{$name}\t$medoid{$cluster_name{$name}}\t$medoid_dist_sum{$cluster_name{$name}}\t$cluster_size{$cluster_name{$name}}\n";
    }
    my $cluster_num = 0;
    foreach my $medoid_name (keys %medoid)  { # go through all the medoids
	$cluster_num++;
	print $medoid_handle "$cluster_num\t$cluster_size{$medoid_name}\t$medoid_name\t$cluster_name{$medoid_name}\n";
    }
    close($medoid_handle);
	
    return;
}

{#main
    &read_input;
    &write_cluster_files;
}
