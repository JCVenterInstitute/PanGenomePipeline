#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $help_text = "This program reads in a file with cluster sizes and a file with paralogous cluster IDs and outputs to standard out the set of cinlge copy core genes.

Input Flags:
-sizes - A tab delimited file with cluster ID in column 1 and cluster size in column two
-paralogs - A file of paralogous cluster IDs separatede by white space
-cutoff - Threshold between 0-100 for definition of what percentage of genomes a core gene must be present in
-help - Outputs this help text";

GetOptions('sizes=s' => \my $sizes,
	'paralogs=s' => \my $paralogs,
	'cutoff=s' => \my $cutoff,
	'help' => \my $help);
	
if($help){
    print("$help_text\n");
    exit;
}
	
if(!$sizes or !$paralogs){
    die("Error: One or more of the required file arguments are missing\n$help_text\n");
}

#Globals
my @core = ();           # index = cluster ID, value = size
my %paralog = ();         # key = cluster ID, value = number of paralogs
my $num_clusters = 0;    # number of clusters
my $num_genomes = 0;     # the maximu size of a cluster which usually will be the number of genomes in the pan-genome


#####################################################################################################
sub read_cluster_sizes    # For each cluster store the size which is the number of genomes the cluster members are present in
{
    my $line;
    my @fields = ();
    my $sizes_file;
    unless (open ($sizes_file, "<", $sizes) )  {
	die ("Cannot open file $sizes!\n");
    }
    while ($line = <$sizes_file>) {
	chomp $line;
	@fields = split(/\t/, $line);
	$core[$fields[0]] = $fields[1];
	if ($fields[1] > $num_genomes) {
	    $num_genomes = $fields[1];
	}
	$num_clusters++;
    }
    $cutoff = ($cutoff * $num_genomes) / 100;
    close($sizes_file);
    return;
}

#####################################################################################################
sub read_paralogs       # For each cluster store the  number of paralogous clusters
{
    my $line;
    my @fields = ();
    my $paralogs_file;
    unless (open ($paralogs_file, "<", $paralogs) )  {
	die ("Cannot open file $paralogs!\n");
    }
    while ($line = <$paralogs_file>) {
	chomp $line;
	@fields = split(/\s/, $line);
	foreach my $cluster (@fields) {
	    $paralog{$cluster} = scalar @fields;
	}
    }
    close($paralogs_file);
    return;
}

#####################################################################################################
sub write_single_cores {  # read in the contigs for a genome and add 100,000 bp or as much as possible to both ends of circular contigs

    foreach my $index (1 .. $#core) {
	if ((!defined $paralog{$index}) && ($core[$index] >= $cutoff)) {
	    print "$index\n";
	}
    }
    return;
}


{#main
    &read_cluster_sizes;
    &read_paralogs;
    &write_single_cores;
}
