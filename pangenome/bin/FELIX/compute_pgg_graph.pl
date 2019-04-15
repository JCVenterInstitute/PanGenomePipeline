#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use File::Compare;

my $single_copy = "single_copy_clusters.txt";
my $core_neighbors = "core_neighbors"; # is the file the core neighbors is stored in
my $genome_path = "";
my $genome_name = "";
my $weights = "cluster_sizes.txt";
my $pgg = "pgg.txt";                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
my $medoids = "medoids.fasta";
my $debug = 0;
my $help = 0;

GetOptions('genome=s' => \ $genome_path,
	   'weights=s' => \ $weights,
	   'name=s' => \ $genome_name,
	   'pgg=s' => \ $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ $medoids,
	   'debug' => \ $debug,
	   'help' => \ $help);
if ($debug) {print "Parameters:\ngenomes: $genome_path\nweights: $weights\npgg: $pgg\nmedoids: $medoids\n";}			
			
			
######################################COMPONENT PROGRAM PATHS################################
my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/medoid_blast_search.pl';
my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_annotate.pl';
#############################################################################################


#############################################################################################
sub compute
# For the genome, run BLAST, run pgg_annotate
{
    my $match_name = ("$genome_name" . "_match.col");
    my $pgg_name = ("$genome_name" . "_pgg.col");
    my $att_name = ("$genome_name" . "_attributes.txt");
    my $blast_name = ("$genome_name" . ".BLAST");
    
    if ($debug) {print "Starting compute ...\n\n";}
    if ($debug) {print "\ngenome_name: $genome_name \t path: $genome_path\n\n";}
    if ($debug) {print "perl $medoid_blast_path -m $medoids -g $genome_path -b $blast_name\n";}  # BLAST genome against medoids
    `perl $medoid_blast_path -m $medoids -g $genome_path -b $blast_name`;  # BLAST genome against medoids
    if ($debug) {print "\nperl $pgg_annotate_path -re -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $genome_name -ro $genome_name -n $core_neighbors -pgg $pgg\n";}
    `perl $pgg_annotate_path -re -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $genome_name -ro $genome_name -n $core_neighbors -pgg $pgg`; 
    if ($debug) {print "\nmatchname: $match_name \t pggname: $pgg_name \n";}
    die ("$match_name doesn't exist \n") unless (-e $match_name);
    die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
    die ("$att_name doesn't exist \n") unless (-e $att_name);
    `rm $blast_name`;
    my $all_files = $genome_name . "_*";
    `mv $all_files ..`;
    `rm $single_copy $core_neighbors`; # remove the links for these
    return;
}

if ($debug) {print "Starting ...\n\n";}
&compute;                                                                                      # for  genome, run blast, run pgg_annotate
exit(0);
