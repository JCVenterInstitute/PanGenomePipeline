#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

#PanGenome Annotation (Vector Approach) [ATTEMPT 1]

# pgg annotation wrapper script
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use File::Compare;

my @genomes = ();
my $single_copy = "single_copy_clusters.txt";
my $stats = "output/cluster_stats.txt";

GetOptions('genomes=s' => \my $genome_list_path,
	   'new_genomes=s' => \my $new_genomes,
	   'attributes=s' => \my $attributes,
	   'weights=s' => \my $weights,
	   'paralogs=s' => \my $paralogs,
	   'pgg=s' => \my $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \my $medoids,
	   'match=s' => \my $matchtable,                                                      # [pangenome_dir]/matchtable.txt
	   'id=i' => \my $id,
	   'help' => \my $help,
	   'debug' => \my $debug);
if ($debug) {print "Parameters:genomes: $genome_list_path\nnew_genomes: $new_genomes\nattributes: $attributes\nweights: $weights\nparalogs: $paralogs\npgg: $pgg\nmedoids: $medoids\nmatch: $matchtable\nid: $id\n";}			
			
######################################COMPONENT PROGRAM PATHS################################
#my $single_copy_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/single_copy_core.pl';
#my $core_neighbor_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/core_neighbor_finder.pl';
#my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/medoid_blast_search.pl';
#my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_annotate.pl';
#my $pgg_multifasta_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_edge_multifasta.pl';
my $single_copy_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/single_copy_core.pl';
my $core_neighbor_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/core_neighbor_finder.pl';
my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/medoid_blast_search.pl';
my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/pgg_annotate.pl';
my $pgg_multifasta_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/pgg_edge_multifasta.pl';
#############################################################################################
`cat $genome_list_path $new_genomes > combined_genome_list`;
&do_core_list;                                                                                 # run single_copy_core
&do_neighbors;                                                                                 # generate pgg_neighborhood data
`mkdir output`;                                                                                # first time - create necessary directories for pgg_edge_multifasta
`mkdir multifasta`;
my $genome_count = &load_genomes;                                                              # read genome list, and store # of genomes
&compute;                                                                                      # for all genomes, run blast, run pgg_annotate, concatenate as we go using paste
&cut_paste;                                                                                    # split cluster_stats.txt into 5 files, 4 of which will be combined with the columns we build into the final file
                                                                                               # combine files 1,3,4,5 described in prep-files along with the columns we build in compute
#############################################################################################
sub do_core_list
# run single_copy_core.pl to generate input for pgg_annotate.pl
{
    `perl $single_copy_path -s $weights -p $paralogs -c $id > $single_copy`;
}
#############################################################################################
sub do_neighbors
# run core_neighbor_finder.pl to generate input for pgg_annotate.pl
{
    `perl $core_neighbor_path -v $pgg -cl $single_copy`;
}
#############################################################################################
sub load_genomes
# read in list of identifiers and genomes paths, store them 
{
    open(GENOMES, "$new_genomes");
    my $count = 0;
    while (my $line = <GENOMES>)
    {
	chomp($line);                                                               # strip newline character
	my @split_line = split(/\t/, $line);                                        # split on tab
	$genomes[$count][0] = $split_line[0];                                       # store identifier in 0
	$genomes[$count][1] = $split_line[1];                                       # store fasta path in 1
	$count++;                                                                   # increment counter
    }
    return $count;
}
#############################################################################################
sub compute
# built matchtable, pgg, and attribute files
# also, build two single column files by counting lines in "new" files, corresponding to uniq_clus and uniq_edge
# these files will be used to rebuild the cluster_stats.txt file
{
    `cp $attributes combined.att`; 
    `cp $matchtable  matchtable.col`; 
    `cp $pgg  pgg.col`; 
    open(CLUS, ">", "uniq_clus");
    open(EDGE, ">", "uniq_edge");
    for (my $j=0; $j < @genomes; $j++)
    {
	my $identifier = $genomes[$j][0];                                                 # get genome name
	my $genome_path = $genomes[$j][1];                                                # get genome path
	`perl $medoid_blast_path -m $medoids -g $genome_path -b $identifier.BLAST`;     # BLAST genome against medoids
	`perl $pgg_annotate_path -b $identifier.BLAST -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg`;
	my $match_name = ("$identifier" . "_match.col");
	my $pgg_name = ("$identifier" . "_pgg.col");
	my $att_name = ("$identifier" . "_attributes.txt");
	my $new_match_name = ("$identifier" . "_match_new.col");
	my $new_pgg_name = ("$identifier" . "_pgg_new.col");
	my $uniq_clus = `wc -l < $new_match_name`;
	my $uniq_edge = `wc -l < $new_pgg_name`;
	die ("$match_name doesn't exist \n") unless (-e $match_name);
	die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
	die ("$att_name doesn't exist \n") unless (-e $att_name);
	print CLUS "$uniq_clus";
	print EDGE "$uniq_edge";
	`paste matchtable.col $match_name > tmp.matchtable.col`;
	`paste pgg.col $pgg_name > tmp.pgg.col`;
	`cat $att_name >> combined.att`; 
	`mv tmp.matchtable.col matchtable.col`;                                            # rename file
	`mv tmp.pgg.col pgg.col;`;                                                         # rename file
    }
    close(CLUS);
    close(EDGE);
    `perl $pgg_multifasta_path -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -A -S`;    # run pgg edge multi_fasta
}
#############################################################################################
sub cut_paste
# Divide cluster_stats.txt into 5 files 
# 1: All lines before first new genome [stats.head]
# 2: All lines corresponding to new genomes [stats.tail]    (this file will be used to generate the next 3)
# 3: Column 1 of just the rows corresponding to the new genomes [stats.tail.col1]       (generate with `cut stats.tail -f1 > stats.tail.col1`)
# 4: Columns 3,4,5 of just the rows corresponding to the new genomes [stats.tail.col345] (generate with `cut stats.tail -f3,4,5 > stats.tail.col345`)
# 5: Columnes 7,8,9 of just the rows corresponding to the new genomes [stats.tail.col7plus] (generate with `cut stats.tail -f7- > stats.tail.col7plus`)
# create the new cluster_stats.txt by combining all the relevant files
{
    `head -n 1 $stats > complete_cluster_stats.txt`;                                         # generate file 1
    `tail -n $genome_count $stats > stats.tail`;                       # generate file 2
    `cut stats.tail -f1 > stats.tail.col1`;                                                           # generate file 3
    `cut stats.tail -f3,4,5 > stats.tail.col345`;                                                     # generate file 4
    `cut stats.tail -f7- > stats.tail.col7plus`;                                                      # generate file 5
    `paste stats.tail.col1  uniq_clus stats.tail.col345 uniq_edge stats.tail.col7plus >> complete_cluster_stats.txt`;
}
