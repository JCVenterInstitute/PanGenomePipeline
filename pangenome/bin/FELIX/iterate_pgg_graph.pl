#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

#PanGenome Annotation (Vector Approach) [ATTEMPT 1]

# pgg re-annotation wrapper script
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use File::Compare;

my @genomes = ();                                                                              
my $max_iterate = 1;                                                                           # if no iteration count is set, default is 1
my $single_copy = "single_copy_clusters.txt";

GetOptions('genomes=s' => \my $genome_list_path,
	   'weights=s' => \my $weights,
	   'paralogs=s' => \my $paralogs,
	   'pgg=s' => \my $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'centroids=s' => \my $centroids,
	   'match=s' => \my $matchtable,                                                      # [pangenome_dir]/matchtable.txt
	   'iterations=i' => $max_iterate,
	   'id=i' => \my $id,
	   'debug' => \my $debug,
	   'help' => \my $help);
			
			
######################################COMPONENT PROGRAM PATHS################################
my $single_copy_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/single_copy_core.pl';
my $core_neighbor_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/core_neighbor_finder.pl';
my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/medoid_blast_search.pl';
my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_annotate.pl';
my $pgg_multifasta_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_edge_multifasta.pl';
#############################################################################################

my $pgg_old = $pgg;



&do_core_list;                                                                                 # run single_copy_core
`mkdir output`;                                                                                # first time - create necessary directories for pgg_edge_multifasta
`mkdir multifasta`;
&load_genomes;                                                                                 # read genome list (we only want to do this once, not each iteration)
&compute;                                                                                      # for all genomes, run blast, run pgg_annotate, concatenate as we go using paste
                                                                       


# get first line of matchtable and pgg
# get core list
# generate pgg neighbor data
# read genome list
# run through genomes one at a time, doing blast then pgg_annotate, and adding match_table and pgg data to files
# check if there is a difference in pgg file
# if so, iterate the pgg_portion


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
# read in list of identifiers and genomes paths, store them so that this file doesn't need to be queried if the re-annotation is iterative
{
    open(GENOMES, "$genome_list_path");
    my $count = 0;
    while (my $line = <GENOMES>)
    {
	chomp($line);                                                               # strip newline character
	my @split_line = split(/\t/, $line);                                        # split on tab
	$genomes[$count][0] = $split_line[0];                                       # store identifier in 0
	$genomes[$count][1] = $split_line[1];                                       # store fasta path in 1
	$count++;                                                                   # increment counter
    }
	
}
#############################################################################################
sub compute
# go through all genomes, run BLAST, run pgg_annotate (building matchtable and pgg edges files as we go), then, see if there is a difference, and re-run if necessary
{
    for (my $i=1; $i <= $max_iterate; $i++)
    {
	if ($debug) {print "Iteration $i\n";}
	&do_neighbors;                                                                                 # run core_neighbor_finder
	`cut $matchtable -f1 > matchtable.col`;                                                        # get first line of existing matchtable file, use that as first column of new file
	`cut $pgg -f1 > pgg.col`;                
	for (my $j=0; $j <= $#genomes; $j++)
	{
	    my $identifier = $genomes[$j][0];                                                 # get genome name
	    my $genome_path = $genomes[$j][1];                                                # get genome path
	    if ($debug) {print "identifier: $identifier \t path: $genome_path \n\n";}
	    if ($debug) {print "perl $medoid_blast_path -m $centroids -g $genome_path -b $identifier.BLAST\n";}  # BLAST genome against centroids
	    `perl $medoid_blast_path -m $centroids -g $genome_path -b $identifier.BLAST`;  # BLAST genome against centroids
	    if ($debug) {print "perl $pgg_annotate_path -re -b $identifier.BLAST -cl $weights -me $centroids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg\n";}
	    `perl $pgg_annotate_path -re -b $identifier.BLAST -cl $weights -me $centroids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg`; 
	    
	    my $match_name = ("$identifier" . "_match.col");
	    my $pgg_name = ("$identifier" . "_pgg.col");
	    my $match_name_new = ("$identifier" . "_match_new.col");
	    my $pgg_name_new = ("$identifier" . "_pgg_new.col");
	    my $att_name = ("$identifier" . "_attributes.txt");
	    my $att_name_new = ("$identifier" . "_attributes_new.txt");
	    if ($debug) {print "matchname: $match_name \t pggname: $pgg_name \n";}
	    die ("$match_name doesn't exist \n") unless (-e $match_name);
	    die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
	    die ("$att_name doesn't exist \n") unless (-e $att_name);
	    `paste matchtable.col $match_name > tmp.matchtable.col`;                           # paste line frome matchtable
	    `paste pgg.col $pgg_name > tmp.pgg.col`;                                           # paste line from edges file
	    die ("tmp.matchtable.col is zero size \n") unless (-s "tmp.matchtable.col");
	    die ("tmp.pgg.col is zero size \n") unless (-s "tmp.pgg.col");
	    `mv tmp.matchtable.col matchtable.col`;                                            # rename file
	    `mv tmp.pgg.col pgg.col`;                                                         # rename file
	    if ($j==0)
	    {
		`cat $att_name > combined.att`;                                       # overwrite combined file from past iteration
	    } else 
	    {
		`cat $att_name >> combined.att`;                                      # add to combined file
	    }
	    # clean up
	    `rm $match_name $pgg_name $att_name $identifier.BLAST`;
	}
	if ($debug) {print "perl $pgg_multifasta_path -s $single_copy -B output -b multifasta -g $genome_list_path -m matchtable.col -a combined.att -p pgg.col -M $centroids -A -S -R\n";}    # run pgg edge multi_fasta
	`perl $pgg_multifasta_path -s $single_copy -B output -b multifasta -g $genome_list_path -m matchtable.col -a combined.att -p pgg.col -M $centroids -A -S -R`;    # run pgg edge multi_fasta
	$centroids = 'output/medoids.fasta';                                                           # after first iteration, we want to use the medoids.fasta file, not the supplied medoids file
	$single_copy = 'output/single_copy_clusters.txt';
	$matchtable = "output/matchtable.txt";                                                         # After first iteration, we need to update location of matchtable and pgg files
	$pgg = 'output/pgg.txt';
	
	if(compare("$pgg","$pgg_old") == 0)
	{
	    print "No differences found in last iteration - PGG is stable!\n";
	    last;
	}
	else
	{
	    `mv output/pgg.txt output/old.pgg.txt`;                                                   # set the current iteration as "old"
	    `mv output/matchtable.txt output/old.matchtable.txt`;                                     # set the current iteration as "old"
	    `mv combined.att output/old.combined.att`;                                                # save a copy of attributes
	    $pgg_old = 'output/old.pgg.txt';                                                          # after iteration 1, "old" is previous run, "new" is current run
	}
    }
}
