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
my $genome_list_path = "";
my $new_genomes = "";
my $attributes = "";
my $weights = "cluster_sizes.txt";
my $paralogs = "";
my $pgg = "pgg.txt";                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
my $medoids = "medoids.fasta";
my $matchtable = "matchtable.txt";                                                      # [pangenome_dir]/matchtable.txt
my $id = 95;
my $debug = 0;
my $help = 0;
my $logfile = "multi_pgg_annotate_new.logfile";
my %old_genomes = ();


GetOptions('genomes=s' => \ $genome_list_path,
	   'new_genomes=s' => \ $new_genomes,
	   'single_copy=s' => \ $single_copy,
	   'attributes=s' => \ $attributes,
	   'weights=s' => \ $weights,
	   'paralogs=s' => \ $paralogs,
	   'pgg=s' => \ $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ $medoids,
	   'match=s' => \ $matchtable,                                                      # [pangenome_dir]/matchtable.txt
	   'id=i' => \ $id,
	   'help' => \ $help,
	   'debug' => \ $debug);
if ($debug) {print "Parameters:\ngenomes: $genome_list_path\nnew_genomes: $new_genomes\nattributes: $attributes\nweights: $weights\nparalogs: $paralogs\npgg: $pgg\nmedoids: $medoids\nmatch: $matchtable\nid: $id\nsingle_copy_clusters: $single_copy\n";}			
			
######################################COMPONENT PROGRAM PATHS################################
my $single_copy_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/single_copy_core.pl';
my $core_neighbor_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/core_neighbor_finder.pl';
my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/medoid_blast_search.pl';
my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_annotate.pl';
my $pgg_multifasta_path = '/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/pgg_edge_multifasta.pl';
#my $single_copy_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/single_copy_core.pl';
#my $core_neighbor_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/core_neighbor_finder.pl';
#my $medoid_blast_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/medoid_blast_search.pl';
#my $pgg_annotate_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/pgg_annotate.pl';
#my $pgg_multifasta_path = '/usr/local/projdata/8520/projects/PANGENOME/testing/pgg_edge_multifasta.pl';
#my $single_copy_path = '/home/gsutton/FELIX/single_copy_core.pl';
#my $core_neighbor_path = '/home/gsutton/FELIX/core_neighbor_finder.pl';
#my $medoid_blast_path = '/home/gsutton/FELIX/medoid_blast_search.pl';
#my $pgg_annotate_path = '/home/gsutton/FELIX/pgg_annotate.pl';
#my $pgg_multifasta_path = '/home/gsutton/FELIX/pgg_edge_multifasta.pl';
#############################################################################################
`echo "Starting" > $logfile`;
if ($debug) {print "Starting ...\n\n";}
if ($paralogs ne "") {
    &do_core_list; # run single_copy_core
} else {
    `cp $single_copy single_copy_clusters.txt`;
}
$single_copy = "single_copy_clusters.txt";
&do_neighbors;                                                                                 # generate pgg_neighborhood data
&read_old_genomes;
my $genome_count = &load_genomes;                                                              # read genome list, and store # of genomes
&compute;                                                                                      # for all genomes, run blast, run pgg_annotate, concatenate as we go using paste
                                                                                               # combine files 1,3,4,5 described in prep-files along with the columns we build in compute
#############################################################################################
sub do_core_list
# run single_copy_core.pl to generate input for pgg_annotate.pl
{
    if ($debug) {print "\n$single_copy_path -s $weights -p $paralogs -c $id > $single_copy\n";}
    `/usr/bin/time -o cpustats -v $single_copy_path -s $weights -p $paralogs -c $id > $single_copy`;
    `cat cpustats >> separate_cpustats`;
}
#############################################################################################
sub do_neighbors
# run core_neighbor_finder.pl to generate input for pgg_annotate.pl
{
    if ($debug) {print "\n$core_neighbor_path -v $pgg -cl $single_copy\n";}
    `/usr/bin/time -o cpustats -v $core_neighbor_path -v $pgg -cl $single_copy >& $logfile`;
    `cat cpustats >> separate_cpustats`;
}
#############################################################################################
sub read_old_genomes
# read in list of identifiers and genomes paths, store them 
{
    if ($debug) {print "Read old genomes list\n\n";}
    open(GENOMES, "$genome_list_path");
    while (my $line = <GENOMES>)
    {
	chomp($line);                                                               # strip newline character
	my @split_line = split(/\t/, $line);                                        # split on tab
	$old_genomes{$split_line[0]} = 1;                                           # store old genome name in hash as a key
    }
    close(GENOMES);
    return;
}
#############################################################################################
sub load_genomes
# read in list of identifiers and genomes paths, store them 
{
    if ($debug) {print "Load new genomes list\n\n";}
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
    close(GENOMES);
    return $count;
}
#############################################################################################
sub compute
# built matchtable, pgg, and attribute files
# also, build two single column files by counting lines in "new" files, corresponding to uniq_clus and uniq_edge
# these files will be used to rebuild the cluster_stats.txt file
{
    if ($debug) {print "Starting compute\n\n";}
    my $duplicate;
    open(GENEANI, ">", "gene_ANI");
    open(REARRANGE, ">", "rearrange");
    # print headers to columns that are new (currently gene_ANI, rearrange, and wgsANI)
    print GENEANI "geneANI\n";
    print REARRANGE "rearrange\n";
    `echo "wgsANI" > wgs_ANI`;
    if ($debug) {print "Starting genome processing\n\n";}
    for (my $j=0; $j < @genomes; $j++)
    {
	if ($debug) {print "Genome $j\n\n";}
	`cp $attributes combined.att`; 
	`cp $matchtable  matchtable.col`; 
	`cp $pgg  pgg.col`; 
	`cp $genome_list_path combined_genome_list`;
	open(CLUS, ">", "uniq_clus");
	open(EDGE, ">", "uniq_edge");
	my $identifier = $genomes[$j][0];                                                 # get genome name
	my $old_identifier = $identifier;
	if (defined $old_genomes{$identifier}) {
	    $duplicate = 1;
	    $identifier .= "_ReDoDuP";
	} else {
	    $duplicate = 0;
	}
	if ($debug) {print "Genome $identifier $old_identifier\n\n";}
	my $genome_path = $genomes[$j][1];                                                # get genome path
	open(GLIST, ">>", "combined_genome_list");
	print GLIST "$identifier\t$genome_path\n";
	close(GLIST);
	if ($debug) {print "\nidentifier: $identifier \t path: $genome_path\n\n";}
	if ($debug) {print "$medoid_blast_path -m $medoids -g $genome_path -b $identifier.BLAST\n";}
	`/usr/bin/time -o cpustats -v $medoid_blast_path -m $medoids -g $genome_path -b $identifier.BLAST >& $logfile`;     # BLAST genome against medoids
	`cat cpustats >> separate_cpustats`;
	if ($debug) {print "\n$pgg_annotate_path -b $identifier.BLAST -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg\n";}
#	`/usr/bin/time -o cpustats -v $pgg_annotate_path -b $identifier.BLAST -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg -reann >& $logfile`;
	`/usr/bin/time -o cpustats -v $pgg_annotate_path -b $identifier.BLAST -cl $weights -me $medoids -g $genome_path -co $single_copy  -t $identifier -ro $identifier -n core_neighbors -pgg $pgg >& $logfile`;
	`cat cpustats >> separate_cpustats`;
	my $match_name = ("$identifier" . "_match.col");
	my $pgg_name = ("$identifier" . "_pgg.col");
	my $gene_ani_name = ("$identifier" . "_geneANI.txt");
	my $rearrange_name = ("$identifier" . "_rearrange.txt");
	my $wgs_ani_name = ("$identifier" . "_wgsANI.txt");
	my $att_name = ("$identifier" . "_attributes.txt");
	my $new_att_name = ("$identifier" . "_attributes_new.txt");
	my $new_match_name = ("$identifier" . "_match_new.col");
	my $new_pgg_name = ("$identifier" . "_pgg_new.col");
	my $uniq_clus_name = ("$identifier" . "_uniq_clus.txt");
	my $uniq_edge_name = ("$identifier" . "_uniq_edge.txt");
	my $anomalies_name = ("output/anomalies.txt");
	my $anomalies_name_genome = ("$identifier" . "_anomalies.txt");
	my $stats_name = "output/cluster_stats.txt";
	my $uniq_clus = `wc -l < $uniq_clus_name`;
	my $uniq_edge = `wc -l < $uniq_edge_name`;
	my $gene_ani = `wc -l < $gene_ani_name`;
	my $rearrange = `wc -l < $rearrange_name`;
	if ($debug) {print "matchname: $match_name \t pggname: $pgg_name \n";}
	die ("$match_name doesn't exist \n") unless (-e $match_name);
	die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
	die ("$att_name doesn't exist \n") unless (-e $att_name);
	print CLUS "$uniq_clus";
	print EDGE "$uniq_edge";
	print GENEANI "$gene_ani";
	print REARRANGE "$rearrange";
	`cat $wgs_ani_name >> wgs_ANI`;                                                    # we don't need to do a line-count here, we just copy over the entire one-line file
	`paste matchtable.col $match_name > tmp.matchtable.col`;
	`paste pgg.col $pgg_name > tmp.pgg.col`;
	`cat $att_name >> combined.att`; 
	`mv tmp.matchtable.col matchtable.col`;                                            # rename file
	`mv tmp.pgg.col pgg.col;`;                                                         # rename file
	`rm $match_name $pgg_name $wgs_ani_name $att_name $new_match_name $new_pgg_name $new_att_name $identifier.BLAST`;
	close(CLUS);
	close(EDGE);
	`mkdir output`;                                                                                # first time - create necessary directories for pgg_edge_multifasta
	`mkdir multifasta`;
	if ($duplicate) {
	    if ($debug) {print "\n$pgg_multifasta_path -I $old_identifier -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $identifier -S -s $single_copy\n";}
	    `/usr/bin/time -o cpustats -v $pgg_multifasta_path -I $old_identifier -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $identifier -S -s $single_copy >& $logfile`;    # run pgg edge multi_fasta
	    `cat cpustats >> separate_cpustats`;
	} else {
	    if ($debug) {print "\n$pgg_multifasta_path -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $identifier -S -s $single_copy\n";}
	    `/usr/bin/time -o cpustats -v $pgg_multifasta_path -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $identifier -S -s $single_copy >& $logfile`;    # run pgg edge multi_fasta
	    `cat cpustats >> separate_cpustats`;
	}
	`cat $anomalies_name > $anomalies_name_genome`;
	`cat $gene_ani_name $rearrange_name $uniq_clus_name $uniq_edge_name >> $anomalies_name_genome`;
	`rm combined.att pgg.col matchtable.col $gene_ani_name $rearrange_name $uniq_clus_name $uniq_edge_name`;
	# Divide cluster_stats.txt into 5 files 
	# 1: All lines before first new genome [stats.head]
	# 2: All lines corresponding to new genomes [stats.tail]    (this file will be used to generate the next 3)
	# 3: Column 1 of just the rows corresponding to the new genomes [stats.tail.col1]       (generate with `cut stats.tail -f1 > stats.tail.col1`)
	# 4: Columns 3,4,5 of just the rows corresponding to the new genomes [stats.tail.col345] (generate with `cut stats.tail -f3,4,5 > stats.tail.col345`)
	# 5: Columnes 7,8,9 of just the rows corresponding to the new genomes [stats.tail.col7plus] (generate with `cut stats.tail -f7- > stats.tail.col7plus`)
	# create the new cluster_stats.txt by combining all the relevant files
	if ($j == 0) {
	    `head -n 1 $stats_name > PGG_stats.txt`;                                                               # generate file 1
	}
	`tail -n 1 $stats_name > stats.tail`;                                                      # generate file 2
	`cut stats.tail -f1 > stats.tail.col1`;                                                           # generate file 3
	`cut stats.tail -f3,4,5 > stats.tail.col345`;                                                     # generate file 4
	`cut stats.tail -f7- > stats.tail.col7plus`;                                                      # generate file 5
	`paste stats.tail.col1  uniq_clus stats.tail.col345 uniq_edge stats.tail.col7plus >> PGG_stats.txt`;
	`rm stats.tail stats.tail.col1 stats.tail.col345 stats.tail.col7plus uniq_clus uniq_edge`;
	`rm -r output multifasta`;
    }
    close(GENEANI);
    close(REARRANGE);
    if ($duplicate) {
	`paste PGG_stats.txt gene_ANI rearrange wgs_ANI | sed -e 's/_ReDoDuP\t/\t/' > tmp.PGG_stats.txt`;                             #add in all columns that contain their own header (new columns)
    } else {
	`paste PGG_stats.txt gene_ANI rearrange wgs_ANI > tmp.PGG_stats.txt`;                             #add in all columns that contain their own header (new columns)
    }
    `mv tmp.PGG_stats.txt PGG_stats.txt`;
    `rm core_neighbors single_copy_clusters.txt combined_genome_list gene_ANI rearrange wgs_ANI cpustats`;
}
