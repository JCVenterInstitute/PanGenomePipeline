#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use File::Compare;

my $bin_directory = "/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/";
my $input_bin_directory = "";
my $single_copy = "single_copy_clusters.txt";
my $core_neighbors = "core_neighbors"; # is the file the core neighbors is stored in
my $genome_path = "";
my $genome_name = "";
my $weights = "cluster_sizes.txt";
my $pgg = "pgg.txt";                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
my $medoids = "medoids.fasta";
my $debug = 0;
my $help = 0;
my $duplicate = 0;
my $dup_genome_name = "";
my $target_genome_name = "";
my $cwd = getcwd;

GetOptions('genome=s' => \ $genome_path,
	   'weights=s' => \ $weights,
	   'name=s' => \ $genome_name,
	   'pgg=s' => \ $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ $medoids,
	   'bin_directory=s' => \ $input_bin_directory,
	   'duplicate=i' => \ $duplicate,
	   'strip_version' => \my $strip_version,
	   'reannotate' => \my $reannotate,
	   'debug' => \ $debug,
	   'help' => \ $help);
if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genome=s' => \ genome_path,
	   'weights=s' => \ weights,
	   'name=s' => \ genome_name,
	   'pgg=s' => \ pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ medoids,
	   'bin_directory=s' => \ input_bin_directory,
	   'duplicate=i' => \ duplicate,
	   'strip_version' => \ strip_version,
	   'reannotate' => \ reannotate,
	   'debug' => \ debug,
	   'help' => \ help);
_EOB_
    exit(0);
}
if ($debug) {print STDERR "Parameters:\ngenomes: $genome_path\nweights: $weights\npgg: $pgg\nmedoids: $medoids\n";}			

if ($duplicate) {
    $dup_genome_name = $genome_name . "_ReDoDuP";
    $target_genome_name = $dup_genome_name;
} else {
    $target_genome_name = $genome_name;
}
			
if ($input_bin_directory) {
    if (-d $input_bin_directory) {
	if (substr($input_bin_directory, 0, 1) ne "/") {
	    $input_bin_directory = $cwd . "/$input_bin_directory";
	}
    } else {
	die "The specified bin directory: $input_bin_directory does not exist!\n";
    }
    $bin_directory = $input_bin_directory;
}
			
######################################COMPONENT PROGRAM PATHS################################
my $medoid_blast_path = "$bin_directory/medoid_blast_search.pl";
my $pgg_annotate_path = "$bin_directory/pgg_annotate.pl";
my $pgg_multifasta_path = "$bin_directory/pgg_edge_multifasta.pl";
#############################################################################################

sub bash_error_check {
    my ($command, $error, $message) = @_;
    if (!$error) {
	return;
    }
    print STDERR "$command FAILED\n";
    if ($error == -1) {
	printf STDERR "failed to execute code(%d): %s\n", $error >> 8, $message;
    } elsif ($error & 127) {
	printf STDERR "child died with code %d signal %d, %s coredump\n", $error >> 8, ($error & 127),  ($error & 128) ? 'with' : 'without';
    } else {
	printf STDERR "child exited with value %d\n", $error >> 8;
    }
    return;
}

#############################################################################################
sub compute
# For the genome, run BLAST, run pgg_annotate
{
    my $match_name = "$genome_name" . "_match.col";
    my $pgg_name = "$genome_name" . "_pgg.col";
    my $att_name = "$genome_name" . "_attributes.txt";
    my $topology_name = "$genome_name" . "_topology.txt";
    my $blast_name = "$genome_name" . ".BLAST";
    my $cpu_name = "$genome_name" . "_cpu_separate_stats";
    
    if ($debug) {print STDERR "Starting compute ...\n\n";}
    if ($debug) {print STDERR "\ngenome_name: $genome_name \t path: $genome_path\n\n";}
    if ($strip_version) {
	if ($debug) {print STDERR "/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -strip_version -topology $topology_name -m $medoids -g $genome_path -b $blast_name\n";}  # BLAST genome against medoids
	`/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -strip_version -topology $topology_name -m $medoids -g $genome_path -b $blast_name`;  # BLAST genome against medoids
	`echo "***$medoid_blast_path***" >> $cpu_name`;
	`cat tmp_cpu_stats >> $cpu_name`;
	`rm tmp_cpu_stats`;
	&bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -strip_version -topology $topology_name -m $medoids -g $genome_path -b $blast_name", $?, $!);
    } else {
	if ($debug) {print STDERR "/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -b $blast_name\n";}  # BLAST genome against medoids
	`/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -b $blast_name`;  # BLAST genome against medoids
	`echo "***$medoid_blast_path***" >> $cpu_name`;
	`cat tmp_cpu_stats >> $cpu_name`;
	`rm tmp_cpu_stats`;
	&bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -b $blast_name", $?, $!);
    }
    if ($strip_version) {
	if ($reannotate) {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg`;
	    `echo "***$pgg_annotate_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg", $?, $!);
	} else {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg`;
	    `echo "***$pgg_annotate_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -strip_version -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg", $?, $!);
	} 
    } else {
	if ($reannotate) {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg`; 
	    `echo "***$pgg_annotate_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -re -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg", $?, $!);
	} else {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg`; 
	    `echo "***$pgg_annotate_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_annotate_path -topology $topology_name -b $blast_name -cl $weights -me $medoids -g $genome_path -co $single_copy -target $target_genome_name -ro $genome_name -n $core_neighbors -pgg $pgg", $?, $!);
	}
    }
    if ($debug) {print STDERR "\nmatchname: $match_name \t pggname: $pgg_name \n";}
    die ("$match_name doesn't exist \n") unless (-e $match_name);
    die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
    die ("$att_name doesn't exist \n") unless (-e $att_name);
    if (!$reannotate) {
	`mkdir output multifasta`;                                                                                # create necessary directories for pgg_edge_multifasta
	open(GLIST, ">>", "combined_genome_list");
	print GLIST "$target_genome_name\t$genome_path\n";
	close(GLIST);
	my $gene_ani_name = "$genome_name" . "_geneANI.txt";
	my $rearrange_name = "$genome_name" . "_rearrange.txt";
	my $wgs_ani_name = "$genome_name" . "_wgsANI.txt";
	my $new_att_name = "$genome_name" . "_attributes_new.txt";
	my $new_match_name = "$genome_name" . "_match_new.col";
	my $new_pgg_name = "$genome_name" . "_pgg_new.col";
	my $uniq_clus_name = "$genome_name" . "_uniq_clus.txt";
	my $uniq_edge_name = "$genome_name" . "_uniq_edge.txt";
	my $anomalies_name = "output/anomalies.txt";
	my $anomalies_name_genome = "$genome_name" . "_anomalies.txt";
	my $stats_name = "output/cluster_stats.txt";
	my $stats_name_genome = "$genome_name" . "_cluster_stats.txt";
	`paste matchtable.col $match_name > tmp.matchtable.col`;
	`paste pgg.col $pgg_name > tmp.pgg.col`;
	`cat $att_name >> combined.att`; 
	`cat $topology_name | sed -e 's/\t/_ReDoDuP\t/' >> full_topology.txt`; 
	`mv tmp.matchtable.col matchtable.col`;                                            # rename file
	`mv tmp.pgg.col pgg.col;`;                                                         # rename file
	if ($duplicate) {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy`;    # run pgg edge multi_fasta
	    `echo "***$pgg_multifasta_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy", $?, $!);
	} else {
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy`;    # run pgg edge multi_fasta
	    `echo "***$pgg_multifasta_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b multifasta -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy", $?, $!);
	}
	`mv $stats_name $stats_name_genome`;
	`cat $anomalies_name $gene_ani_name $rearrange_name $uniq_clus_name $uniq_edge_name > $anomalies_name_genome`;
	`rm combined.att pgg.col matchtable.col full_topology.txt combined_genome_list`;
    }
    `rm $blast_name`;
    my $all_files = $genome_name . "_*";
    `rm $topology_name $single_copy $core_neighbors`; # remove the links for these
    `mv $all_files ..`;
    return;
}

{#main
if ($debug) {print STDERR "Starting ...\n\n";}
&compute;                                                                                      # for  genome, run blast, run pgg_annotate
exit(0);
}
