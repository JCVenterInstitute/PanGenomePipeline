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

my $commandline = join (" ", @ARGV);
print STDERR "$commandline\n";
my $blast_directory = "";
my $ld_load_directory = "";
my $blast_task = "blastn";
my $muscle_path = "";
my $rscript_path = "";
my $bin_directory = "/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/";
my $input_bin_directory = "";
my $single_copy = "single_copy_clusters.txt";
my $core_neighbors = "core_neighbors"; # is the file the core neighbors is stored in
my $genome_path = "";
my $genome_name = "";
my $weights = "cluster_sizes.txt";
my $pgg = "pgg.txt";
my $medoids = "medoids.fasta";
my $debug = 0;
my $help = 0;
my $duplicate = 0;
my $reannotate = 0;
my $strip_version = 0;
my $dup_genome_name = "";
my $target_genome_name = "";
my $cwd = getcwd;
my $multifastadir = "multifasta";
my $keep_divergent_alignments = "";
my $engdb = "";
my $nrdb = "";
my $pggdb = "";
my $no_MSA = 0;
my $no_filter_anomalies = 0;
my $less_memory = 0;
my $project = "8520";
my $combine_topology_ids = 0;
my $use_existing_db = 0;
my $codon_opt = 0;
my $soft_mask_id = "";
my $pggdb_topology_file = "";
my $use_local_disk;
my $orf_cluster_medoids = "";

GetOptions('genome=s' => \ $genome_path,
	   'weights=s' => \ $weights,
	   'name=s' => \ $genome_name,
	   'pgg=s' => \ $pgg,
	   'medoids=s' => \ $medoids,
	   'orf_cluster_medoids=s' => \ $orf_cluster_medoids,
	   'PGGdb_topology=s' => \ $pggdb_topology_file,
	   'pggdb=s' => \ $pggdb,
	   'engdb=s' => \ $engdb,
	   'nrdb=s' => \ $nrdb,
	   'project=s' => \ $project,
	   'bin_directory=s' => \ $input_bin_directory,
	   'blast_directory=s' => \ $blast_directory,
	   'ld_load_directory=s' => \ $ld_load_directory,
	   'blast_task=s' => \ $blast_task,
	   'soft_mask_id=s' => \ $soft_mask_id,
	   'muscle_path=s' => \ $muscle_path,
	   'rscript_path=s' => \ $rscript_path,
	   'multifastadir=s' => \ $multifastadir,
	   'alignments=s' => \ $keep_divergent_alignments,
	   'duplicate=i' => \ $duplicate,
	   'use_local_disk' => \ $use_local_disk,
	   'strip_version' => \ $strip_version,
	   'no_MSA' => \ $no_MSA,
	   'no_filter_anomalies' => \ $no_filter_anomalies,
	   'less_memory' => \ $less_memory,
	   'reannotate' => \ $reannotate,
	   'combine_topology_ids' => \ $combine_topology_ids,
	   'use_existing_db' => \ $use_existing_db,
	   'codon_opt' => \ $codon_opt,
	   'debug' => \ $debug,
	   'help' => \ $help);

if ($muscle_path) {
    if (-x $muscle_path) {
	if (substr($muscle_path, 0, 1) ne "/") {
	    $muscle_path = $cwd . "/$muscle_path";
	}
    } else {
	print STDERR "Error with -muscle_path $muscle_path\n";
	$help = 1;
    }
} else {
    $muscle_path = "";
}

if ($rscript_path) {
    if (-x $rscript_path) {
	if (substr($rscript_path, 0, 1) ne "/") {
	    $rscript_path = $cwd . "/$rscript_path";
	}
    } else {
	print STDERR "Error with -rscript_path $rscript_path\n";
	$help = 1;
    }
} else {
    $rscript_path = "";
}

if ($blast_directory) {
    if (-d $blast_directory) {
	if (substr($blast_directory, -1, 1) ne "/") {
	    $blast_directory .= "/";
	}
	if (substr($blast_directory, 0, 1) ne "/") {
	    $blast_directory = $cwd . "/$blast_directory";
	}
    } else {
	print STDERR "Error with -blast_directory $blast_directory\n";
	$help = 1;
    }
} else {
    $blast_directory = "";
}

if ($ld_load_directory) {
    if (-d $ld_load_directory) {
	if (substr($ld_load_directory, -1, 1) ne "/") {
	    $ld_load_directory .= "/";
	}
	if (substr($ld_load_directory, 0, 1) ne "/") {
	    $ld_load_directory = $cwd . "/$ld_load_directory";
	}
    } else {
	print STDERR "Error with -ld_load_directory $ld_load_directory\n";
	$help = 1;
    }
} else {
    $ld_load_directory = "";
}

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genome=s' => \ genome_path,
	   'weights=s' => \ weights,
	   'name=s' => \ genome_name,
	   'pgg=s' => \ pgg,
	   'medoids=s' => \ medoids,
	   'orf_cluster_medoids=s' => \ orf_cluster_medoids,
	   'PGGdb_topology=s' => \ pggdb_topology_file,
	   'pggdb=s' => \ pggdb,
	   'engdb=s' => \ engdb,
	   'nrdb=s' => \ nrdb,
	   'project=s' => \ project,
	   'bin_directory=s' => \ input_bin_directory,
	   'blast_directory=s' => \ blast_directory,
	   'ld_load_directory=s' => \ ld_load_directory,
	   'blast_task=s' => \ blast_task,
	   'soft_mask_id=s' => \ soft_mask_id,
	   'muscle_path=s' => \ muscle_path,
	   'rscript_path=s' => \ rscript_path,
	   'multifastadir=s' => \ multifastadir,
	   'alignments=s' => \ keep_divergent_alignments,
	   'duplicate=i' => \ duplicate,
	   'use_local_disk' => \ use_local_disk,
	   'strip_version' => \ strip_version,
	   'no_MSA' => \ no_MSA,
	   'no_filter_anomalies' => \ no_filter_anomalies,
	   'less_memory' => \ less_memory,
	   'reannotate' => \ reannotate,
	   'combine_topology_ids' => \ combine_topology_ids,
	   'use_existing_db' => \ use_existing_db,
	   'codon_opt' => \ codon_opt,
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
			
if (($orf_cluster_medoids ne "") && (substr($orf_cluster_medoids, 0, 1) ne "/")) {
    $orf_cluster_medoids = $cwd . "/$orf_cluster_medoids";
}
if ($keep_divergent_alignments) {
    if (-d $keep_divergent_alignments) {
	if (substr($keep_divergent_alignments, 0, 1) ne "/") {
	    $keep_divergent_alignments = $cwd . "/$keep_divergent_alignments";
	}
    } else {
	die "The specified alignments directory: $keep_divergent_alignments does not exist!\n";
    }
}

######################################COMPONENT PROGRAM PATHS################################
my $medoid_blast_path = "$bin_directory/medoid_blast_search.pl";
my $pgg_annotate_path = "$bin_directory/pgg_annotate.pl";
my $pgg_multifasta_path = "$bin_directory/pgg_edge_multifasta.pl -P $project ";
my $filter_anomalies_path = "$bin_directory/filter_anomalies.pl";
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
    $pgg_multifasta_path .= " -N $bin_directory ";
    if ($orf_cluster_medoids ne "") {
	$pgg_multifasta_path .= " -o $orf_cluster_medoids ";
    }
    if (!$no_MSA) {
	$pgg_multifasta_path .= " -l ";
    }
    if ($muscle_path ne "") {
	$pgg_multifasta_path .= " -C $muscle_path ";
    }
    if ($rscript_path ne "") {
	$pgg_multifasta_path .= " -r $rscript_path ";
    }
    if ($keep_divergent_alignments) {
	$pgg_multifasta_path .= " -k $keep_divergent_alignments ";
    }
    if ($codon_opt) {
	$pgg_multifasta_path .= " -O ";
    }	
    if ($strip_version) {
	$filter_anomalies_path .= " -strip_version ";
    }
    if ($combine_topology_ids) {
	$filter_anomalies_path .= " -combine_topology_ids ";
    }	
    if ($soft_mask_id) {
	$filter_anomalies_path .= " -soft_mask_id $soft_mask_id ";
    }	
    if ($use_existing_db) {
	$filter_anomalies_path .= " -use_existing_db ";
    }	
    if ($use_local_disk) {
	$medoid_blast_path .= " -use_local_disk ";
	$filter_anomalies_path .= " -use_local_disk ";
    }
    if ($blast_directory) {
	$medoid_blast_path .= " -blast_directory $blast_directory ";
	$filter_anomalies_path .= " -blast_directory $blast_directory ";
    }	
    if ($blast_task) {
	$medoid_blast_path .= " -blast_task $blast_task ";
	$filter_anomalies_path .= " -blast_task $blast_task ";
    }	
    if ($ld_load_directory) {
	$medoid_blast_path .= " -ld_load_directory $ld_load_directory ";
	$filter_anomalies_path .= " -ld_load_directory $ld_load_directory ";
    }	
    if ($strip_version) {
	$medoid_blast_path .= " -strip_version ";
    }	
    if ($debug) {print STDERR "/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -blastout $blast_name\n";}  # BLAST genome against medoids
    `/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -blastout $blast_name`;  # BLAST genome against medoids
    `echo "***$medoid_blast_path***" >> $cpu_name`;
    `cat tmp_cpu_stats >> $cpu_name`;
    `rm tmp_cpu_stats`;
    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $medoid_blast_path -topology $topology_name -m $medoids -g $genome_path -blastout $blast_name", $?, $!);
    if ($strip_version) {
	$pgg_annotate_path .= " -strip_version ";
    }	
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
	my $ce_sizes_name = "output/ce_sizes.txt";
	my $ce_sizes_name_genome = "$genome_name" . "_ce_sizes.txt";
	my $stats_name = "output/cluster_stats.txt";
	my $stats_name_genome = "$genome_name" . "_cluster_stats.txt";
	`paste matchtable.col $match_name > tmp.matchtable.col`;
	`paste pgg.col $pgg_name > tmp.pgg.col`;
	`cat $att_name >> combined.att`; 
	`mv tmp.matchtable.col matchtable.col`;                                            # rename file
	`mv tmp.pgg.col pgg.col;`;                                                         # rename file
	if ($strip_version) {
	    $pgg_multifasta_path .= " -V ";
	}
	if (!$no_filter_anomalies) {
	    if ($pggdb_topology_file eq "") {
		`cp full_topology.txt PGGdb_topology.txt`; # need this if doing filter_anomalies here
		$pggdb_topology_file = "PGGdb_topology.txt";
	    }
	}
	if ($duplicate) {
	    my $dup_topology_name = "dup_" . $topology_name;
	    `sed -e 's/\t/_ReDoDuP\t/' $topology_name > $dup_topology_name`; 
	    `cat $dup_topology_name >> full_topology.txt`; 
	    if ($less_memory) {
		$pgg_multifasta_path .= " -F -e $att_name -i $dup_topology_name ";
	    }
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy`;    # run pgg edge multi_fasta
	    `echo "***$pgg_multifasta_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    `rm $dup_topology_name`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -I $genome_name -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $dup_genome_name -S -s $single_copy", $?, $!);
	} else {
	    `cat $topology_name >> full_topology.txt`; 
	    if ($less_memory) {
		$pgg_multifasta_path .= " -F -e $att_name -i $topology_name ";
	    }
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy`;    # run pgg edge multi_fasta
	    `echo "***$pgg_multifasta_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $pgg_multifasta_path -T full_topology.txt -B output -b $multifastadir -g combined_genome_list -m matchtable.col -a combined.att -p pgg.col -t $genome_name -S -s $single_copy", $?, $!);
	}
	`mv $ce_sizes_name $ce_sizes_name_genome`;
	`cat $anomalies_name $gene_ani_name $rearrange_name $uniq_clus_name $uniq_edge_name > $anomalies_name_genome`;
	`rm combined.att pgg.col matchtable.col full_topology.txt combined_genome_list`;
	
	if (!$no_filter_anomalies) {
	    my $filter_genomes_name = "$genome_name" . "_filter_genomes.txt";
	    my $filter_features_name = "$genome_name" . "_FEATURES";
	    open(FGLIST, ">", $filter_genomes_name);
	    print FGLIST "$genome_name\t$genome_path\t$topology_name\t$anomalies_name_genome\n";
	    close(FGLIST);
	    if ($debug) {print STDERR "\n/usr/bin/time -o tmp_cpu_stats -v $filter_anomalies_path -bin_directory $bin_directory -PGGdb_topology $pggdb_topology_file -genomes $filter_genomes_name -engdb $engdb -nrdb $nrdb -pggdb $pggdb\n";}
	    `/usr/bin/time -o tmp_cpu_stats -v $filter_anomalies_path -bin_directory $bin_directory -PGGdb_topology $pggdb_topology_file -genomes $filter_genomes_name -engdb $engdb -nrdb $nrdb -pggdb $pggdb`;
	    `echo "***$filter_anomalies_path***" >> $cpu_name`;
	    `cat tmp_cpu_stats >> $cpu_name`;
	    `rm tmp_cpu_stats`;
	    &bash_error_check("/usr/bin/time -o tmp_cpu_stats -v $filter_anomalies_path -bin_directory $bin_directory -PGGdb_topology $pggdb_topology_file -genomes $filter_genomes_name -engdb $engdb -nrdb $nrdb -pggdb $pggdb", $?, $!);
	    `paste $stats_name $filter_features_name > $stats_name_genome`;
	    `rm $filter_features_name $filter_genomes_name`;
	} else {
	    `mv $stats_name $stats_name_genome`;
	}
	if ($codon_opt) {
	    my $codon_opt_name = "output/$genome_name" . "_codon_opt_mutation.txt";
	    my $codon_opt_name_genome = "$genome_name" . "_codon_opt_mutation.txt";
	    `mv $codon_opt_name $codon_opt_name_genome`;
	    $codon_opt_name = "output/$genome_name" . "_codon_opt_insertion.txt";
	    $codon_opt_name_genome = "$genome_name" . "_codon_opt_insertion.txt";
	    `mv $codon_opt_name $codon_opt_name_genome`;
	}
	if ($orf_cluster_medoids ne "") {
	    my $stop_codon_name = "output/$genome_name" . "_stop_codon.txt";
	    my $stop_codon_name_genome = "$genome_name" . "_stop_codon.txt";
	    `mv $stop_codon_name $stop_codon_name_genome`;
	}
	`rm -r output multifasta`;
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
