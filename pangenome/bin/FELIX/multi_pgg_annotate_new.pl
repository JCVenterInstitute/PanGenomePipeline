#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

#PanGenome Annotation (Vector Approach) [ATTEMPT 1]

# pgg annotation wrapper script
use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use File::Compare;

my $bin_directory = "/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/";
my $input_bin_directory = "";
my @genomes = ();
my $input_single_copy = "";
my $single_copy = "single_copy_clusters.txt";
my $core_neighbors = "core_neighbors"; # is the file the core neighbors is stored in
my $genome_list_path = "";
my $new_genomes = "";
my $attributes = "";
my $weights = "cluster_sizes.txt";
my $project = "8520";
my $paralogs = "";
my $pgg = "pgg.txt";                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
my $medoids = "medoids.fasta";
my $matchtable = "matchtable.txt";                                                      # [pangenome_dir]/matchtable.txt
my $id = 95;
my $debug = 0;
my $help = 0;
my $logfile = "multi_pgg_annotate_new.logfile";
my $topology_file = "topology.txt";
my $cwd = getcwd;
my %old_genomes = ();


GetOptions('genomes=s' => \ $genome_list_path,
	   'new_genomes=s' => \ $new_genomes,
	   'topology=s' => \ $topology_file,
	   'single_copy=s' => \ $input_single_copy,
	   'bin_directory=s' => \ $input_bin_directory,
	   'attributes=s' => \ $attributes,
	   'weights=s' => \ $weights,
	   'project=s' => \ $project,
	   'paralogs=s' => \ $paralogs,
	   'pgg=s' => \ $pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ $medoids,
	   'match=s' => \ $matchtable,                                                      # [pangenome_dir]/matchtable.txt
	   'id=i' => \ $id,
	   'strip_version' => \my $strip_version,
	   'help' => \ $help,
	   'debug' => \ $debug);
if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genomes=s' => \ genome_list_path,
	   'new_genomes=s' => \ new_genomes,
	   'topology=s' => \ topology_file,
	   'single_copy=s' => \ input_single_copy,
	   'bin_directory=s' => \ input_bin_directory,
	   'attributes=s' => \ attributes,
	   'weights=s' => \ weights,
	   'project=s' => \ project,
	   'paralogs=s' => \ paralogs,
	   'pgg=s' => \ pgg,                                                               # [pangenome_dir]/0_core_adjacency_vector.txt
	   'medoids=s' => \ medoids,
	   'match=s' => \ matchtable,                                                      # [pangenome_dir]/matchtable.txt
	   'id=i' => \ id,
	   'strip_version' => \ strip_version,
	   'help' => \ help,
	   'debug' => \ debug);
_EOB_
    exit(0);
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

if ($paralogs && $input_single_copy) {
    die "You can only specify a paralogs file or a single copy clusters file but not both!\n";
} elsif (!$paralogs && !$input_single_copy) {
    die "You must specify either a paralogs file or a single copy clusters file!\n";
} elsif ($input_single_copy) {
    `cp $input_single_copy $single_copy`;
}
if (substr($weights, 0, 1) ne "/") {
    $weights = $cwd . "/$weights";
}
if (substr($medoids, 0, 1) ne "/") {
    $medoids = $cwd . "/$medoids";
}
if (substr($pgg, 0, 1) ne "/") {
    $pgg = $cwd . "/$pgg";
}
if ($debug) {print STDERR "Parameters:\ngenomes: $genome_list_path\nnew_genomes: $new_genomes\nattributes: $attributes\nweights: $weights\nparalogs: $paralogs\npgg: $pgg\nmedoids: $medoids\nmatch: $matchtable\nid: $id\nsingle_copy_clusters: $single_copy\n";}
			
######################################COMPONENT PROGRAM PATHS################################
my $single_copy_path = "$bin_directory/single_copy_core.pl";
my $core_neighbor_path = "$bin_directory/core_neighbor_finder.pl";
my $medoid_blast_path = "$bin_directory/medoid_blast_search.pl";
my $pgg_annotate_path = "$bin_directory/pgg_annotate.pl";
my $pgg_multifasta_path = "$bin_directory/pgg_edge_multifasta.pl";
my $compute_path = "$bin_directory/compute_pgg_graph.pl";
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

sub launch_grid_job {
# Given a shell script, launch it via qsub.

    my ( $name, $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue, $job_array_max ) = @_;

    my $qsub_command = "qsub -V -P $project_code -o $stdoutdir -e $stderrdir -wd $working_dir -r n -N $name";
    $qsub_command .= " -l $queue" if $queue;
    $qsub_command .= " -t 1-$job_array_max" if $job_array_max;

    $qsub_command .= " $shell_script";

    my $response = `$qsub_command`;
    my $job_id;

    if ($response =~ (/Your job (\d+) \(.*\) has been submitted/) || $response =~ (/Your job-array (\d+)\./)) {

        $job_id = $1;

    } else {
        die "Problem submitting the job!: $response";
    }

    return $job_id;

}


sub wait_for_grid_jobs {
    # Given a hash of job ids wait until hash is reduced to number of jobs specified and return number of jobs; name is the job name
    
    my ( $name, $number, $job_ids ) = @_;
    my $size = scalar( keys %{$job_ids} );

    while ( $size > $number ) {
	sleep 60;
	my $response = `qacct -j $name 2>&1`;
	&parse_response( $response, $job_ids );
	$size = scalar( keys %{$job_ids} );
    }
    return ($size);
}


sub parse_response {
# NOT INTENDED TO BE CALLED DIRECTLY.
# Given a qacct response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $job_ids ) = @_;
    my @qacct_array = split ( /=+\n/, $response );
    if (scalar(@qacct_array) <= 1) {
	return; # jobs haven't hit the grid yet
    }
    shift @qacct_array; # get rid of empty record at beginning.

    for my $record ( @qacct_array ) {

        next if ( $record =~ /error: ignoring invalid entry in line/ );

        chomp $record;

        my @rec_array = split ( "\n", $record );

        for my $line (@rec_array) {

            $line =~ s/(.*\S)\s+$/$1/;
            my ( $key, $value ) = split ( /\s+/, $line, 2 );
	    if ($key eq "jobnumber") {
		if ( defined $job_ids->{$value} ) {
		    delete ( $job_ids->{$value} )
		}
	    }
	}
    }
    return;
}
                                                                                               # combine files 1,3,4,5 described in prep-files along with the columns we build in compute
#############################################################################################
sub do_core_list
# run single_copy_core.pl to generate input for pgg_annotate.pl
{
    if ($debug) {print STDERR "\n$single_copy_path -s $weights -p $paralogs -c $id > $single_copy\n";}
    `/usr/bin/time -o cpustats -v $single_copy_path -s $weights -p $paralogs -c $id > $single_copy`;
    `echo "***$single_copy_path***" >> overhead_cpustats`;
    `cat cpustats >> overhead_cpustats`;
    `rm cpustats`;
    &bash_error_check("$single_copy_path -s $weights -p $paralogs -c $id > $single_copy", $?, $!);
}
#############################################################################################
sub do_neighbors
# run core_neighbor_finder.pl to generate input for pgg_annotate.pl
{
    if ($debug) {print STDERR "\n$core_neighbor_path -v $pgg -cl $single_copy\n";}
    `/usr/bin/time -o cpustats -v $core_neighbor_path -v $pgg -cl $single_copy >& $logfile`;
    `echo "***$core_neighbor_path***" >> overhead_cpustats`;
    `cat cpustats >> overhead_cpustats`;
    `rm cpustats`;
    &bash_error_check("$core_neighbor_path -v $pgg -cl $single_copy >& $logfile", $?, $!);
}
#############################################################################################
sub read_old_genomes
# read in list of identifiers and genomes paths, store them 
{
    if ($debug) {print STDERR "Read old genomes list\n\n";}
    open(GENOMES, "<", "$genome_list_path");
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
    if ($debug) {print STDERR "Load new genomes list\n\n";}
    open(GENOMES, "<", "$new_genomes");
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
######################################################################################################################################################################
sub read_topology {

    unless (open (CIRCFILE, "<", "$topology_file") )  {
	die ("ERROR: can not open contig topology file $topology_file.\n");
    }
    my $cur_tag = "";
    while (<CIRCFILE>) {
	my $tag = "";
	my $asmbl_id = "";
	my $type = "";

	($tag, $asmbl_id, $type) = split(/\t/, $_);  # split the scalar $line on tab
	if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	    die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the contig topology file $topology_file.\nLine:\n$_\n");
	}
	$cur_tag = $tag;
	
	unless (open (TOPFILE, ">", $cur_tag . "_topology.txt") )  {
	    die ("ERROR: can not open contig topology file $cur_tag" . "_topology.txt.\n");
	}
	print TOPFILE $_;
	last;
    }
    while (<CIRCFILE>) {
	my $tag = "";
	my $asmbl_id = "";
	my $type = "";

	($tag, $asmbl_id, $type) = split(/\t/, $_);  # split the scalar $line on tab
	if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	    die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the contig topology file $topology_file.\nLine:\n$_\n");
	}
	if ($tag ne $cur_tag) {
	    close (TOPFILE);
	    $cur_tag = $tag;
	    unless (open (TOPFILE, ">", $cur_tag . "_topology.txt") )  {
		die ("ERROR: can not open contig topology file $cur_tag" . "_topology.txt.\n");
	    }
	}
	print TOPFILE $_;
    }
    close (TOPFILE);
    close (CIRCFILE);
    return;
}
#############################################################################################
sub compute
# build matchtable, pgg, and attribute files
# also, build two single column files by counting lines in "new" files, corresponding to uniq_clus and uniq_edge
# these files will be used to rebuild the cluster_stats.txt file
{
    if ($debug) {print STDERR "Starting compute\n\n";}
    my $job_name = "cpgg_" . $$; #use a common job name so that qacct can access all of them together
    my %job_ids = ();
    my $num_jobs = 0;
    my $duplicate;
    if ($debug) {print STDERR "Starting grid genome processing\n\n";}
    for (my $j=0; $j <= $#genomes; $j++)
    {
	my $identifier = $genomes[$j][0];                                                 # get genome name
	my $genome_path = $genomes[$j][1];                                                # get genome path
	my $cpu_name = $cwd . "/$identifier" . "_cpu_total_stats";
	if (substr($genome_path, 0, 1) ne "/") {
	    $genome_path = $cwd . "/$genome_path";
	}
	if (defined $old_genomes{$identifier}) {
	    $duplicate = 1;
	} else {
	    $duplicate = 0;
	}
	my $shell_script = "/usr/bin/time -o $cpu_name -v $compute_path -duplicate $duplicate -name $identifier -genome $genome_path -weights $weights -medoids $medoids -pgg $pgg -debug";
	if ($strip_version) {
	    $shell_script .= " -strip_version";
	}
	my $stdoutfile = $cwd . "/" . $identifier . "_stdout";
	my $stderrfile = $cwd . "/" . $identifier . "_stderr";
	my $working_dir = $cwd . "/TMP_" . $identifier;
	my $queue = "himem";
	my $match_name = "$identifier" . "_match.col";
	my $pgg_name = "$identifier" . "_pgg.col";
	my $att_name = "$identifier" . "_attributes.txt";
	my $topology_name = "$identifier" . "_topology.txt";
	if ((-e $match_name) && (-e $pgg_name) && (-e $att_name)){
	    next; #we have already annotated this genome in a previous aborted run
	}
	`mkdir $working_dir`;
	`ln $topology_name $single_copy $core_neighbors $working_dir`;
	`cp $attributes $working_dir/combined.att`; 
	`cp $topology_file $working_dir/full_topology.txt`; 
	`cp $matchtable  $working_dir/matchtable.col`; 
	`cp $pgg $working_dir/pgg.col`; 
	`cp $genome_list_path $working_dir/combined_genome_list`;
	if ($debug) {print STDERR "\nidentifier: $identifier \t path: $genome_path\n\n";}
	if ($debug) {print STDERR "qsub $shell_script\n";}
	$job_ids{&launch_grid_job($job_name, $project, $working_dir, $shell_script, $stdoutfile, $stderrfile, $queue)} = 1;
	$num_jobs++;
	if ($num_jobs >= 50) {
	    $num_jobs = &wait_for_grid_jobs($job_name, 40, \%job_ids);
	}
    }
    &wait_for_grid_jobs($job_name, 0, \%job_ids);
    `rm -r TMP_*`;
    if ($debug) {print STDERR "removed TMP directories\n";}
    
    $num_jobs = 0;
    for (my $j=0; $j <= $#genomes; $j++)
    {
	my $identifier = $genomes[$j][0];                                                 # get genome name
	my $genome_path = $genomes[$j][1];                                                # get genome path
	my $match_name = "$identifier" . "_match.col";
	my $pgg_name = "$identifier" . "_pgg.col";
	my $att_name = "$identifier" . "_attributes.txt";
	if (!(-e $match_name) || !(-e $pgg_name) || !(-e $att_name)){
	    $num_jobs++;
	}
    }
    if ($debug) {print STDERR "$num_jobs FAILED resubmitting\n";}
    if ($num_jobs > 50) {
	die "Too many grid jobs failed $num_jobs\n";
    } elsif ($num_jobs > 0) {
	for (my $k=0; $k <= 2; $k++){ #try a maximum of 3 times on failed jobs
	    if ($debug) {print STDERR "Resubmit $num_jobs jobs Iteration $k\n";}
	    %job_ids = ();
	    $num_jobs = 0;
	    for (my $j=0; $j <= $#genomes; $j++)
	    {
		my $identifier = $genomes[$j][0];                                                 # get genome name
		my $genome_path = $genomes[$j][1];                                                # get genome path
		my $cpu_name = $cwd . "/$identifier" . "_cpu_total_stats";
		if (substr($genome_path, 0, 1) ne "/") {
		    $genome_path = $cwd . "/$genome_path";
		}
		my $match_name = "$identifier" . "_match.col";
		my $pgg_name = "$identifier" . "_pgg.col";
		my $att_name = "$identifier" . "_attributes.txt";
		if (!(-e $match_name) || !(-e $pgg_name) || !(-e $att_name)){
		    if (defined $old_genomes{$identifier}) {
			$duplicate = 1;
		    } else {
			$duplicate = 0;
		    }
		    my $shell_script = "/usr/bin/time -o $cpu_name -v $compute_path -duplicate $duplicate -name $identifier -genome $genome_path -weights $weights -medoids $medoids -pgg $pgg -debug";
		    if ($strip_version) {
			$shell_script .= " -strip_version";
		    }
		    my $stdoutfile = $cwd . "/" . $identifier . "_stdout";
		    my $stderrfile = $cwd . "/" . $identifier . "_stderr";
		    my $working_dir = $cwd . "/TMP_" . $identifier;
		    my $topology_name = "$identifier" . "_topology.txt";
		    my $queue = "himem";
		    `mkdir $working_dir`;
		    `ln $topology_name $single_copy $core_neighbors $working_dir`;
		    `cp $attributes $working_dir/combined.att`; 
		    `cp $topology_file $working_dir/full_topology.txt`; 
		    `cp $matchtable  $working_dir/matchtable.col`; 
		    `cp $pgg $working_dir/pgg.col`; 
		    `cp $genome_list_path $working_dir/combined_genome_list`;
		    if ($debug) {print STDERR "\nidentifier: $identifier \t path: $genome_path\n\n";}
		    if ($debug) {print STDERR "resubmit qsub $shell_script\n";}
		    $job_ids{&launch_grid_job($job_name, $project, $working_dir, $shell_script, $stdoutfile, $stderrfile, $queue)} = 1;
		    $num_jobs++;
		}
	    }
	    if ($num_jobs == 0) {
		last; # no failed jobs
	    }
	    if ($debug) {print STDERR "$num_jobs relaunched\n";}
	    &wait_for_grid_jobs($job_name, 0, \%job_ids);
	    `rm -r TMP_*`;
	    if ($debug) {print STDERR "removed resubmitted TMP directories\n";}
	}
    }

    if ($debug) {print STDERR "Starting genome processing\n\n";}
    # print headers to columns that are new (currently gene_ANI, rearrange, and wgsANI)
    `echo "geneANI" > gene_ANI`;
    `echo "rearrange" > rearrange`;
    `echo "wgsANI" > wgs_ANI`;
    for (my $j=0; $j <= $#genomes; $j++)
    {
	if ($debug) {print STDERR "Genome $j\n\n";}
	my $identifier = $genomes[$j][0];                                                 # get genome name
	if (defined $old_genomes{$identifier}) {
	    $duplicate = 1;
	} else {
	    $duplicate = 0;
	}
	if ($debug) {print STDERR "Genome $identifier\n\n";}
	my $match_name = "$identifier" . "_match.col";
	my $pgg_name = "$identifier" . "_pgg.col";
	my $gene_ani_name = "$identifier" . "_geneANI.txt";
	my $rearrange_name = "$identifier" . "_rearrange.txt";
	my $wgs_ani_name = "$identifier" . "_wgsANI.txt";
	my $att_name = "$identifier" . "_attributes.txt";
	my $new_att_name = "$identifier" . "_attributes_new.txt";
	my $new_match_name = "$identifier" . "_match_new.col";
	my $new_pgg_name = "$identifier" . "_pgg_new.col";
	my $uniq_clus_name = "$identifier" . "_uniq_clus.txt";
	my $uniq_edge_name = "$identifier" . "_uniq_edge.txt";
	my $anomalies_name_genome = "$identifier" . "_anomalies.txt";
	my $stats_name = "$identifier" . "_cluster_stats.txt";
	if ($debug) {print STDERR "matchname: $match_name \t pggname: $pgg_name \n";}
	die ("$match_name doesn't exist \n") unless (-e $match_name);
	die ("$pgg_name doesn't exist \n") unless (-e $pgg_name);
	die ("$att_name doesn't exist \n") unless (-e $att_name);
	`wc -l < $uniq_clus_name > uniq_clus`;
	`wc -l < $uniq_edge_name > uniq_edge`;
	`wc -l < $gene_ani_name >> gene_ANI`;
	`wc -l < $rearrange_name >> rearrange`;
	`cat $wgs_ani_name >> wgs_ANI`;                                                    # we don't need to do a line-count here, we just copy over the entire one-line file
	`rm $match_name $pgg_name $wgs_ani_name $att_name $new_match_name $new_pgg_name $new_att_name`;
	`rm $gene_ani_name $rearrange_name $uniq_clus_name $uniq_edge_name`;
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
	`rm stats.tail stats.tail.col1 stats.tail.col345 stats.tail.col7plus uniq_clus uniq_edge $stats_name`;
    }
    if ($duplicate) {
	`paste PGG_stats.txt gene_ANI rearrange wgs_ANI | sed -e 's/_ReDoDuP\t/\t/' > tmp.PGG_stats.txt`;                             #add in all columns that contain their own header (new columns)
    } else {
	`paste PGG_stats.txt gene_ANI rearrange wgs_ANI > tmp.PGG_stats.txt`;                             #add in all columns that contain their own header (new columns)
    }
    `mv tmp.PGG_stats.txt PGG_stats.txt`;
    `rm core_neighbors $single_copy gene_ANI rearrange wgs_ANI *_topology.txt`;
}

############################################### main

{#main
    `echo "Starting" > $logfile`;
    if ($debug) {print STDERR "Starting ...\n\n";}
    if ($paralogs ne "") {
	&do_core_list; # run single_copy_core
    }
    &do_neighbors;                                                                                 # generate pgg_neighborhood data
    &read_old_genomes;
    my $genome_count = &load_genomes;                                                              # read genome list, and store # of genomes
    &read_topology;
    &compute;                                                                                      # for all genomes, run blast, run pgg_annotate, concatenate as we go using paste
}
