#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $cwd = getcwd;
my $help_text = "This program BLASTs a FASTA file of medoids against itself.
This work is done in a folder called M_BLAST_TMP, which is deleted when the 
program finishes. This program also takes a cluster sizes/weights file and a
percentage cutoff to determine a list of core genes. If core genes have high
identity matches to other core genes then they are not considered single copy
core genes which are output to standard out.

Input Flags:
-medoids - The nucleotide multiFASTA file containing the medoids centroids.fasta for PanOCT (required)
-sizes - A two column (tab delimitted) file first column is cluster number and second column is cluster size (required)
-number - The number of genomes (required)
-blast_directory - directory name for where blast executables are located - default is not to use a directory
-ld_load_directory - directory name for where blast libraries are located - default is not to use a directory
-help - Outputs this help text";

GetOptions('medoids=s' => \my $medoids,
	'sizes=s' => \my $cluster_sizes,
	'blast_directory=s' => \my $blast_directory,
	'ld_load_directory=s' => \my $ld_load_directory,
	'number=i' => \my $genome_number,
	'percentage=f' => \my $percentage,
	'help' => \my $help);
	
if($help){
    print("$help_text\n");
    exit;
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
	print("$help_text\n");
	exit;
    } # if no value for option s (seq_file), quit with help menu
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
	$blast_directory = 'export LD_LIBRARY_PATH=' . $ld_load_directory . ':$LD_LIBRARY_PATH; ' . $blast_directory;
    } else {
	print STDERR "Error with -ld_load_directory $ld_load_directory\n";
	print("$help_text\n");
	exit;
    }
} else {
    $ld_load_directory = "";
}
	
if(!$medoids or !$cluster_sizes or !$genome_number or !$percentage){
    die("Error: One or more of the required arguments are missing\n$help_text\n");
}

#Globals
my @core = ();     # 1 if core based on percentage 0 otherwise
my @single = ();   # 1 if core and single copy 0 otherwise
my $num_clusters = 0;
my $cutoff;        # cut off for core given number of genomes and percentage

#####################################################################################################
sub read_cluster_sizes    # For each cluster store the size which is the number of genomes the cluster members are present in
{
    my $line;
    my @fields = ();
    my $sizes_file;
    unless (open ($sizes_file, "<", $cluster_sizes) )  {
	die ("Cannot open file $cluster_sizes!\n");
    }
    $cutoff = ($percentage * $genome_number) / 100;
    $core[0] = 0;
    $single[0] = 0;
    while ($line = <$sizes_file>) {
	chomp $line;
	@fields = split(/\t/, $line);
	$single[$fields[0]] = $core[$fields[0]] = $fields[1] >= $cutoff ? 1 : 0;
	$num_clusters++;
    }
    close($sizes_file);
    return;
}


sub mod_blast { # eliminate blast matches to the added regions but keep one copy of overlapping matches
    my $blastin = shift (@_);
    my @btab_line = (); # array variable to store split btab lines
    my $qid = ""; # query id (cluster id)
    my $sid = ""; # subject id (contig from genome)
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $evalue; # blast evalue
    my $pid = ""; # percent identity
    my $score = ""; # BLAST bit score
    my $qlength = ""; # length of query sequence
    my $slength = ""; # length of subject sequence
    my $line = ""; # raw input line
    my $blast_in;
    
    unless (open($blast_in,"<", $blastin)) {
	die ("cannot open file $blastin!\n");
    }
    while ($line = <$blast_in>) {
	if ($line =~ /^#/)                                           # Skip header line
	{
	    next;
	}
	chomp $line;
	@btab_line = split(/\t/, $line);
	# ========================================================
	# btab output from NCBI blast+ blastn customized: -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore\"
	# column number Description
	# 0      Query_id
	# 1	 subject_id (Hit from db)
	# 2	 % Identity
	# 3	 start of alignment on query (5' nucleotide match in query)
	# 4	 end of alignment on query (3' nucleotide match in query)
	# 5	 query length
	# 6	 start of alignment on subject (5' for query)
	# 7	 end of alignment on subject (3' for query)
	# 8	 subject length
	# 9      e-value
	# 10     score (bits)
	# ========================================================
	$qid = $btab_line[0];
	$qid =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	$sid = $btab_line[1];
	$sid =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	if ($qid == $sid) {
	    next;
	}
	$pid = $btab_line[2];
	$qbegin = $btab_line[3];
	$qend = $btab_line[4];
	$qlength = $btab_line[5];
	$sbegin = $btab_line[6];
	$send = $btab_line[7];
	$slength = $btab_line[8];
	$evalue = $btab_line[9];
	$score = $btab_line[10];
	if ($core[$qid] && $core[$sid] && ((abs($qend - $qbegin) / $qlength) >= 0.5) && ((abs($send - $sbegin) / $slength) >= 0.5) && ($pid >= 95)) {
	    $single[$qid] = $single[$sid] = 0;
	}
    }
    close ($blast_in);
    foreach my $cluster (1 .. $num_clusters) {
	if ($single[$cluster]) {
	    print STDOUT "$cluster\n";
	}
    }
    return;
}

{#main
    &read_cluster_sizes;
    `mkdir M_BLAST_TMP`;
    `cp $medoids M_BLAST_TMP/temp_fasta.ftmp`;
    my $makeblastdb = $blast_directory . "makeblastdb";
    `$makeblastdb -in M_BLAST_TMP/temp_fasta.ftmp -dbtype nucl -out M_BLAST_TMP/temp_fasta.ftmp`;
    my $blastn = $blast_directory . "blastn";
    `$blastn -query $medoids -db M_BLAST_TMP/temp_fasta.ftmp -out M_BLAST_TMP/temp_results.ftmp -task blastn -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore\"`;
    &mod_blast("M_BLAST_TMP/temp_results.ftmp");
    `rm -r M_BLAST_TMP`;
}
