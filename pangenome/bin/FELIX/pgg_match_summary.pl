#!/usr/bin/env perl
#Copy (C) 2013  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Granger Sutton, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes
my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use warnings;
use Getopt::Std;
getopts ('Dhm:c:g:');
our ($opt_D,$opt_h,$opt_m,$opt_g,$opt_c);

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver1.0";
my $matchtable_file;
my $DEBUG;
my $genomes_file_name;
my @genome_array = ();
my $genome_number;
my $clusters_file;
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if (($opt_c) && (-s "$opt_c")) {$clusters_file = $opt_c;} else { print STDERR "Error with -c $opt_c\n"; &option_help; } # if no value for option c (clusters input file), quit with help menu
if (($opt_m) && (-s "$opt_m")) {$matchtable_file = $opt_m;} else { print STDERR "Error with -m $opt_m\n"; &option_help; } # if no value for option m (matchtable input file), quit with help menu
if (($opt_g) && (-s "$opt_g")) {$genomes_file_name = $opt_g;} else { print STDERR "Error with -g\n"; &option_help; } # if no value for option g (genome tagsnames input file), quit with help menu

my %include_cluster = ();      # key = cluster number, value is defined if cluster is to be included otherwise undefiend
my @genome_totals = ();        # index is order of genome in matchtable, value is number of the specified clusters in that genome

######################################################################################################################################################################

sub get_genomes {  # obtain list of genomes - must be in the same order as the matchtable columns - and the mulitfasta contigs file for the genomes
   
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: cannot open file $genomes_file_name\n");
    print STDERR "Order of genomes in $genomes_file_name with array index\n";
    my $target_found = 0;
    while (my $line1 = <$infile>)  {
	my $name = $line1;
	chomp $name;

	push (@genome_array, $name); # populate the genome_array in the order of the genome file
	print STDERR "$name\t$genome_number\n";
	$genome_totals[$genome_number] = 0;
	$genome_number++;
    }
    close($infile);
    print STDERR "$genome_number genomes\n\n";
}

sub process_matchtable {

    unless (open (TABLEFILE, "<", "$matchtable_file") )  {
	die ("ERROR: cannot open file $matchtable_file.\n");
    }
    my $cluster_num = 0;
    print "Cluster";
    foreach my $genome_tag (@genome_array) {
	print "\t$genome_tag";
    }
    print "\n";
    while (my $line = <TABLEFILE>) {
	chomp $line;
	my @feat_names = split(/\t/, $line);  # split the scalar $line on tab
	my $cluster_id = shift @feat_names;
	if (!defined $include_cluster{$cluster_id}) {
	    next; # skip clusters not in the clusters file
	}
	print "$line\n";
	my $genome_count = 0;
	foreach my $feat_name (@feat_names) {
	    if (($feat_name eq "----------") || ($feat_name eq "")) { #this is a placeholder and can be skipped
	    } else {
		$genome_totals[$genome_count]++;
	    }
	    $genome_count++;
	}
	$cluster_num++;
    }
    close (TABLEFILE);
    print "Total($cluster_num)";
    foreach my $total (@genome_totals) {
	print "\t$total";
    }
    print "\n";
    return;
}

sub read_clusters_file                                               # Read in the set of clusters which are single copy core, renumber and output
{
    unless (open(CLUSTER_CORES, "<", $clusters_file)) {
	die ("cannot open cores file: $clusters_file!\n");
    }
    while (my $line = <CLUSTER_CORES>) {
	chomp $line;
	$line =~ s/\s*//g; # remove all whitespace characters
	$line =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	$include_cluster{$line} = 1;
	#print "$line:$include_cluster{$line}\n";
    }
    close(CLUSTER_CORES);
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - Pan-genome Graph, multifasta cluster and edge files
Copy (C) 2018  The J. Craig Venter Institute (JCVI).  All rights reserved

License:   This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

Citation:  Derrick E. Fouts, Lauren Brinkac, Erin Beck, Jason Inman, and Granger Sutton (2011) "PanOCT: Automated Clustering of Orthologs using 
           Conserved Gene Neighborhood for Pan-Genomic Analysis of Bacterial Strains and Closely Related Species" Nucleic Acids Res. 2012 Dec;40(22).

  Usage: $prog <options>
Example: $prog -g panoct_genome_tag_file -m panoct_matchtable_file -c clusters_file > output_file
Version: $version
 Option:
     -h: print this help page
     -g: one genome_tag per line as used in panoct in same order as panoct
     -m: panoct matchtable file (input)
     -c: cluster file one cluster ID per line
     -D: DEBUG MODE (DEFAULT = off)
 Output: lines from the matchtable for the clusters specified plus a summary line

 Authors: Granger Sutton, Ph.D.
  Date: May 15, 2019; last revised   May 15, 2019
  Input: The input files are primarily files generated by a panoct run and must be consistent with and across that run.
$prog requires several input files to specify the cluster genome tags,  matchtable, and clusters.

The input files are:

A genome tags file:  has the genome identifier(tag) used in the panoct run,

A matchtable produced by PanOCT with the first column being a cluster identifier and subsequent columns being gene identifiers or a
placeholder indicating no gene for that genome

A clusters file

The output goes to STDOUT and is the matchtable lines for the specified clusters plus a summary line.
_EOB_
    exit;
}

########################################  M A I N  #####################################################
print STDERR "Getting genome namesfrom $genomes_file_name\n";
&get_genomes;
print STDERR "Reading clusters from $clusters_file\n";
&read_clusters_file;
print STDERR "Reading matchtable from $matchtable_file and outputting to stdout\n";
&process_matchtable;
exit(0);
