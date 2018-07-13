#!/usr/bin/env perl
#Copy (C) 2011-2012  The J. Craig Venter Institute (JCVI).  All rights reserved
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


my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
getopts ('DhM:m:C:c:I:');#M is max, m is min, C is centroids file, c is cluster weights file - all required unless -I is used which is a file of cluster_ids
our ($opt_D, $opt_h, $opt_M, $opt_m, $opt_C, $opt_c, $opt_I);

my $max_cutoff;
my $min_cutoff;
my $centroid_file_name;
my $cluster_file_name;
my $cluster_ids_file_name;
my %cluster_size = ();
my %cluster_ids =();
my $use_cluster_ids = 0;
my $DEBUG = 0;

if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if (($opt_C) && (-s "$opt_C")) {
    $centroid_file_name = $opt_C;
} else {
    print STDERR "Error with -C\n";
    &option_help;
}
if ($opt_I) {
    if (-s "$opt_I") {
	$use_cluster_ids = 1;
	$cluster_ids_file_name = $opt_I;
    } else {
	print STDERR "Error with -I\n";
	&option_help;
    }
} else {
    if (($opt_c) && (-s "$opt_c")) {
	$cluster_file_name = $opt_c;
    } else {
	print STDERR "Error with -c\n";
	&option_help;
    }
    if ($opt_M) {
	if ($DEBUG) {
	    print STDERR "max cutoff $opt_M\n";
	}
	if (!(looks_like_number($opt_M))) {
	    die ("ERROR: $opt_M is not a number for max cutoff (-M)!\n");
	}
	if ($opt_M <= 0) {
	    die ("ERROR: $opt_M is not > 0 for max cutoff (-M)!\n");
	}
	$max_cutoff = $opt_M;
    } else {
	die ("ERROR: max cutoff must be specified using -M");
    }
    if ($opt_m) {
	if ($DEBUG) {
	    print STDERR "min cutoff $opt_m\n";
	}
	if (!(looks_like_number($opt_m))) {
	    die ("ERROR: $opt_m is not a number for min cutoff (-m)!\n");
	}
	if ($opt_m < 0) {
	    die ("ERROR: $opt_m is not >= 0 for min cutoff (-m)!\n");
	}
	$min_cutoff = $opt_m;
    } else {
	die ("ERROR: min cutoff must be specified using -m");
    }
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - Choose centroids takes a multifasta file of centroids (output from PanOCT) and a cluster weights file (output from PanOCT)
         and outputs a multifasta file of centroids whose cluster size is greater than or equal to the minimum cutoff and
         less than or equal to the maximum cutoff.

Copy (C) 2013  The J. Craig Venter Institute (JCVI).  All rights reserved

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

Citation:  

  Usage: $prog <options>
Example: choose_centroids.pl -M 20 -m 2 -C centroids.fasta -c cluster_weights.txt > selected_centroids.fasta
 Option:
     -h: print this help page
     -I: cluster ids file with one cluster id per line this overrides the -M, -m, and -c opptions
     -M: max threshold, cluster size must be less than or equal to this to be selected
     -m: min threshold, cluster size must be greater than or equal to this to be selected
     -C: centroids multifasta file from PanOCT, format is >centroid_7 for the centroid of cluster number 7
     -c: cluster weights file from PanOCT output, tab delimited, column 1 is cluster number, column 2 is cluster size
     -D: DEBUG MODE (DEFAULT = off)
 Output: the selected centroids are output as a multifasta file to stdout
 Author: Granger Sutton, Ph.D.
 Date: 11/12/13
_EOB_
   exit:
}

sub get_cluster_ids {  # obtain list of genomes to compareread cluster ids into a hash
   
    open (my $infile, "<", "$cluster_ids_file_name") || die ("ERROR: can't open file $cluster_ids_file_name\n");
    while (<$infile>)  {
	chomp;
	my $cluster_id = $_;
	$cluster_id =~ s/\r$//; #strip ^M for Windows files
	$cluster_id =~ s/\s+$//; #strip trailing white space
	$cluster_id =~ s/\r$//; #do it again in case
	if (defined $cluster_ids{$cluster_id})  {
               die ("ERROR:  You have more than one occurance of $cluster_id in $cluster_ids_file_name!\n");
	}
	else  {
	    $cluster_ids{$cluster_id} = 1; # used to be genome_hash
	}
    }  
    close($infile);
    return;
}

sub get_cluster_sizes {

    my $number = "";
    my $size = "";

    unless (open (CLUSTERFILE, $cluster_file_name) )  {
	die ("ERROR: can not open file $cluster_file_name.\n");
    }
    while (<CLUSTERFILE>) {
	my @cluster_line = ();
	chomp;
	@cluster_line = split(/\t/, $_);  # split the scalar $cluster_line on tab
	$number = $cluster_line[0];
	if ($number eq "") {
	    die ("cluster number must not be empty/null in cluster file\n");
	}
	if (!(looks_like_number($number))) {
	    die ("ERROR: $number is not a number for cluster number in cluster file!\n");
	}
	if ($number <= 0) {
	    die ("ERROR: $number is not > 0 for cluster number in cluster file!\n");
	}
	$size = $cluster_line[1];
	if (!(looks_like_number($size))) {
	    die ("ERROR: $size is not a number for cluster size in cluster file!\n");
	}
	if ($size <= 0) {
	    die ("ERROR: $size is not > 0 for cluster size in cluster file!\n");
	}
	if ($size eq "") {
	    die ("cluster size must not be empty/null in cluster file\n");
	}
	$cluster_size{$number} = $size;
    }
    close (CLUSTERFILE);
    return;
}

sub select_centroids {

  my @line = ();
  my @centroid_id;
  my $centroid_num;
  my $size;
  my $id;
  my $title = "";
  my $sequence = "";

  unless (open (CENTROIDFILE, $centroid_file_name) )  {
    die ("can't open file $centroid_file_name.\n");
  }
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<CENTROIDFILE>) {
    ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool) also removes the leading > and new line from the header
    @line = split(/\s+/, $title);  # split the scalar $title on space or tab (to separate the identifier from the header and store in array @line
    $id = $line[0]; #centroid_number expected here
    @centroid_id = split('_', $id);  # split the scalar $id at the _ character
    if ($centroid_id[0] ne "centroid") {
	die ("ERROR: centroid fasta header line id: $id - not in expected format of >centroid_\n");
    }
    $centroid_num = $centroid_id[1];
    if ($use_cluster_ids) {
	if (defined $cluster_ids{$centroid_num})  {
	    print STDOUT ">$title\n";
	    print STDOUT $sequence
	}
    } else {
	$size = $cluster_size{$centroid_num};
	if (($size >= $min_cutoff) && ($size <= $max_cutoff)) {
	    print STDOUT ">$title\n";
	    print STDOUT $sequence
	}
    }
    $title = ""; # clear the title for the next round.
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (CENTROIDFILE);
  return;
}

########################################  M A I N  #####################################################
if ($use_cluster_ids) {
    print STDERR "Getting cluster ids from $cluster_ids_file_name\n";
    &get_cluster_ids;
    print STDERR "Selecting centroids from $centroid_file_name\n";
} else {
    print STDERR "Getting cluster sizes from $cluster_file_name\n";
    &get_cluster_sizes;
    print STDERR "Selecting centroids from $centroid_file_name with Max = $max_cutoff and Min = $min_cutoff\n";
}
&select_centroids;
exit(0);
