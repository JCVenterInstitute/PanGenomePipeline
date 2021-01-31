#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

##### Attribute Compare Script

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;
use List::Util qw[min max];

my @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there are STRAND, USED, and SPECIES fields 2) coordinates are now smallest then largest, not start then stop 3) there is an STRAND field to indicate strand rather than STOP being smaller than START
my @sorted_annotations = (); # Same as above but sorted by CONTIG then START
my %species_order = (); # Key is species name, value is order for output
my $blank = "\t\t\t\t\t\t";
my $highest;

# CONSTANTS #
use constant CONTIG => 0;
use constant LOCUS => 1;
use constant START => 2;
use constant STOP => 3;
use constant ANNOTATION => 4;
use constant GENOME => 5;
use constant PID => 6;
use constant STRAND => 7;
use constant USED => 8;
use constant SPECIES => 9;
use constant BEST => 10;
use constant VALUE => 11;
use constant ORDER => 12;
# END CONSTANTS #

GetOptions('att_files=s' => \my $att_files,
	'help' => \my $help,
	'debug' => \my $debug,
	'percent=f' => \my $percent,
	'name=s' => \my $name);
if ($help || !$att_files) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('att_files=s' => \ att_files,
	'help' => \ help,
	'debug' => \ debug,
	'percent=f' => percent,
	'name=s' => \ name);
_EOB_
    exit(0);
}
my $fraction;	
if (!$percent) {
    $fraction = 0.5;
} else {
    $fraction = $percent / 100;
}
# First, read in attribute files
# read file which specifies the species name, and attribute file
my $index = 0;
my $count = 0;
open (my $infile, "<", $att_files) || die ("ERROR: cannot open list of attribute files $att_files\n");
while (my $list_line = <$infile>)  {
    chomp $list_line;
    (my $species, my $species_att_file) = split(/\t/, $list_line);  # split the scalar $line on tab
    if (defined $species_order{$species}) {
	die ("ERROR: $species occurs more than once as species name in $att_files\n");
    }
    my $order = $index;
    $highest = $order;
    $species_order{$species} = $index++;

    # read in attribute file
    open(my $species_att_fd, "<", $species_att_file) || die ("Couldn't open $species_att_file\n");
    while (my $line = <$species_att_fd>) {
	chomp($line);
	my @split_line = split(/\t/,$line);
	$annotations[$count][CONTIG] = $split_line[0];     # contig
	$annotations[$count][LOCUS] = $split_line[1];     # locus_id
	if ($split_line[2] < $split_line[3]) {
	    $annotations[$count][START] = $split_line[2]; # start
	    $annotations[$count][STOP] = $split_line[3]; # stop
	    $annotations[$count][STRAND] = "+";              # invert
	} else {
	    $annotations[$count][START] = $split_line[3]; # start
	    $annotations[$count][STOP] = $split_line[2]; # stop
	    $annotations[$count][STRAND] = "-";              # invert
	}
	$annotations[$count][ANNOTATION] = $split_line[4];     # annotation
	$annotations[$count][GENOME] = $split_line[5];     # genome
	$annotations[$count][PID] = $split_line[6];     # genome
	$annotations[$count][USED] = 0;  # initialize to indicate this attribute line has not been used yet
	$annotations[$count][SPECIES] = $species;  # species
	$annotations[$count][ORDER] = $order;  # order of species
	$annotations[$count][BEST] = -1;  # initialize BEST match to impossible array index
	$annotations[$count][VALUE] = 0;  # initialize VALUE to BEST match overlap
	$count++;
    }
    close($species_att_fd);
}
close($infile);

if ($debug) {
    print "DEBUG***annotations\n";
    for (my $j=0; $j < @annotations; $j++) {
	print ("$j: $annotations[$j][CONTIG]\t$annotations[$j][LOCUS]\t$annotations[$j][START]\t$annotations[$j][STOP]\t$annotations[$j][STRAND]\t$annotations[$j][ANNOTATION]\t$annotations[$j][GENOME]\t$annotations[$j][PID]\t$annotations[$j][USED]\t$annotations[$j][SPECIES]\t$annotations[$j][ORDER]\t$annotations[$j][BEST]\t$annotations[$j][VALUE]\n");
    }
}

# Second, sort attribute files by contig then by start, store in ordered data-structure. 

@sorted_annotations = sort { $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[ORDER] <=> $b->[ORDER] } @annotations; # sort on contig, then on start, then on species order

for (my $i=0; $i < @sorted_annotations; $i++) {
    if ($sorted_annotations[$i][USED]) {
	next;
    }
    print "$sorted_annotations[$i][CONTIG]";
    for (my $order = 0; $order < $sorted_annotations[$i][ORDER]; $order++) {
	print $blank;
    }
    print "\t$sorted_annotations[$i][LOCUS]\t$sorted_annotations[$i][START]\t$sorted_annotations[$i][STOP]\t$sorted_annotations[$i][STRAND]\t$sorted_annotations[$i][ANNOTATION]\t$sorted_annotations[$i][PID]";
    $sorted_annotations[$i][USED] = 1;
    my $best_pid = $sorted_annotations[$i][PID];
    my $best_species = $sorted_annotations[$i][SPECIES];
    for (my $order = $sorted_annotations[$i][ORDER] + 1; $order <= $highest; $order++) {
	for (my $j=$i+1; $j < @sorted_annotations; $j++) {
	    if ($sorted_annotations[$i][STOP] < $sorted_annotations[$j][START]) {
		last;
	    }
	    if ($sorted_annotations[$i][CONTIG] lt $sorted_annotations[$j][CONTIG]) {
		last;
	    }
	    if ($sorted_annotations[$j][USED]) {
		next;
	    }
	    if ($sorted_annotations[$j][ORDER] != $order) {
		next;
	    }
	    my $overlap;
	    if ($sorted_annotations[$i][STOP] <= $sorted_annotations[$j][STOP]) {
		$overlap = ($sorted_annotations[$i][STOP] - $sorted_annotations[$j][START]) + 1;
	    } else {
		$overlap = ($sorted_annotations[$j][STOP] - $sorted_annotations[$j][START]) + 1;
	    }
	    my $value = (($overlap / (($sorted_annotations[$j][STOP] - $sorted_annotations[$j][START]) + 1)) + ($overlap / (($sorted_annotations[$i][STOP] - $sorted_annotations[$i][START]) + 1))) / 2;
	    if ($value >= $fraction) {
		if ($value > $sorted_annotations[$i][VALUE]) {
		    $sorted_annotations[$i][VALUE] = $value;
		    $sorted_annotations[$i][BEST] = $j;
		}
		if ($value > $sorted_annotations[$j][VALUE]) {
		    $sorted_annotations[$j][VALUE] = $value;
		    $sorted_annotations[$j][BEST] = $i;
		}
	    }
	}
	if ($sorted_annotations[$i][BEST] >= 0) {
	    my $j = $sorted_annotations[$i][BEST];
	    if ($j <= $i) {
		print STDERR "WARNING: unexpected best overlap to previous attribute - probably a split gene: $j:$i\n$sorted_annotations[$j][CONTIG]\t$sorted_annotations[$j][LOCUS]\t$sorted_annotations[$j][START]\t$sorted_annotations[$j][STOP]\t$sorted_annotations[$j][STRAND]\t$sorted_annotations[$j][ANNOTATION]\t$sorted_annotations[$j][PID]\n$sorted_annotations[$i][CONTIG]\t$sorted_annotations[$i][LOCUS]\t$sorted_annotations[$i][START]\t$sorted_annotations[$i][STOP]\t$sorted_annotations[$i][STRAND]\t$sorted_annotations[$i][ANNOTATION]\t$sorted_annotations[$i][PID]\n";
	    }
	    print "\t$sorted_annotations[$j][LOCUS]\t$sorted_annotations[$j][START]\t$sorted_annotations[$j][STOP]\t$sorted_annotations[$j][STRAND]\t$sorted_annotations[$j][ANNOTATION]\t$sorted_annotations[$j][PID]";
	    $sorted_annotations[$j][USED] = 1;
	    if ($sorted_annotations[$j][PID] > $best_pid) {
		$best_pid = $sorted_annotations[$j][PID];
		$best_species = $sorted_annotations[$j][SPECIES];
	    }
	} else {
	    print $blank;
	}
    }
    print "\t$best_species\n";
}

if ($debug) {
    print "DEBUG***sorted_annotations\n";
    for (my $j=0; $j < @sorted_annotations; $j++) {
	print ("$j: $sorted_annotations[$j][CONTIG]\t$sorted_annotations[$j][LOCUS]\t$sorted_annotations[$j][START]\t$sorted_annotations[$j][STOP]\t$sorted_annotations[$j][STRAND]\t$sorted_annotations[$j][ANNOTATION]\t$sorted_annotations[$j][GENOME]\t$sorted_annotations[$j][PID]\t$sorted_annotations[$j][USED]\t$sorted_annotations[$j][SPECIES]\t$sorted_annotations[$j][ORDER]\t$sorted_annotations[$j][BEST]\t$sorted_annotations[$j][VALUE]\n");
    }
}


exit(0);




