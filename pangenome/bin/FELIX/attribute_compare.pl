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

my @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
my @ordered = (); # Same as above but sorted by CONTIG then START

# CONSTANTS #
use constant OLD => 0;
use constant NEW => 1;
use constant CONTIG => 0;
use constant LOCUS => 1;
use constant START => 2;
use constant STOP => 3;
use constant ANNOTATION => 4;
use constant GENOME => 5;
use constant INVERT => 6;
use constant BEST => 7;
use constant VALUE => 8;
use constant TYPE => 9;
# END CONSTANTS #

GetOptions('old=s' => \my $old,
	'new=s' => \my $new,
	'help' => \my $help,
	'debug' => \my $debug,
	'percent=f' => \my $percent,
	'name=s' => \my $name);

my $fraction;	
if (!$percent) {
    $fraction = 0.5;
} else {
    $fraction = $percent / 100;
}
# First, read in both attribute files
open(OLD_FILE, "<", $old) || die ("Couldn't open $old\n");
my $count = 0;
while (my $line = <OLD_FILE>) {
    chomp($line);
    my @split_line = split(/\t/,$line);
    $annotations[$count][CONTIG] = $split_line[0];     # contig
    $annotations[$count][LOCUS] = $split_line[1];     # locus_id
    if ($split_line[2] < $split_line[3]) {
	$annotations[$count][START] = $split_line[2]; # start
	$annotations[$count][STOP] = $split_line[3]; # stop
	$annotations[$count][INVERT] = 0;              # invert
    } else {
	$annotations[$count][START] = $split_line[3]; # start
	$annotations[$count][STOP] = $split_line[2]; # stop
	$annotations[$count][INVERT] = 1;              # invert
    }
    $annotations[$count][ANNOTATION] = $split_line[4];     # annotation
    $annotations[$count][GENOME] = $split_line[5];     # genome
    $annotations[$count][BEST] = -1;  # initialize BEST match to impossible array index
    $annotations[$count][VALUE] = 0;  # initialize VALUE to BEST match overlap
    $annotations[$count][TYPE] = OLD;  # initialize TYPE to OLD
    $count++;
}
my $old_count = $count;
close(OLD_FILE);
open(NEW_FILE, "<", $new) || die ("Couldn't open $new\n");
while (my $line = <NEW_FILE>) {
    chomp($line);
    my @split_line = split(/\t/,$line);
    $annotations[$count][CONTIG] = $split_line[0];     # contig
    $annotations[$count][LOCUS] = $split_line[1];     # locus_id
    if ($split_line[2] < $split_line[3]) {
	$annotations[$count][START] = $split_line[2]; # start
	$annotations[$count][STOP] = $split_line[3]; # stop
	$annotations[$count][INVERT] = 0;              # invert
    } else {
	$annotations[$count][START] = $split_line[3]; # start
	$annotations[$count][STOP] = $split_line[2]; # stop
	$annotations[$count][INVERT] = 1;              # invert
    }
    $annotations[$count][ANNOTATION] = $split_line[4];     # annotation
    $annotations[$count][GENOME] = $split_line[5];     # genome
    $annotations[$count][BEST] = -1;  # initialize BEST match to impossible array index
    $annotations[$count][VALUE] = 0;  # initialize VALUE to BEST match overlap
    $annotations[$count][TYPE] = NEW;  # initialize TYPE to NEW
    $count++;
}
my $new_count = $count - $old_count;
my $total_count = $count;
close(NEW_FILE);

if ($debug) {
    print "DEBUG***annotations\n";
    for (my $j=0; $j < @annotations; $j++) {
	print ("$j: $annotations[$j][CONTIG]\t$annotations[$j][LOCUS]\t$annotations[$j][START]\t$annotations[$j][STOP]\t$annotations[$j][ANNOTATION]\t$annotations[$j][GENOME]\t$annotations[$j][INVERT]\t$annotations[$j][BEST]\t$annotations[$j][VALUE]\t$annotations[$j][TYPE]\n");
    }
}

# Second, sort attribute files by contig then by start, store in ordered data-structure. 

@ordered = sort { $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] } @annotations; # sort on contig, then on start

for (my $i=0; $i < @ordered; $i++) {
    for (my $j=$i+1; $j < @ordered; $j++) {
	if ($ordered[$i][TYPE] eq $ordered[$j][TYPE]) {
	    next;
	}
	if ($ordered[$i][CONTIG] lt $ordered[$j][CONTIG]) {
	    last;
	}
	if ($ordered[$i][CONTIG] gt $ordered[$j][CONTIG]) {
	    die ("ERROR: Annotations are not properly sorted!\n");
	}
	if ($ordered[$i][STOP] < $ordered[$j][START]) {
	    last;
	}
	my $overlap;
	if ($ordered[$i][STOP] <= $ordered[$j][STOP]) {
	    $overlap = ($ordered[$i][STOP] - $ordered[$j][START]) + 1;
	} else {
	    $overlap = ($ordered[$j][STOP] - $ordered[$j][START]) + 1;
	}
	my $value = (($overlap / (($ordered[$j][STOP] - $ordered[$j][START]) + 1)) + ($overlap / (($ordered[$i][STOP] - $ordered[$i][START]) + 1))) / 2;
	if ($value >= $fraction) {
	    if ($value > $ordered[$i][VALUE]) {
		$ordered[$i][VALUE] = $value;
		$ordered[$i][BEST] = $j;
	    }
	    if ($value > $ordered[$j][VALUE]) {
		$ordered[$j][VALUE] = $value;
		$ordered[$j][BEST] = $i;
	    }
	}
    }
}

if ($debug) {
    print "DEBUG***ordered\n";
    for (my $j=0; $j < @ordered; $j++) {
	print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][ANNOTATION]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\t$ordered[$j][BEST]\t$ordered[$j][VALUE]\t$ordered[$j][TYPE]\n");
    }
}

my $invert = 0;
my $both = 0;
my $both_invert = 0;
my $added = 0;
my $removed = 0;
my $start = 0;
my $start_invert = 0;
my $stop = 0;
my $stop_invert = 0;
my $same = 0;

open(CHANGED, ">", "$name.CHANGED") || die ("Couldn't open $name.CHANGED\n");
open(SAME, ">", "$name.SAME") || die ("Couldn't open $name.SAME\n");

for (my $i=0; $i < @ordered; $i++) {
    my $j = $ordered[$i][BEST];
    if (($j >= 0) && ($ordered[$j][BEST] == $i)) {
	if ($j < $i) { # we've already processed this one
	    next;
	}
	if ($ordered[$i][TYPE] eq NEW) { # always want the OLD annotation first
	    my $tmp_swap = $j;
	    $j = $i;
	    $i = $tmp_swap;
	}
	if (($ordered[$i][START] == $ordered[$j][START]) && $ordered[$i][STOP] == $ordered[$j][STOP]) {           # if coordinates are an exact match...
	    if ($ordered[$i][INVERT] == $ordered[$j][INVERT]) {                                                   # and if there is no invert...
		print SAME "SAME: $ordered[$i][LOCUS], $ordered[$j][LOCUS]\n";
		$same++;
	    } else {
		if ($ordered[$i][INVERT]) {
		    print CHANGED "INVERT: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t";
		} else {
		    print CHANGED "INVERT: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t";
		}
		if ($ordered[$j][INVERT]) {
		    print CHANGED "INVERT: $ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		} else {
		    print CHANGED "INVERT: $ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		}
		$invert++;
	    }
	} elsif ($ordered[$i][START] != $ordered[$j][START]) {                                                     # if start is off...
	    if ($ordered[$i][STOP] != $ordered[$j][STOP]) {                                                        # and if end is off...
		if ($ordered[$i][INVERT] != $ordered[$j][INVERT]) {
		    if ($ordered[$i][INVERT]) {
			print CHANGED "BOTH + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		    } else {
			print CHANGED "BOTH + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		    }
		    $both_invert++;
		} else {
		    if ($ordered[$i][INVERT]) {
			print CHANGED "BOTH: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n"; 
		    } else {
			print CHANGED "BOTH: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n"; 
		    }
		    $both++;
		}
	    } else {
		if ($ordered[$i][INVERT] != $ordered[$j][INVERT]) {
		    if ($ordered[$i][INVERT]) {
			print CHANGED "DIFF_5P + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		    } else {
			print CHANGED "DIFF_5P + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		    }
		    $start_invert++;
		} else {
		    if ($ordered[$i][INVERT]) {
			print CHANGED "DIFF_5P: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		    } else {
			print CHANGED "DIFF_5P: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		    }
		    $start++;
		}
	    }
	} else {
	    if ($ordered[$i][INVERT] != $ordered[$j][INVERT]) {
		if ($ordered[$i][INVERT]) {
		    print CHANGED "DIFF_3P + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		} else {
		    print CHANGED "DIFF_3P + INVERT: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		}
		$stop_invert++;
	    } else {
		if ($ordered[$i][INVERT]) {
		    print CHANGED "DIFF_3P: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\t$ordered[$j][LOCUS]: ($ordered[$j][STOP],$ordered[$j][START])\n";
		} else {
		    print CHANGED "DIFF_3P: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\t$ordered[$j][LOCUS]: ($ordered[$j][START],$ordered[$j][STOP])\n";
		}
		$stop++;
	    }
	}
    } else {
	# IF AT THIS POINT, NO EXACT MATCH. NOW WE CHECK IF THERE IS AN OMISSION                  # if new has no match something was added		
	if ($ordered[$i][TYPE] eq NEW) { # always want the OLD annotation first
	    if ($ordered[$i][INVERT]) {
		print CHANGED "NEW: N/A\t$ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\n";
	    } else {
		print CHANGED "NEW: N/A\t$ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\n";
	    }
	    $added++;
	} else {                                                                                  # if old has no match something was removed
	    if ($ordered[$i][INVERT]) {
		print CHANGED "REMOVED: $ordered[$i][LOCUS]: ($ordered[$i][STOP],$ordered[$i][START])\tN/A\n";
	    } else {
		print CHANGED "REMOVED: $ordered[$i][LOCUS]: ($ordered[$i][START],$ordered[$i][STOP])\tN/A\n";
	    }
	    $removed++;
	}
    }
}

# PRINT THE STATISTICS HERE
print "Original Feature Count: $old_count\n";
print "Re-annoated Feature Count: $new_count\n";
print "Unchanged Co-ords, Inverted: $invert\n";
print "Both Co-ords Changed: $both\n";
print "Both Co-ords Changed, Inverted: $both_invert\n";
print "Added: $added\n";
print "Removed: $removed\n";
print "5 Prime Changed: $start\n";
print "5 Prime Changed, Inverted: $start_invert\n";
print "3 Prime Changed: $stop\n";
print "3 Prime Changed, Inverted: $stop_invert\n";
print "Unchanged: $same\n";
exit(0);




