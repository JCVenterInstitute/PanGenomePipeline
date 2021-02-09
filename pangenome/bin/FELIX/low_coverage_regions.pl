#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in a 3 column tab delimitted file from the coverage option file. The 3 columns are the contig name, the contig position, and the coverage depth - missing contig positions are assumed to have coverage 0. This program outputs to standard output a 3 column tab delimitted file: column 1 is contig name, column 2 is start of 0 coverage region, and column 3 is stop of 0 coverage region - final row for a contig will have a stop of end representing the end coordinate of the contig which is unknown to the program.

Input Flags:
-threshold - used for column 3 to see if two genomes should be clustered <= (default 0.02)
-coverage - file name of coverage data
-help - outputs this help text";

GetOptions('threshold=s' => \my $threshold,
	   'coverage=s' => \my $coverage_file,
	   'help' => \my $help);
	
if(!$threshold) {
    $threshold = 0;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my $last_high = 0; # position of last high coverage base pair
my $current_contig = ""; # current contig name

#####################################################################################################
my $coverage_handle;
unless (open ($coverage_handle, "<", $coverage_file) )  {
    die ("cannot open file $coverage_file!\n");
}
my $line;
while ($line = <$coverage_handle>) {
    (my $contig, my $position, my $coverage) = split(/\t/, $line, 3);
    if ($contig ne $current_contig) {
	if ($current_contig ne "") {
	    $last_high++;
	    print "$current_contig\t$last_high\tend\n";
	}
	$last_high = 0;
	$current_contig = $contig;
    }
    if ($coverage > $threshold) {
	if ($position > ($last_high + 1)) {
	    $last_high++;
	    my $stop = $position - 1;
	    print "$current_contig\t$last_high\t$stop\n";
	}
	$last_high = $position;
    }
}
if ($current_contig ne "") {
    $last_high++;
    print "$current_contig\t$last_high\tend\n";
}
close($coverage_handle);
