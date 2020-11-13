#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in a 3 column tab delimitted file from standard in. The 3 columns are the first three from mash dist runs: name1, name2, distance. This program outputs to standard out the names of genomes that are redundant with other genomes.

Input Flags:
-threshold - used for column 3 to see if two genomes are redundant <= (default 0.0002)
-help - outputs this help text";

GetOptions('threshold=s' => \my $threshold,
	   'help' => \my $help);
	
if(!$threshold) {
    $threshold = 0.0002;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my %removed = ();         # key = name, value doesn't matter - existence shows that it has been removed as redundant so cannot be used to remove others.

#####################################################################################################
sub read_input    # Read in the 3 column input (name1, name2, distance)
{
    my $line;
    while ($line = <STDIN>) {
	(my $name1, my $name2, my $distance) = split(/\t/, $line, 3);
	if (($name1 ne $name2) && (!defined($removed{$name1}) && !defined($removed{$name2})) && ($distance <= $threshold)) {
	    print STDOUT "$name2\n";
	    $removed{$name2} = 1;
	}
    }
    return;
 }

{#main
    &read_input;
}
