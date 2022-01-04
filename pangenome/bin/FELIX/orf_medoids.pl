#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $help_text = "This program reads in a fasta file with cluster medoids and outputs to standard out a tab delimited file with which reading frames (1-6) have no stop codons other than a terminating stop codon (1 - no stop codons, 0 - one or more stop codons).

Input Flags:
-medoids - A fasta file of cluster medoids
-help - Outputs this help text";

GetOptions('medoids=s' => \my $medoids,
	'help' => \my $help);
	
if($help){
    print("$help_text\n");
    exit;
}
	
if(!$medoids){
    die("Error: One or more of the required file arguments are missing\n$help_text\n");
}

#Globals


###################################################################################################
sub read_medoids {  # read in the medoids for the PGG

    my $medoidfile;
    unless (open ($medoidfile, "<", $medoids) )  {
	die ("cannot open medoid fasta file: $medoids!\n");
    }
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    my @frames;
    while ($line = <$medoidfile>) {
	(my $title, my $tmp_seq) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	my $sequence = uc($tmp_seq); # convert to uppercase
	$sequence =~ s/[^A-Z]//g; # remove any non-alphabet characters
	my $revcomp = reverse($sequence); # reverse
	$revcomp =~ tr/AGCTYRWSKMDVHB/TCGARYWSMKHBDV/; # complement
	$frames[0] = 1;
	$frames[1] = 1;
	$frames[2] = 1;
	$frames[3] = 1;
	$frames[4] = 1;
	$frames[5] = 1;
	my $last_index = length($sequence) - 4; # don't want a terminating stop codon
	foreach my $frame (0, 1, 2) {
	    my $index = $frame;
	    while ($index <= $last_index) {
		my $codon = substr($sequence, $index, 3);
		$index += 3;
		if (($codon eq "TAA") || ($codon eq "TAG") || ($codon eq "TGA")) {
		    $frames[$frame] = 0;
		    last;
		}
	    }
	}
	foreach my $frame (0, 1, 2) {
	    my $index = $frame;
	    while ($index <= $last_index) {
		my $codon = substr($revcomp, $index, 3);
		$index += 3;
		if (($codon eq "TAA") || ($codon eq "TAG") || ($codon eq "TGA")) {
		    $frames[$frame + 3] = 0;
		    last;
		}
	    }
	}
	my $outline = join("\t", @frames);
	print "$outline\n";
	$title = ""; # clear the title for the next medoid
	$tmp_seq = ""; #clear out the sequence for the next medoid
    }
    $/ = $save_input_separator; # restore the input separator
    close ($medoidfile);
    return;
}



{#main
    &read_medoids;
}
