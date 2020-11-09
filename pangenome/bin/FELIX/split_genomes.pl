#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in a multifasta file of complete genomes where the header identifies which genome a contig belongs to as standard input. The genomes are output as separate multifasta files in the current directory with names based on the fasta header.

Input Flags:
-name_fields - comma separates list of the field indices from the fasta header (separated by |) to use for genome names (starts at 0) (required)
-segment_field - index of field from fasta header that is the segment identifier
-help - Outputs this help text";

GetOptions('name_fields=s' => \my $name_fields,
	   'segment_field=s' => \my $segment_field,
	   'help' => \my $help);
	
if(!$name_fields || !$segment_field){
    print("$help_text\n");
    exit;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my %segments_seen = ();         # key = segment name based on fasta header, value doesn't matter - existence shows that it has been seen already and should be ignored.
(my @header_indices) = split(/,/, $name_fields);
my $first_index = shift (@header_indices);

#####################################################################################################
sub read_genomes    # Read in the multifasta file, determine if a segment is a duplicate, then write the segment to a file if not.
{
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    while ($line = <STDIN>) {
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	$title =~ s/^>//; # remove leading > from first sequence
	$title =~ s/^\s+//; # remove leading spaces
	my @fields = split(/\|/, $title);  # split the fasta header on | which separates the expected fields
	my $name = $fields[$first_index];
	foreach my $index (@header_indices) {
	    $name .= "_" . $fields[$index];
	}
	$name =~ tr{/}{:};
	my $segment_key = $name . "_" . $fields[$segment_field];
	if (!defined $segments_seen{$segment_key}) {
	    $segments_seen{$segment_key} = 1;
	    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	    &print_fasta($name, $title, $sequence);
	}
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    return;
 }

#####################################################################################################
# subroutine to print contigs out in fasta format
sub print_fasta { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file_name, $seq_name, $seq) = @_;
    my $file_handle;
    unless (open ($file_handle, ">>", $file_name) )  {
	die ("cannot open file $file_name!\n");
    }
    print $file_handle ">$seq_name\n";
    my $pos;
    my $seq_len = length($seq);
    for ( $pos = 0 ; $seq_len > 60 ; $pos += 60 ) {
	print $file_handle substr($seq, $pos, 60), "\n";
	$seq_len -= 60;
    }
    print $file_handle substr($seq, $pos, $seq_len), "\n";
    return;
}

{#main
    &read_genomes;
}
