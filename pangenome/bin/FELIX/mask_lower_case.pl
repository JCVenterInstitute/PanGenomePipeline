#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

##### Attribute Compare Script

use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;
use List::Util qw[min max];

my $cwd = getcwd;
my $commandline = join (" ", @ARGV);
print STDERR "$commandline\n";
my %anchor_coords = (); #hash of arrays containing anchor coordinates for each contig, key = contig_id

my $genomes;
my $help;
my $debug;
my $strip_version;
my $out_file;
my $combine_ids = 0;

GetOptions('genomes=s' => \ $genomes,
	   'out_file=s' => \ $out_file,
	   'combine_ids' => \ $combine_ids,
	   'strip_version' => \ $strip_version,
	   'help' => \ $help,
	   'debug' => \ $debug);


if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genomes=s' => \ genomes,
	   'out_file=s' => \ out_file,
	   'combine_ids' => \ combine_ids,
	   'strip_version' => \ strip_version,
	   'help' => \ help,
	   'debug' => \ debug);
_EOB_
    exit(0);
}

# subroutine to print contig segments out in fasta format
sub print_fasta { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file_handle, $seq_name, $seq, $seq_len) = @_;
    print $file_handle ">$seq_name\n";
    my $pos;
    for ( $pos = 0 ; $seq_len > 60 ; $pos += 60 ) {
	print $file_handle substr($seq, $pos, 60), "\n";
	$seq_len -= 60;
    }
    print $file_handle substr($seq, $pos, $seq_len), "\n";
    return;
}

my $output_handle;
open($output_handle, ">", $out_file) || die ("Couldn't open $out_file\n");

# read file which specifies the output file prefix, assembly fasta file, and anchors coordinates file
open (my $infile, "<", $genomes) || die ("ERROR: cannot open file $genomes\n");
while (my $genome_info = <$infile>)  {
    # clear data structures for the next genome
    %anchor_coords = (); #hash of arrays containing anchor coordinates for each contig, key = contig_id
    
    chomp $genome_info;
    (my $output, my $genome, my $anchor_file) = split(/\t/, $genome_info);  # split the scalar $line on tab

    # read in anchors file
    my $anchor_handle;
    open($anchor_handle, "<", $anchor_file) || die ("Couldn't open $anchor_file\n");
    while (my $line = <$anchor_handle>) {
	chomp $line;
	(my $contig_id, my $start, my $stop) = split(/\t/, $line);  # split the scalar $line on tab
	if ($start > $stop) {
	    my $tmp = $stop;
	    $stop = $start;
	    $start = $tmp;
	}
	if ($strip_version) {
	    $contig_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	if (!defined $anchor_coords{$contig_id}) {
	    $anchor_coords{$contig_id} = [];
	}
	push (@{ $anchor_coords{$contig_id} }, ($start . ":" . $stop));
    }
    close($anchor_handle);
     
    # read in contigs from genome file
    my $contigfile;
    unless (open ($contigfile, "<", $genome) )  {
	die ("cannot open genome file: $genome!\n");
    }
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    while ($line = <$contigfile>) {
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$title =~ s/^>//; # remove leading >
	$id =~ s/^>//; # remove leading >
	if ($strip_version) {
	    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	$sequence =~ tr/A-Z/a-z/; # everything lower case
	my $contig_length = length($sequence);
	#print STDERR "$id\t$contig_length\n";
	foreach my $coords (@{ $anchor_coords{$id} }) {
	    (my $start, my $stop) = split(/:/, $coords);  # split on :
	    $start--;
	    my $len = $stop - $start;
	    my $anchor = substr($sequence, $start, $len);
	    $anchor =~ tr/a-z/A-Z/; # upper case the anchor
	    substr($sequence, $start, $len) = $anchor;
	}
	if ($combine_ids) {
	    $title = $output . "_" . $title;
	}
	&print_fasta($output_handle, $title, $sequence, $contig_length);
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($contigfile);
}
close ($output_handle);

exit(0);
