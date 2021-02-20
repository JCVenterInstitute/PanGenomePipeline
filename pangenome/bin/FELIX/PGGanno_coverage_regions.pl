#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in 3 files: a read depth coverage file (column 1 contig name, column 2 contig size, column 3 depth of read coverage), a PGG attributes file, and a contig fasta file. Two output files must be specified: main for those contigs with good PGG annotation coverage, and left for the other contigs.

Input Flags:
-threshold - given as a percentage (0-100) for how much of the contig must be spanned by PGG attributes to be considered good
-read_coverage - file with contig size and read coverage
-PGG_coverage - file with PGG attributes
-contigs - contig fasta file
-main - output file for good contigs
-left - output file for the rest of the contigs
-help - outputs this help text";

GetOptions('threshold=s' => \my $threshold,
	   'read_coverage=s' => \my $read_coverage,
	   'PGG_coverage=s' => \my $PGG_coverage,
	   'contigs=s' => \my $contigs,
	   'main=s' => \my $main,
	   'left=s' => \my $left,
	   'help' => \my $help);
	
if(!$threshold) {
    $threshold = 50;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my %contig_size = (); # size of contigs
my %contig_rdepth = (); # read depth of contigs
my %contig_perc = (); # percentage of contig covered by PGG attributes
my $start; # start coordinate of PGG attribute
my $stop; # stop coordinate of PGG attribute
my $prev_stop = 0; # stop coordinate of previous PGG attribute
my $cur_contig; # current contig name
my $prev_contig = ""; # previous contig name
my $PGG_anno_bp = 0; # number of contig base pairs covered by PGG attributes
my $rcfile; # file descriptor for read coverage file
my $PGGfile; # file descriptor for PGG coverage file
my $contigfile; # file descriptor for contig fasta file
my $mainfile; # file descriptor for main output file
my $leftfile; # file descriptor for left output file
my $line; # buffer for reading lines of a file

#####################################################################################################
unless (open ($rcfile, "<", $read_coverage) )  {
    die ("cannot open read coverage file: $read_coverage!\n");
}
while ($line = <$rcfile>) {
    (my $contig_name, my $size, my $read_depth) = split(/\t/, $line, 3);
    $contig_size{$contig_name} = $size;
    $contig_rdepth{$contig_name} = $read_depth;
}
close($rcfile);

unless (open ($PGGfile, "<", $PGG_coverage) )  {
    die ("cannot open PGG coverage file: $PGG_coverage!\n");
}
while ($line = <$PGGfile>) {
    my $cluster;
    ($cur_contig, $cluster, $start, $stop) = split(/\t/, $line, 4);
	if (!defined $contig_size{$cur_contig}) {
	    die ("ERROR $cur_contig in PGG coverage file $PGG_coverage but not in read coverage file $read_coverage!\n");
	}
    if ($prev_contig eq "") {
	$prev_contig = $cur_contig;
    } elsif ($cur_contig ne $prev_contig) {
	$contig_perc{$prev_contig} = 100 * ($PGG_anno_bp / $contig_size{$prev_contig});
	$PGG_anno_bp = 0;
	$prev_stop = 0;
	$prev_contig = $cur_contig;
    }
    if ($start > $stop) {
	my $tmp = $start;
	$start = $stop;
	$stop = $tmp;
    }
    if ($start <= $prev_stop) {
	$start = $prev_stop + 1;
    }
    $PGG_anno_bp += ($stop - $start) + 1;
    $prev_stop = $stop;
}
close($PGGfile);

unless (open ($contigfile, "<", $contigs) )  {
    die ("cannot open input contig fasta file: $contigs!\n");
}
unless (open ($mainfile, ">", $main) )  {
    die ("cannot open main output fasta file: $main!\n");
}
unless (open ($leftfile, ">", $left) )  {
    die ("cannot open left output fasta file: $left!\n");
}
my $is_main = 0; # flag to indicate current contig is main or left
while ($line = <$contigfile>) {
    if ($line =~ /^>/) {
	chomp $line;
	my @fields = split(/\s+/, $line);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$id =~ s/>\s*//; # remove leading > and spaces
	if (!defined $contig_size{$id}) {
	    die ("ERROR $id in contig fasta file $contigs but not in read coverage file $read_coverage!\n");
	}
	if ((defined $contig_size{$id}) && ($contig_perc{$id} > $threshold)) {
	    $is_main = 1;
	    print $mainfile "$line size_$contig_size{$id} rdepth_$contig_rdepth{$id} PGGperc_$contig_perc{$id}\n";
	} else {
	    $is_main = 0;
	    print $leftfile "$line size_$contig_size{$id} rdepth_$contig_rdepth{$id} PGGperc_$contig_perc{$id}\n";
	}
    } else {
	if ($is_main) {
	    print $mainfile $line;
	} else {
	    print $leftfile $line;
	}
    }
}
close($contigfile);
close($mainfile);
close($leftfile);
