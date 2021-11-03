#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

##### Parse SARS-CoV-2 genomes from GenBank .gb and .fasta files to generate PanOCT attribute and fasta feature files Script

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
my @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there are BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
my @ordered = ();      # Same as above but sorted by CONTIG then START

# CONSTANTS #
use constant CONTIG => 0;
use constant LOCUS => 1;
use constant START => 2;
use constant STOP => 3;
use constant TYPE => 4;
use constant GENOME => 5;
use constant INVERT => 6;
use constant LENGTH => 7;
use constant DELETE => 8;
use constant ANNOTATION => 9;
# END CONSTANTS #

my $genomes;
my $help;
my $debug;
my $strip_version;
my $att_dir = "";
my $fasta_dir = "";
my $genome_dir = "";

GetOptions('genomes=s' => \ $genomes,
	   'att_dir=s' => \ $att_dir,
	   'fasta_dir=s' => \ $fasta_dir,
	   'genome_dir=s' => \ $genome_dir,
	   'strip_version' => \ $strip_version,
	   'help' => \ $help,
	   'debug' => \ $debug);

if ($att_dir) {
    if (-d $att_dir) {
	if (substr($att_dir, -1, 1) ne "/") {
	    $att_dir .= "/";
	}
	if (substr($att_dir, 0, 1) ne "/") {
	    $att_dir = $cwd . "/$att_dir";
	}
    } else {
	print STDERR "Error with -att_dir $att_dir\n";
	$help = 1;
    }
} else {
    $att_dir = "";
}

if ($fasta_dir) {
    if (-d $fasta_dir) {
	if (substr($fasta_dir, -1, 1) ne "/") {
	    $fasta_dir .= "/";
	}
	if (substr($fasta_dir, 0, 1) ne "/") {
	    $fasta_dir = $cwd . "/$fasta_dir";
	}
    } else {
	print STDERR "Error with -fasta_dir $fasta_dir\n";
	$help = 1;
    }
} else {
    $fasta_dir = "";
}

if ($genome_dir) {
    if (-d $genome_dir) {
	if (substr($genome_dir, -1, 1) ne "/") {
	    $genome_dir .= "/";
	}
	if (substr($genome_dir, 0, 1) ne "/") {
	    $genome_dir = $cwd . "/$genome_dir";
	}
    } else {
	print STDERR "Error with -genome_dir $genome_dir\n";
	$help = 1;
    }
} else {
    $genome_dir = "";
}

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genomes=s' => genomes,
	   'att_dir=s' => att_dir,
	   'fasta_dir=s' => fasta_dir,
	   'genome_dir=s' => genome_dir,
	   'help' => help,
	   'debug' => debug,
	   'strip_version' => strip_version);
_EOB_
    exit(0);
}

#############################################################################################


# subroutine to print contig segments out in fasta format
sub print_fasta { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file_handle, $seq_name, $seq, $beg, $end) = @_;
    #my $full_len = length($seq);
    #print STDERR "$seq_name $beg $end $full_len\n";
    $beg--;
    $end--;
    print $file_handle ">$seq_name\n";
    my $pos;
    my $seq_len = ($end - $beg) + 1;
    for ( $pos = $beg ; $seq_len > 60 ; $pos += 60 ) {
	print $file_handle substr($seq, $pos, 60), "\n";
	$seq_len -= 60;
    }
    print $file_handle substr($seq, $pos, $seq_len), "\n";
    return;
}


# read file which specifies the output file prefix, assembly fasta file, assembly topology file, and anomalies file
open (my $infile, "<", $genomes) || die ("ERROR: cannot open file $genomes\n");
while (my $line = <$infile>)  {
    # clear data structures for the next genome
    @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
    @ordered = ();      # Same as above but sorted by CONTIG then START
    
    chomp $line;
    (my $genome_id, my $gb_file) = split(/\t/, $line);  # split the scalar $line on tab

    my $line;
    my $contig_name = "";

    #read in annotation inforamtion for genome from .gb file
    unless (open (GBFILE, "<", "$gb_file") )  {
	die ("ERROR: can not open .gb file $gb_file.\n");
    }
    my $count = -1;
    my $locus_id = "";
    my $cur_locus_id = "";
    my $type = "";
    my $coord5p;
    my $coord3p;
    my $anno = "";
    my $bp_count;
    my $max_coord = 0;
    my $first = 1;
    while (my $line = <GBFILE>) {
	chomp($line);
	if ($first) {
	    if ($line =~ /^\s*LOCUS\s+(\S+)\s+(\d+)\s+bp/) {
		$locus_id = $1;
		$locus_id =~ s/\.\d+$//; # remove trailing version number if it exists!
		$bp_count = $2;
		$first = 0;
		next;
	    }
	    die ("ERROR: line in .gb file $gb_file is not in expected LOCUS line format:\n$line\n");
	}
 	if ($line =~ /^\s*ORIGIN/) {
	    last;
	}
 	if ($line =~ /^\s*VERSION\s+(\S+)/) {
	    $contig_name = $1;
	    if ($strip_version) {
		$contig_name =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    }
	    next;
	}
	if ($line =~ /^\s*gene\s+(\d+)\.\.(\d+)/) {
	    if ($type ne "") {
		$count++;
		$cur_locus_id = $locus_id . "_" . $count;
		$annotations[$count][CONTIG] = $contig_name;     # contig
		$annotations[$count][ANNOTATION] = $anno;       # annotation
		$annotations[$count][TYPE] = $type;       # type
		$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
		$annotations[$count][GENOME] = $genome_id;     # genome
		$annotations[$count][DELETE] = 0;     # subsumed by another segment
		if ($coord5p <= $coord3p) {
		    $annotations[$count][START] = $coord5p; # start
		    $annotations[$count][STOP] = $coord3p; # stop
		    $annotations[$count][INVERT] = 0;              # invert
		    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
		} else {
		    $annotations[$count][START] = $coord3p; # start
		    $annotations[$count][STOP] = $coord5p; # stop
		    $annotations[$count][INVERT] = 1;              # invert
		    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
		}
	    }
	    $type = "gene";
	    $coord5p = $1;
	    $coord3p = $2;
	    if ($coord5p > $max_coord) {
		$max_coord = $coord5p;
	    }
	    if ($coord3p > $max_coord) {
		$max_coord = $coord3p;
	    }
	    if ($bp_count < $max_coord) {
		die ("ERROR: contig length in LOCUS line ($bp_count) is less than a feature coordinate in .gb file $gb_file:\n$line\n");
	    }
	    next;
	}
	if ($line =~ /^\s*mat_peptide\s+(\d+)\.\.(\d+)/) {
	    if ($type ne "") {
		$count++;
		$cur_locus_id = $locus_id . "_" . $count;
		$annotations[$count][CONTIG] = $contig_name;     # contig
		$annotations[$count][ANNOTATION] = $anno;       # annotation
		$annotations[$count][TYPE] = $type;       # type
		$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
		$annotations[$count][GENOME] = $genome_id;     # genome
		$annotations[$count][DELETE] = 0;     # subsumed by another segment
		if ($coord5p <= $coord3p) {
		    $annotations[$count][START] = $coord5p; # start
		    $annotations[$count][STOP] = $coord3p; # stop
		    $annotations[$count][INVERT] = 0;              # invert
		    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
		} else {
		    $annotations[$count][START] = $coord3p; # start
		    $annotations[$count][STOP] = $coord5p; # stop
		    $annotations[$count][INVERT] = 1;              # invert
		    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
		}
	    }
	    $type = "mat_peptide";
	    $coord5p = $1;
	    $coord3p = $2;
	    if ($coord5p > $max_coord) {
		$max_coord = $coord5p;
	    }
	    if ($coord3p > $max_coord) {
		$max_coord = $coord3p;
	    }
	    if ($bp_count < $max_coord) {
		die ("ERROR: contig length in LOCUS line ($bp_count) is less than a feature coordinate in .gb file $gb_file:\n$line\n");
	    }
	    next;
	}
	if ($line =~ /^\s*stem_loop\s+(\d+)\.\.(\d+)/) {
	    if ($type ne "") {
		$count++;
		$cur_locus_id = $locus_id . "_" . $count;
		$annotations[$count][CONTIG] = $contig_name;     # contig
		$annotations[$count][ANNOTATION] = $anno;       # annotation
		$annotations[$count][TYPE] = $type;       # type
		$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
		$annotations[$count][GENOME] = $genome_id;     # genome
		$annotations[$count][DELETE] = 0;     # subsumed by another segment
		if ($coord5p <= $coord3p) {
		    $annotations[$count][START] = $coord5p; # start
		    $annotations[$count][STOP] = $coord3p; # stop
		    $annotations[$count][INVERT] = 0;              # invert
		    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
		} else {
		    $annotations[$count][START] = $coord3p; # start
		    $annotations[$count][STOP] = $coord5p; # stop
		    $annotations[$count][INVERT] = 1;              # invert
		    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
		}
	    }
	    $type = "stem_loop";
	    $coord5p = $1;
	    $coord3p = $2;
	    if ($coord5p > $max_coord) {
		$max_coord = $coord5p;
	    }
	    if ($coord3p > $max_coord) {
		$max_coord = $coord3p;
	    }
	    if ($bp_count < $max_coord) {
		die ("ERROR: contig length in LOCUS line ($bp_count) is less than a feature coordinate in .gb file $gb_file:\n$line\n");
	    }
	    next;
	}
	if ($line =~ /^\s*5'UTR\s+(\d+)\.\.(\d+)/) {
	    if ($type ne "") {
		$count++;
		$cur_locus_id = $locus_id . "_" . $count;
		$annotations[$count][CONTIG] = $contig_name;     # contig
		$annotations[$count][ANNOTATION] = $anno;       # annotation
		$annotations[$count][TYPE] = $type;       # type
		$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
		$annotations[$count][GENOME] = $genome_id;     # genome
		$annotations[$count][DELETE] = 0;     # subsumed by another segment
		if ($coord5p <= $coord3p) {
		    $annotations[$count][START] = $coord5p; # start
		    $annotations[$count][STOP] = $coord3p; # stop
		    $annotations[$count][INVERT] = 0;              # invert
		    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
		} else {
		    $annotations[$count][START] = $coord3p; # start
		    $annotations[$count][STOP] = $coord5p; # stop
		    $annotations[$count][INVERT] = 1;              # invert
		    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
		}
	    }
	    $type = "5putr";
	    $anno = "5'UTR";
	    $coord5p = $1;
	    $coord3p = $2;
	    if ($coord5p > $max_coord) {
		$max_coord = $coord5p;
	    }
	    if ($coord3p > $max_coord) {
		$max_coord = $coord3p;
	    }
	    if ($bp_count < $max_coord) {
		die ("ERROR: contig length in LOCUS line ($bp_count) is less than a feature coordinate in .gb file $gb_file:\n$line\n");
	    }
	    next;
	}
	if ($line =~ /^\s*3'UTR\s+(\d+)\.\.(\d+)/) {
	    if ($type ne "") {
		$count++;
		$cur_locus_id = $locus_id . "_" . $count;
		$annotations[$count][CONTIG] = $contig_name;     # contig
		$annotations[$count][ANNOTATION] = $anno;       # annotation
		$annotations[$count][TYPE] = $type;       # type
		$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
		$annotations[$count][GENOME] = $genome_id;     # genome
		$annotations[$count][DELETE] = 0;     # subsumed by another segment
		if ($coord5p <= $coord3p) {
		    $annotations[$count][START] = $coord5p; # start
		    $annotations[$count][STOP] = $coord3p; # stop
		    $annotations[$count][INVERT] = 0;              # invert
		    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
		} else {
		    $annotations[$count][START] = $coord3p; # start
		    $annotations[$count][STOP] = $coord5p; # stop
		    $annotations[$count][INVERT] = 1;              # invert
		    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
		}
	    }
	    $type = "3putr";
	    $anno = "3'UTR";
	    $coord5p = $1;
	    $coord3p = $2;
	    if ($coord5p > $max_coord) {
		$max_coord = $coord5p;
	    }
	    if ($coord3p > $max_coord) {
		$max_coord = $coord3p;
	    }
	    if ($bp_count < $max_coord) {
		die ("ERROR: contig length in LOCUS line ($bp_count) is less than a feature coordinate in .gb file $gb_file:\n$line\n");
	    }
	    next;
	}
	if ($line =~ /^\s*\/gene="([^"]*)"/) {
	    $anno = $1;
	    next;
	}
	if ($line =~ /^\s*\/product="([^"]*)"/) {
	    if (($type eq "gene") || ($type eq "mat_peptide")) {
		$anno .= " " . $1;
	    }
	    next;
	}
	if ($line =~ /^\s*\/note="([^"]*)"/) {
	    if ($type eq "stem_loop") {
		$anno .= " " . $1;
	    }
	    next;
	}
    }
    if ($type ne "") {
	$count++;
	$cur_locus_id = $locus_id . "_" . $count;
	$annotations[$count][CONTIG] = $contig_name;     # contig
	$annotations[$count][ANNOTATION] = $anno;       # annotation
	$annotations[$count][TYPE] = $type;       # type
	$annotations[$count][LOCUS] = $cur_locus_id;      # locus_id
	$annotations[$count][GENOME] = $genome_id;     # genome
	$annotations[$count][DELETE] = 0;     # subsumed by another segment
	if ($coord5p <= $coord3p) {
	    $annotations[$count][START] = $coord5p; # start
	    $annotations[$count][STOP] = $coord3p; # stop
	    $annotations[$count][INVERT] = 0;              # invert
	    $annotations[$count][LENGTH] = ($coord3p - $coord5p) + 1;     # length
	} else {
	    $annotations[$count][START] = $coord3p; # start
	    $annotations[$count][STOP] = $coord5p; # stop
	    $annotations[$count][INVERT] = 1;              # invert
	    $annotations[$count][LENGTH] = ($coord5p - $coord3p) + 1;     # length
	}
    }
    my $contig_sequence = "";
    while (my $line = <GBFILE>) {# this is after the ORIGIN line to get the contig sequence
	chomp($line);
	if ($line =~ /^\/\//) {
	    last;
	}
	$line =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	$line =~ tr/a-z/A-Z/; # convert to upper case
	$contig_sequence .= $line
    }
    close (GBFILE);
    my $contig_length = length($contig_sequence);
    if ($contig_length != $bp_count) {
	    die ("ERROR: contig length in LOCUS line ($bp_count) not equal to length of contig sequence after ORIGIN line ($contig_length) in .gb file $gb_file.\n");
    }
    
    if ($debug) {
	print STDERR "DEBUG***annotations\n";
	for (my $j=0; $j < @annotations; $j++) {
	    print STDERR ("$j: $annotations[$j][CONTIG]\t$annotations[$j][LOCUS]\t$annotations[$j][START]\t$annotations[$j][STOP]\t$annotations[$j][TYPE]\t$ordered[$j][ANNOTATION]\t$annotations[$j][GENOME]\t$annotations[$j][INVERT]\n");
	}
    }

    # sort anomalies by contig then by start, store in ordered data-structure. 

    @ordered = sort { $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[TYPE] cmp $b->[TYPE] } @annotations; # sort on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***ordered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][ANNOTATION]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }

    # eliminate overlaps between annotated segments by truncation - delete those that are contained by others - order of precedence: IDENTICAL, CONSERVED, GAPPED, DIVERGED
    for (my $i=0; $i < @ordered; $i++) {
	if ($ordered[$i][DELETE]) {
	    next;
	}
	for (my $j=$i+1; $j < @ordered; $j++) {
	    if ($ordered[$j][DELETE]) {
		last;
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
	    if ($ordered[$i][STOP] < $ordered[$j][STOP]) { #overlap but not contained
		my $overlap = ($ordered[$i][STOP] - $ordered[$j][START]) + 1;
		if (($overlap >= 50) || ($overlap > ($ordered[$i][LENGTH] / 2)) || ($overlap > ($ordered[$j][LENGTH] / 2))) { #too much overlap to ignore
		    if ($ordered[$i][LENGTH] < $ordered[$j][LENGTH]) {
			$ordered[$j][START] = $ordered[$i][STOP] + 1;
			$ordered[$j][LENGTH] -= $overlap;
			if ($ordered[$j][LENGTH] < 25) { #too short to keep
			    $ordered[$j][DELETE] = 1;
			}
		    } else {
			$ordered[$i][STOP] = $ordered[$j][START] - 1;
			$ordered[$i][LENGTH] -= $overlap;
			if ($ordered[$i][LENGTH] < 25) { #too short to keep
			    $ordered[$i][DELETE] = 1;
			    last;
			}
		    }
		}
	    } else { #overlap contained
		if (($ordered[$i][TYPE] eq "gene") && ($ordered[$j][TYPE] eq "mat_peptide")) { #remove gene in favor of mat_peptide
		    $ordered[$i][DELETE] = 1;
		    last;
		} else {
		    $ordered[$j][DELETE] = 1;
		}
	    }
	}
    }

    # resort anomalies by deleted or not, then by contig, then by start, store in ordered data-structure. 

    @ordered = sort { $a->[DELETE] <=> $b->[DELETE] || $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[ANNOTATION] cmp $b->[ANNOTATION] } @ordered; # sort on deleted or not, then on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***filtered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][ANNOTATION]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }
    my $genome_file = $genome_dir . $genome_id . ".fasta";
    my $genome_fp;
    unless (open ($genome_fp, ">", $genome_file) )  {
	die ("cannot open file $genome_file!\n");
    }
    &print_fasta($genome_fp, $contig_name, $contig_sequence, 1, $contig_length); 
    close($genome_fp);
    my $att_file = $att_dir . $genome_id . ".natt";
    my $att_fp;
    unless (open ($att_fp, ">", $att_file) )  {
	die ("cannot open file $att_file!\n");
    }
    my $fasta_file = $fasta_dir . $genome_id . ".nuc";
    my $fasta_fp;
    unless (open ($fasta_fp, ">", $fasta_file) )  {
	die ("cannot open file $fasta_file!\n");
    }
    for (my $j=0; $j < @ordered; $j++) {
	if ($ordered[$j][DELETE]) {
	    last;
	}
	if ((($ordered[$j][STOP] - $ordered[$j][START]) + 1) != $ordered[$j][LENGTH]) {
	    die ("In the ordered array, the stop and start coordinates are not consistent with the length for entry $j ($ordered[$j][STOP],$ordered[$j][START],$ordered[$j][LENGTH])\n");
	}
	my $feature_seq = substr($contig_sequence, ($ordered[$j][START] - 1), (($ordered[$j][STOP] - $ordered[$j][START]) + 1));
	if ($ordered[$j][INVERT]) {
	    my $swap = $ordered[$j][START];
	    $ordered[$j][START] = $ordered[$j][STOP];
	    $ordered[$j][STOP] = $swap;
	    my $tmp_seq = reverse($feature_seq);
	    $feature_seq = $tmp_seq;
	    $feature_seq =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
	}
	&print_fasta($fasta_fp, $ordered[$j][LOCUS], $feature_seq, 1, $ordered[$j][LENGTH]); 
	print $att_fp "$ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][ANNOTATION]\t$ordered[$j][GENOME]\n";
    }
    close($att_fp);
    close ($fasta_fp);
}
close ($infile);
exit(0);
