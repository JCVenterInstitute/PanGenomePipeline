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
getopts ('Dhm:M:p:P:T:d:G:');
our ($opt_D, $opt_h, $opt_m, $opt_M, $opt_p, $opt_P, $opt_T, $opt_d, $opt_G);

my $matchtable_file_name;
my $genomes_file_name;
my $hist_dir;
my @genome_array = ();
my @genome_hists = ();
my @genome_min = ();
my @genome_max = ();
my $genome_number;
my $DEBUG = 0;
my $min_thresh = undef;
my $max_thresh = undef;
my $min_perc;
my $max_perc;

if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if (($opt_T) && (-s "$opt_T")) {
    $matchtable_file_name = $opt_T;
} else {
    print STDERR "Error with -T\n";
    &option_help;
}
if (($opt_G) && (-s "$opt_G")) {
    $genomes_file_name = $opt_G;
} else {
    print STDERR "Error with -G\n";
    &option_help;
}
if ($opt_d) {
    $hist_dir = $opt_d;
} else {
    print STDERR "Error with -d\n";
    &option_help;
}
if ($opt_p) {
    if (!(looks_like_number($opt_p))) {
	die ("ERROR: $opt_p is not a number(-p)!\n");
    }
    if (($opt_p > 100) || ($opt_p < 0)) {
	die ("ERROR: $opt_p is not between 0-100 (-p)!\n");
    }
    $min_perc = $opt_p;
} else {
    $min_perc = 5;#default
}
if ($opt_P) {
    if (!(looks_like_number($opt_P))) {
	die ("ERROR: $opt_P is not a number(-P)!\n");
    }
    if (($opt_P > 100) || ($opt_P < 0)) {
	die ("ERROR: $opt_P is not between 0-100 (-P)!\n");
    }
    $max_perc = $opt_P;
} else {
    $max_perc = 95;#default
}
if ($opt_m) {
    if (!(looks_like_number($opt_m))) {
	die ("ERROR: $opt_m is not a number(-m)!\n");
    }
    if ($opt_m < 0) {
	die ("ERROR: $opt_m is not > 0 (-m)!\n");
    }
    $min_thresh = $opt_m;
}
if ($opt_M) {
    if (!(looks_like_number($opt_M))) {
	die ("ERROR: $opt_M is not a number(-M)!\n");
    }
    if ($opt_M < 0) {
	die ("ERROR: $opt_M is not > 0 (-M)!\n");
    }
    $max_thresh = $opt_M;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - cluster_size.pl reads a PanOCT matchtable_0_1 file (-T) and a genome tag file (-G) and outputs a
         summary of cluster sizes to stdout for the genomes based on the thresholds specified as either
	 numbers (-m, -M) or percentages (-p, -P). A histogram of cluster size is also output for each
	 genome in a separate file in a user specified directory (-d).

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
Example: cluster_size.pl -G tags -d cluster_hists -T matchtable_0_1.txt -p 10 -P 90 > cluster_size.txt
 Option:
     -h: print this help page
     -G: genome tage file name
     -T: matchtable_0_1 file name
     -d: output file directory
     -M: maximum threshold number
     -m: minimum threshold number
     -P: maximum threshold pecentage: if -M is specified this is ignored (default 95)
     -p: minimum threshold pecentage: if -m is specified this is ignored (default 5)
     -D: DEBUG MODE (DEFAULT = off)
 Output: summary of cluster sizes and hitogram of cluster sizes for each genome
 Author: Granger Sutton, Ph.D.
 Date: 2/5/14
_EOB_

exit;
}

sub get_genomes {  # obtain list of genomes - must be in the same order as the matchtable columns
   
    my %temp_hash = ();
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: can't open file $genomes_file_name\n");
    print "Order of genomes in $genomes_file_name with array index\n";
    while (<$infile>)  {
	chomp;
	if (length($_) > 0) {
	    my $name = $_;
	    $name =~ s/\s+$//;

            if (defined $temp_hash{$name})  {
               die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
            } else  {
		$temp_hash{$name} = 1;
		push (@genome_array, $name); # populate the genome_array in the order of the genome file
		$genome_hists[$genome_number] = [];
		$genome_min[$genome_number] = 0;
		$genome_max[$genome_number] = 0;
		print "$name\t$genome_number\n";
		$genome_number++;
	    }
	}  
    }
    close($infile);
    print "$genome_number genomes\n\n";
    if (!defined $max_thresh) {
	$max_thresh = ($genome_number * $max_perc) / 100;
    }
    if (!defined $min_thresh) {
	$min_thresh = ($genome_number * $min_perc) / 100;
    }
    for (my $index1 = 0; $index1 <= $genome_number; $index1++) {
	for (my $index2 = 0; $index2 <= $genome_number; $index2++) {
	    $genome_hists[$index1]->[$index2] = 0;
	}
    }
}

sub get_matchtable_rows {

    unless (open (MATCHTABLEFILE, $matchtable_file_name) )  {
	die ("ERROR: can not open file $matchtable_file_name.\n");
    }
    while (<MATCHTABLEFILE>) {
	my $number = "";
	my @matchtable_row = ();
	my $cluster_size = 0;
	chomp;
	($number, @matchtable_row) = split(/\t/, $_);  # split the scalar $matchtable_row on tab
	if ($number eq "") {
	    die ("there must be a cluster number on each line of the matchtable file: $matchtable_file_name\n");
	}
	if (!(looks_like_number($number))) {
	    die ("ERROR: $number is not a number for cluster number in matchtable file!\n");
	}
	if ($number <= 0) {
	    die ("ERROR: $number is not > 0 for cluster number in matchtable file!\n");
	}
	if ($#matchtable_row != $#genome_array) {
	    die ("ERROR: genome columns in $matchtable_file_name is not equal to genomes in $genomes_file_name!\n");
	}
	for (my $index = 0; $index <= $#matchtable_row; $index++) {
	    $cluster_size += $matchtable_row[$index];
	}
	for (my $index = 0; $index <= $#matchtable_row; $index++) {
	    if ($matchtable_row[$index] > 0) {
		$genome_hists[$index]->[$cluster_size]++;
		if ($cluster_size >= $max_thresh) {
		    $genome_max[$index]++;
		}
		if ($cluster_size <= $min_thresh) {
		    $genome_min[$index]++;
		}
	    }
	}
    }
    close (MATCHTABLEFILE);
    return;
}

sub print_summary {

    print "Number of clusters per genome >= $max_thresh sorted in ascending order\n\n";
    foreach my $index (sort {$genome_max[$a] <=> $genome_max[$b]} (0 .. $#genome_max))  {
	print "$genome_array[$index]\t$genome_max[$index]\t$genome_min[$index]\n";
    }
    print "\nNumber of clusters per genome <= $min_thresh sorted in descending order\n";
    foreach my $index (sort {$genome_min[$b] <=> $genome_min[$a]} (0 .. $#genome_min))  {
	print "$genome_array[$index]\t$genome_min[$index]\t$genome_max[$index]\n";
    }
    $hist_dir =~ s/\/*$//; #remove trailing / if needed
    if (-d "$hist_dir") { #directory already exists
    } elsif (-e "$hist_dir") { #file instead of directory
	die ("ERROR: $hist_dir exists but is not a directory\n");
    } else { #try to create a directory
	if (!mkdir($hist_dir)) { #failed to create directory
	    die ("ERROR: could not create directory $hist_dir\n");
	}
    }

    for (my $index1 = 0; $index1 < $genome_number; $index1++) {
	my $file_name = $hist_dir . "/" . $genome_array[$index1];
	open (my $outfile, ">", "$file_name") || die ("ERROR: can't open file $file_name\n");
	for (my $index2 = 1; $index2 <= $genome_number; $index2++) {
	    print $outfile "$index2\t$genome_hists[$index1]->[$index2]\n";
	}
	close ($outfile);
    }
    return;
}


########################################  M A I N  #####################################################
print STDERR "Getting genome names from $genomes_file_name\n";
&get_genomes;
print STDERR "Getting matchtable rows from $matchtable_file_name\n";
&get_matchtable_rows;
print STDERR "Printing cluster size summary to STDOUT\n";
&print_summary;
exit(0);
