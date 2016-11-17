#!/usr/bin/env perl
#Copyright (C) 2014-2015  The J. Craig Venter Institute (JCVI).  All rights reserved
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
getopts ('DhVM:P:');# M is matchtable_0_1 file, P is paralogs file - both required
our ($opt_h, $opt_D, $opt_M, $opt_P, $opt_V);

my $matchtable_file_name;
my $paralog_file_name;
my %paralog_cluster = ();
my %paralog_matchtable_row = ();
my $DEBUG = 0;
my $version = "1.0";

if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_V) {die "$prog version $version\n";}
if (($opt_M) && (-s "$opt_M")) {
    $matchtable_file_name = $opt_M;
} else {
    print STDERR "Error with -M\n";
    &option_help;
}
if (($opt_P) && (-s "$opt_P")) {
    $paralog_file_name = $opt_P;
} else {
    print STDERR "Error with -P\n";
    &option_help;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog:

           Takes a PanOCT matchtable_0_1 and a paralogs file (output from PanOCT)
           and outputs a parlogs matchtable where instead of 0 or 1 for presence or 
           absence of the cluster a number is given indicating the number of paralogs. 
           This is accomplished by adding together the rows in the matchtable_0_1, 
           which are paralogs from the paralogs file. Only the lowest numbered paralog 
           cluster is output.

Copyright (C) 2014-2015  The J. Craig Venter Institute (JCVI).  All rights reserved

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

Citation:  Agnes P. Chan, Granger Sutton, Jessica DePew, Radha Krishnakumar, Yongwook Choi, Xiao-Zhe Huang,
           Erin Beck, Derek M. Harkins, Maria Kim, Emil P. Lesho, Mikeljon P. Nikolich and Derrick E. Fouts 
           (2015) "A novel method of consensus pan-chromosome assembly and large-scale comparative analysis 
           reveal the highly flexible pan-genome of Acinetobacter baumannii" Genome Biol. 16(143):1-28.

  Usage: $prog <options>
Example: $prog -M matchtable_0_1.txt -P paralogs.txt > matchtable_paralog.txt
Version: $version
Options:
     -h: print this help page
     -M: matchtable_0_1 file from PanOCT, format is tab delimited with cluster number followed by 0s or 1s
     -P: paralogs file from PanOCT, tab delimited, cluster numbers - one line for each group of paralogs
     -D: DEBUG MODE (DEFAULT = off)
 Output: rows for each paralog group from the matchtable_0_1 file are added together and output
     -V: print just (V)ersion information
 Author: Granger Sutton, Ph.D.
 Date: 09/04/15
_EOB_
   exit;
}

sub get_paralogs {

    my $number = "";
    my $second = "";

    unless (open (PARALOGFILE, $paralog_file_name) )  {
	die ("ERROR: can not open file $paralog_file_name.\n");
    }
    while (<PARALOGFILE>) {
	my @paralog_line = ();
	chomp;
	@paralog_line = split(/\t/, $_);  # split the scalar $paralog_line on tab
	$number = $paralog_line[0];
	$second = $paralog_line[1];
	if (($number eq "") || ($second eq "")) {
	    die ("there must be at least two paralog cluster numbers per line in paralog file: $paralog_file_name\n");
	}
	if (!(looks_like_number($number))) {
	    die ("ERROR: $number is not a number for paralog cluster number in paralog file!\n");
	}
	if ($number <= 0) {
	    die ("ERROR: $number is not > 0 for paralog cluster number in paralog file!\n");
	}
	foreach my $paralog (@paralog_line) {
	    $paralog_cluster{$paralog} = $number; #have all paralog clusters in the group point to the first paralog cluster
	}
    }
    close (PARALOGFILE);
    return;
}

sub get_matchtable_rows {

    my $number = "";

    unless (open (MATCHTABLEFILE, $matchtable_file_name) )  {
	die ("ERROR: can not open file $matchtable_file_name.\n");
    }
    while (<MATCHTABLEFILE>) {
	my @matchtable_line = ();
	chomp;
	($number, @matchtable_line) = split(/\t/, $_);  # split the scalar $matchtable_line on tab
	if ($number eq "") {
	    die ("there must be a cluster number on each line of the matchtable file: $matchtable_file_name\n");
	}
	if (!(looks_like_number($number))) {
	    die ("ERROR: $number is not a number for cluster number in matchtable file!\n");
	}
	if ($number <= 0) {
	    die ("ERROR: $number is not > 0 for cluster number in matchtable file!\n");
	}
	if (!defined $paralog_cluster{$number}) {
	    next;
	}
	$number = $paralog_cluster{$number};
	if (!defined $paralog_matchtable_row{$number}) {
	    $paralog_matchtable_row{$number} = [];
	    for (my $index = 0; $index <= $#matchtable_line; $index++) {
		$paralog_matchtable_row{$number}->[$index] = $matchtable_line[$index];
	    }
	} else {
	    for (my $index = 0; $index <= $#matchtable_line; $index++) {
		$paralog_matchtable_row{$number}->[$index] += $matchtable_line[$index];
	    }
	}
    }
    close (MATCHTABLEFILE);
    return;
}

sub print_paralog_matchtable {

    my $number = "";

    unless (open (MATCHTABLEFILE, $matchtable_file_name) )  {
	die ("ERROR: can not open file $matchtable_file_name.\n");
    }
    while (<MATCHTABLEFILE>) {
	my @matchtable_line = ();
	chomp;
	($number, @matchtable_line) = split(/\t/, $_);  # split the scalar $matchtable_line on tab
	if (!defined $paralog_cluster{$number}) {
	    print "$_\n";
	} elsif (defined $paralog_matchtable_row{$number}) {
	    my $matchtable_line = join("\t", @{ $paralog_matchtable_row{$number} });
	    print "$number\t$matchtable_line\n";
	}
    }
    close (MATCHTABLEFILE);
    return;
}


########################################  M A I N  #####################################################
print STDERR "Getting paralog clusters from $paralog_file_name\n";
&get_paralogs;
print STDERR "Getting matchtable rows from $matchtable_file_name\n";
&get_matchtable_rows;
print STDERR "Printing paralog matchtable to STDOUT\n";
&print_paralog_matchtable;
exit(0);
