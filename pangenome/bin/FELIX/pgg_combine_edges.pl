#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use strict;
my $genome_number;
my %pgg_edges = ();  # key = pgg core edge, value = array of 0/1 of size number of genomes
my %name_to_cluster = (); # key = gene name, value = cluster number gene is in

my $help_text = "THIS IS PLACEHOLDER TEXT\n";                             #<------------------------------------------------------------------------ TODO: write help text


GetOptions('index=s' => \my $index_file,
	   'help' => \my $help);

if ($help)
{
    print STDERR "$help_text";
    exit;
}

sub read_index {

############ Input a two column tab delimited file with gene_name and cluster number ###############
    
    open (my $indexfile, "<", $index_file) || die ("ERROR: can't open index file $index_file\n");
    while (my $line = <$indexfile>) {
	chomp $line;
	(my $gene_name, my $cluster) = split(/\t/, $line);  # split the scalar $line on tab
	$name_to_cluster{$gene_name} = $cluster;
	#print "$gene_name $cluster $name_to_cluster{$gene_name}\n";
    }
    close ($indexfile);
    return;
}

sub get_edges {  # obtain list of edge files - order will determine column order in PGG output
   
    $genome_number = 0;     # total number of genomes to be processed
    my @zero_vector = ();   # vector of 0/1 for PGG columns
    my @edgefiles = ();     # array of edgefile names

    while (my $edgefile = <STDIN>)  {
	chomp $edgefile;
	$zero_vector[$genome_number] = 0;
	$edgefiles[$genome_number] = $edgefile;
	$genome_number++;
    }
    my $index = 0;
    foreach my $edgefile (@edgefiles) {
	my $edgehandle;
	unless (open ($edgehandle, "<", $edgefile) )  {
	    die ("ERROR: cannot open file $edgefile.\n");
	}
	while (my $edge = <$edgehandle>) {
	    chomp $edge;
	    my $cluster1;
	    my $cluster2;
	    my $whichend1;
	    my $whichend2;
	    if ($edge =~ /\((.+)_([35]),(.+)_([35])\)/) {
		$cluster1 = $1;
		$cluster2 = $3;
		$whichend1 = $2;
		$whichend2 = $4;
	    } else {
		die ("ERROR: Bad edge formatting $edge in file $edgefile.\n");
	    }
	    if (defined $name_to_cluster{$cluster1}) {
		$cluster1 = $name_to_cluster{$cluster1};
	    }
	    if (defined $name_to_cluster{$cluster2}) {
		$cluster2 = $name_to_cluster{$cluster2};
	    }
	    $edge = "(" . $cluster1 . "_" . $whichend1 . "," . $cluster2 . "_" . $whichend2 . ")";
	    if (defined $pgg_edges{$edge}) {
		$pgg_edges{$edge}[$index] = 1;
	    } else {
		$pgg_edges{$edge} = [];
		@{ $pgg_edges{$edge} } = @zero_vector;
		$pgg_edges{$edge}[$index] = 1;
	    }
	}
	close ($edgehandle);
	$index++;
    }
    return;
}

sub print_edges {  # sort and print out the combined PGG edges
    my $sort_edges = sub { # sort by cluster1 then cluster1 end then cluster2 then cllustr2 end
	if ($a =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
	    my $a_cluster1 = $1;
	    my $a_cluster2 = $3;
	    my $a_end1 = $2;
	    my $a_end2 = $4;
	    if ($b =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
		my $b_cluster1 = $1;
		my $b_cluster2 = $3;
		my $b_end1 = $2;
		my $b_end2 = $4;
		if ($a_cluster1 <=> $b_cluster1) {
		    return ($a_cluster1 <=> $b_cluster1);
		} elsif ($a_end1 <=> $b_end1) {
		    return ($a_end1 <=> $b_end1);
		} elsif ($a_cluster2 <=> $b_cluster2) {
		    return ($a_cluster2 <=> $b_cluster2);
		} else {
		    return ($a_end2 <=> $b_end2);
		}
	    } else {
		die ("ERROR: Bad edge formatting $b.\n");
	    }
	} else {
	    die ("ERROR: Bad edge formatting $a.\n");
	}
    };

    foreach my $edge (sort $sort_edges (keys %pgg_edges)) {
	print STDOUT "$edge";
	foreach my $value ( @{ $pgg_edges{$edge} } ) {
	    print STDOUT "\t$value";
	}
	print STDOUT "\n";
    }
    return;
}

########################################  M A I N  #####################################################
if ($index_file) {
    &read_index;
}
&get_edges;
&print_edges;
exit(0);
