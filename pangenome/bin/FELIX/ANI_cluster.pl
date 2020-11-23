#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;

my $help_text = "This program reads in a 3 column tab delimitted file from standard in. The 3 columns are the first three from mash dist runs: name1, name2, distance. This program outputs to files in the current directory single-linkage clusters using the threshold. The files are named after the founding genome and have suffix .cluster.

Input Flags:
-threshold - used for column 3 to see if two genomes should be clustered <= (default 0.02)
-help - outputs this help text";

GetOptions('threshold=s' => \my $threshold,
	   'help' => \my $help);
	
if(!$threshold) {
    $threshold = 0.02;
}
	
if($help){
    print("$help_text\n");
    exit;
}
	
#Globals
my %cluster_name = ();         # key = name, value is the name of the founding cluster - with adjustments when clusters are merged.

#####################################################################################################
sub read_input    # Read in the 3 column input (name1, name2, distance)
{
    my $line;
    while ($line = <STDIN>) {
	(my $name1, my $name2, my $distance) = split(/\t/, $line, 3);
	if (($name1 ne $name2) && ($distance <= $threshold)) {
	    if (!defined($cluster_name{$name1})) {
		if (!defined($cluster_name{$name2})) {
		    $cluster_name{$name1} = $cluster_name{$name2} = $name1;
		} else {
		    $cluster_name{$name1} = $cluster_name{$name2};
		}
	    } else {
		if (!defined($cluster_name{$name2})) {
		    $cluster_name{$name2} = $cluster_name{$name1};
		} else {
		    if ($cluster_name{$name1} ne $cluster_name{$name2}) {
			my $head1 = $cluster_name{$name1};
			while ($head1 ne $cluster_name{$head1}) {
			    $head1 = $cluster_name{$head1};
			}
			my $head2 = $cluster_name{$name2};
			while ($head2 ne $cluster_name{$head2}) {
			    $head2 = $cluster_name{$head2};
			}
			$cluster_name{$head2} = $head1;
		    }
		}
	    }
	}
    }
	
    return;
}

#####################################################################################################
sub write_cluster_files    # Write out the single-linkage cluster files with names of the head genome and suffix .cluster
{
    foreach my $name (keys %cluster_name)  { # go through all the nodes of the cluster trees
	my $head = $cluster_name{$name};
	while ($head ne $cluster_name{$head}) {
	    $head = $cluster_name{$head};
	}
	my $cluster_file = $head . ".cluster";
	my $cluster_handle;
	unless (open ($cluster_handle, ">>", $cluster_file) )  {
	    die ("cannot open file $cluster_file!\n");
	}
	print $cluster_handle "$name\n";
	close($cluster_handle);
    }
	
    return;
}

{#main
    &read_input;
    &write_cluster_files;
}
