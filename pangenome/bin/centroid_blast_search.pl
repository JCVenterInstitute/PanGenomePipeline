#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my %blast_results = ();
my $help_text = "This program BLASTs a FASTA file of centroids against a nucleotide space,
groups the results by the molecules of the nucleotide space, and then sorts
these groups by the midpoint of the matched region from closest to furthest.
This work is done in a folder called C_BLAST_TMP, which is deleted when the 
program finishes.

Input Flags:
-centroids - The FASTA file containing the centroids (required)
-nspace - A FASTA file with representing the nucleotide space (required)
-help - Outputs this help text";

GetOptions('centroids=s' => \my $centroids,
	'nspace=s' => \my $nucleotide_space,
	'help' => \ my $help);
	
if($help){
	print("$help_text\n");
}
	
elsif(!$nucleotide_space or !$centroids){
	print("Error: One or both of the required flags are missing\n$help_text\n");
}

else{
	`mkdir C_BLAST_TMP`;
	`cp $nucleotide_space C_BLAST_TMP/temp_fasta.ftmp`;
	`makeblastdb -in C_BLAST_TMP/temp_fasta.ftmp -dbtype nucl -out C_BLAST_TMP/temp_fasta.ftmp`;
	`blastn -query $centroids -db C_BLAST_TMP/temp_fasta.ftmp -out C_BLAST_TMP/temp_results.ftmp -outfmt 6`;

	open(BLAST, "C_BLAST_TMP/temp_results.ftmp");

	while(my $line = <BLAST>){
		#print("$line");
		my @blast_data = split(/\t/, $line);
		my $midpoint = int((($blast_data[9] + $blast_data[8]) / 2));
		$line =~ s/\n/\t$midpoint\n/;
		
		if(!($blast_results{$blast_data[1]})){
			$blast_results{$blast_data[1]} = ();
		}
		
		$blast_results{$blast_data[1]}{$line} = $midpoint;
	}

	close(BLAST);

	open(OUT, ">blast_results");

	print OUT "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, s. midpoint\n";
	foreach my $molecule(sort keys(%blast_results)){
		my %molecule_hash = %{$blast_results{$molecule}};
		foreach my $blast_info(sort{$molecule_hash{$a} <=> $molecule_hash{$b}} keys(%molecule_hash)){
			print OUT "$blast_info";
		}
	}

	`rm -r C_BLAST_TMP`;

	close(OUT);
}