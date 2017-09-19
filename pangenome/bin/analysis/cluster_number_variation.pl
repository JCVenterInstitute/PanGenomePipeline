#!/usr/bin/env perl

###############################################################################
#                                                                             #
#       Copyright (C) 2016-2017 J. Craig Venter Institute (JCVI).             #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

##---------------------------------------------------------------------------##
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
##---------------------------------------------------------------------------##

my $dir        = shift;
my $genome_dir = shift;
my %genomes;
opendir( DIR, $genome_dir ) || die "Can't open $genome_dir\n";

while ( my $file = readdir(DIR) ) {
	next unless $file =~ /(\S+)\.pep/;
	my $genome = $1;
	my $full_name = File::Spec->catfile( $genome_dir, $file );
	open( FILE, $full_name ) || die "Can't open $full_name\n";

	while ( my $line = <FILE> ) {
		next unless $line =~ /^\>(\S+)\s+/;
		my $tagname = $1;

		if ( $tagname =~ /\|/ ) {
			my @f = split( /\|/, $tagname );
			$tagname = $f[0];
		}

		if ( exists $genomes{$tagname} ) {
			die "Duplicate $tagname exists in $file\n";
		}
		$genomes{$tagname}->{$genome} = 1;
	}
	close(FILE);
}
close(DIR);

opendir( DIR, $dir ) || die "Can't open $dir\n";
print "Method\tNum_cluster\tNum_in_cluster\tDelta\n";

while ( my $file = readdir(DIR) ) {
	next unless $file =~ /^parsed_(\S+)\.txt$/;
	my $method = $1;
	my $full_path = File::Spec->catfile( $dir, $file );
	my ( $num_cluster, $total_number ) = read_data($full_path);
	print "$method\t$num_cluster\t$total_number\t",
		( $total_number - $num_cluster ), "\n";

}
close(DIR);

sub read_data {
	my $file = shift;
	open( FILE, $file ) || die "Can't open $file\n";
	my $number_cluster = 0;
	my $total_number   = 0;

	while ( my $line = <FILE> ) {
		chomp $line;
#		$number_cluster++;
		my ( $cluster, @f ) = split( /\t/, $line );
		if ( scalar(@f) != 4 ) { next; }
		my %num_genomes;

		foreach my $i (@f) {
			die "$i not found in $file" unless exists $genomes{$i};
			$num_genomes{ $genomes{$i} } = 1;
		}

		if ( keys(%num_genomes) == 4 ) {

			$total_number += ( scalar(@f) );
			$number_cluster++;
		}
	}
	close(FILE);
	return $number_cluster, $total_number;
}

