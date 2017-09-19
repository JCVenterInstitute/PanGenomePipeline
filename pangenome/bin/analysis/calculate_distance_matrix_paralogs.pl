#!/usr/bin/env perl

###############################################################################
#                                                                             #
#       Copyright (C) 2010 J. Craig Venter Institute (JCVI).                  #
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

my ( $ref_genome, $blast_result_file, @files ) = @ARGV;

my $data;
my $clusters;
my $results;
my $blast_result;

read_blast_result($blast_result_file);

foreach my $file (@files) {
	read_data($file);
}

my @tags = keys( %{$data} );
open( my $fh, $ref_genome ) || die "Can't open $ref_genome\n";

while ( my $line = <$fh> ) {
	next unless $line =~ /^\>(\S+)\|/;
	my $gi = $1;
	die "Can't parse gi" unless $gi;
	calculate_similarity($gi);
}
close($fh);

print "Method1\tMethod2\tScore\tMatched_cluster\tSimilarity\n";

foreach my $i (@tags) {
	foreach my $j (@tags) {
		print $i, "\t", $j, "\t",
			sprintf( "%.2f", $results->{$i}->{$j}->{result} ), "\t",
			$results->{$i}->{$j}->{number}, "\t",
			sprintf( "%.3f",
					       $results->{$i}->{$j}->{result}
						 / $results->{$i}->{$j}->{number} ),
			"\n";
	}
}

sub read_blast_result {
	my $file = shift;
	open( my $infile, $file ) || die "Can't open $file\n";

	while ( my $line = <$infile> ) {
		next if $line =~ /^\#/;
		chomp $line;
		my @f = split( /\t/, $line );
		my $q;
		my $s;

		if ( $f[0] =~ /(\S+)\|/ ) {
			$q = $1;
		}
		else {
			$q = $f[0];
		}

		if ( $f[1] =~ /(\S+)\|/ ) {
			$s = $1;
		}
		else {
			$s = $f[1];
		}
		my $score = $f[11];
		$blast_result->{$q}->{$s} = $score;
	}
}

#sub calculate_similarity {
#	my $gi = shift;
#
#	for (my $i = 0; $i < scalar(@tags); $i++) {
#		for (my $j = 0; $j <scalar (@tags); $j++){
#		my ($num_1, $num_2);
#
#		if (!exists $data->{$tags[$i]}->{$gi} || !exists $data->{$tags[$j]}->{$gi}) {
#			next;
#		}
#		my $cluster_1 = get_clusters($tags[$i], $gi);
#		my $cluster_2 = get_clusters ($tags[$j], $gi);
#		my $sim = get_similar_count ($cluster_1, $cluster_2);
#		my $count_1 = scalar(keys %{$cluster_1});
#		my $count_2 = scalar(keys %{$cluster_2});
#		my $result = (($sim/$count_1) + ($sim/$count_2))/2;
#		$results->{$tags[$i]}->{$tags[$j]}->{result} += $result;
#		$results->{$tags[$i]}->{$tags[$j]}->{number}++;
#
#	}
#	}
#}

sub calculate_similarity {
	my $gi = shift;

	unless ( exists $blast_result->{$gi} ) {
		print STDERR "Can't find blast_result for $gi\n";
		return;
	}

	die "Can't find self score for $gi\n"
		unless exists $blast_result->{$gi}->{$gi};

	my $self_score = $blast_result->{$gi}->{$gi};
	for ( my $i = 0; $i < scalar(@tags); $i++ ) {

		for ( my $j = 0; $j < scalar(@tags); $j++ ) {
			my ( $num_1, $num_2 );

			if ( !exists $data->{ $tags[$i] }->{$gi} ) {
				next;
			}
			my $cluster_1 = get_clusters( $tags[$i], $gi );
			my $cluster_2 = get_clusters( $tags[$j], $gi );
			my ( $diff_1, $diff_2 )
				= get_difference( $cluster_1, $cluster_2 );
			my $avg_par_score1 = get_avg_par_score( $diff_1, $gi );
			my $avg_par_score2 = get_avg_par_score( $diff_2, $gi );

			#		my $sim = get_similar_count ($cluster_1, $cluster_2);
			#		my $count_1 = scalar(keys %{$cluster_1});
			#		my $count_2 = scalar(keys %{$cluster_2});
			my $result = ( $avg_par_score1 + $avg_par_score2 ) / 2;
			$results->{ $tags[$i] }->{ $tags[$j] }->{result} += $result;
			$results->{ $tags[$i] }->{ $tags[$j] }->{number}++;

		}
	}
}

sub get_avg_par_score {
	my ( $diff_1, $gi ) = @_;
	my $count       = 0;
	my $total_score = 0;

	foreach my $i ( %{$diff_1} ) {
		if ( exists $blast_result->{$gi}->{$i} ) {
			$total_score += $blast_result->{$gi}->{$i};
			$count++;
		}
	}
	if ( $count == 0 ) { return 0; }
	return $total_score / $count;
}

sub get_difference {
	my ( $array1, $array2 ) = @_;
	my $diff_1;
	my $diff_2;

	foreach my $i ( keys %{$array1} ) {
		if ( exists $array2->{$i} ) {
			next;
		}
		else {
			$diff_1->{$i} = 1;
		}
	}

	foreach my $i ( keys %{$array2} ) {
		if ( exists $array1->{$i} ) {
			next;
		}
		else {
			$diff_2->{$i} = 1;
		}
	}
	return $diff_1, $diff_2;
}

sub get_similar_count {
	my ( $array1, $array2 ) = @_;
	my $match_count = 0;

	foreach my $i ( keys %{$array1} ) {
		if ( exists $array2->{$i} ) {
			$match_count++;
		}
	}
	return $match_count;

}

sub get_clusters {
	my ( $tag, $gi ) = @_;
	my @clusters = keys( %{ $data->{$tag}->{$gi} } );
	my @super;
	my %r;

	foreach my $c (@clusters) {
		push @super, @{ $clusters->{$tag}->{$c} };
	}

	foreach my $s (@super) {
		$r{$s} = 1;
	}
	return \%r;
}

sub read_data {
	my $file = shift;
	my $tag  = $file;
	print STDERR $tag, "\n";
	$tag =~ s/\S*\///;
	$tag =~ /parsed\_(\S+)\.txt/;
	my $tag_name = $1;
	die "Could not parse the tag from filename" unless $tag_name;
	open( my $infile, $file ) || die "Can't open $file\n";

	while ( my $line = <$infile> ) {
		chomp $line;
		my ( $cluster, @genes ) = split( /\t/, $line );
		$clusters->{$tag_name}->{$cluster} = \@genes;

		foreach my $g (@genes) {
			if ( exists $data->{$tag_name}->{$g} ) {
				print STDERR
					"$g found in more than one cluster in $tag_name\n";
				$data->{$tag_name}->{$g}->{$cluster} = 1;
			}
			else {
				$data->{$tag_name}->{$g}->{$cluster} = 1;
			}
		}
	}
	close($infile);
}

