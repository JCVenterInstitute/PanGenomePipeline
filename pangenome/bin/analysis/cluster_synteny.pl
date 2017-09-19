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
use SeqToolBox::Interval::Tree;
use SeqToolBox::Interval::Range;
use Getopt::Long;
use Data::Dumper;
use DBI;

my $coord_db    = "";
my $cluster_dir = "";
my $peptide_dir = "";

GetOptions( "coord|o=s"    => \$coord_db,
			"clusters|c=s" => \$cluster_dir,
			"genomes|g=s"  => \$peptide_dir
) || die "Missing parameter\n";

my $coord;
my %genomes;
my $cluster_members;
my $clusters;
my $clusters_only_4;
my $dbh = DBI->connect( "dbi:SQLite:dbname=$coord_db", "", "",
						{ AutoCommit => 1, RaiseError => 1 } );
my $left_node
	= $dbh->prepare(
	'select previous_node_name from synteny where organism = ? and node_name = ?'
	);
my $right_node
	= $dbh->prepare(
	 'select next_node_name from synteny where organism = ? and node_name = ?'
	);

#my $coord_genome_dir = shift;

#print STDERR "INPUT: $coord_genome_dir\t$cluster_dir\t$peptide_dir\n";
#opendir( DIR, $coord_genome_dir ) || die "Can't open $coord_genome_dir\n";
#print STDERR "After opening\n";
#
#while ( my $file = readdir(DIR) ) {
#	print STDERR $file, "\n";
#	next unless $file =~ /(\S+)\.coords/;
#	my $basename = $1;
#	$coord->{$basename}
#		= read_coord( File::Spec->catfile( $coord_genome_dir, $file ) );
#		print STDERR "KEYS: " , keys %{$coord}, "\n";
#
#}
#close DIR;
#
#sub read_coord {
#	my $filename = shift;
#	open( FILE, $filename ) || die "Can't open $filename\n";
#	my $tree = SeqToolBox::Interval::Tree->new();
#
#	while ( my $line = <FILE> ) {
#		chomp $line;
#		my ( $name, $start, $end, $dummy ) = split( /\t/, $line );
#
#		if ( $start > $end ) {
#			my $tmp = $start;
#			$start = $end;
#			$end   = $tmp;
#
#		}
#		my $range = SeqToolBox::Interval::Range->new( $name, $start, $end );
#		$tree->insert($range);
#	}
#
#	return $tree;
#}

#my $dir        = shift;
#my $genome_dir = shift;

opendir( DIR, $peptide_dir ) || die "Can't open $peptide_dir\n";

while ( my $file = readdir(DIR) ) {
	next unless $file =~ /(\S+)\.pep/;
	my $genome = $1;
	my $full_name = File::Spec->catfile( $peptide_dir, $file );
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
		$genomes{$tagname} = $genome;
	}
	close(FILE);
}
close(DIR);

opendir( DIR, $cluster_dir ) || die "Can't open $cluster_dir\n";

#print "Method\tNum_cluster\tNum_in_cluster\tDelta\n";

while ( my $file = readdir(DIR) ) {
	next unless $file =~ /^parsed_(\S+)\.txt$/;
	my $method       = $1;
	my $full_path    = File::Spec->catfile( $cluster_dir, $file );
	my $total_number = read_data( $full_path, $method );

	#print "$method\t$num_cluster\t$total_number\t",
	#	( $total_number - $num_cluster ), "\n";
	#	print $method, "\t", $total_number, "\n";
}
close(DIR);

foreach my $m ( keys %{$clusters_only_4} ) {
	my $total_count = 0;

	foreach my $c ( keys %{ $clusters_only_4->{$m} } ) {

		#		print STDERR "$m\t$c\n";

		my @proximal_genes;
		my @distal_genes;

		foreach my $g ( keys %{ $clusters_only_4->{$m}->{$c} } ) {
			my $genome = get_genome_name($g);
			next unless $genome;
			my $proximal = get_proximal_node( $genome, $g );
			next unless $proximal;

			#			my $proximal_cluster_name = get_cluster_name ($proximal,$m);
			#			next unless $proximal_cluster_name;

			push @proximal_genes, $proximal;
			my $distal = get_distal_node( $genome, $g );
			next unless $distal;
			push @distal_genes, $distal;

		}

		my $max_minus_one = get_count( \@proximal_genes, $m );
		my $max_plus_one  = get_count( \@distal_genes,   $m );
		$total_count += $max_minus_one;
		$total_count += $max_plus_one;
	}

	print "RESULT: ", $m, "\t", $total_count, "\n";
}

$left_node->finish();
$right_node->finish();
$dbh->disconnect();

sub get_count {
	my ( $array, $method ) = @_;
	my $max_count = 0;

	foreach my $g ( @{$array} ) {
		my $cluster = get_cluster_name( $g, $method );
		next unless $cluster;
		my $members = get_cluster_members( $method, $cluster );
		next unless $members;
		my $count = find_member_count( $members, $array );

		if ( $count > $max_count ) {
			$max_count = $count;
		}
	}
	return $max_count;
}

sub find_member_count {
	my ( $members, $array ) = @_;
	my $count = 0;

	foreach my $i ( @{$array} ) {
		if ( exists $members->{$i} ) {
			$count++;
		}
	}
	return $count;
}

sub get_cluster_name {
	my $name = shift;

	if ( $name eq 'NA' ) {
		return;
	}
	my $method = shift;

	if ( exists $cluster_members->{$method}->{$name} ) {
		return $cluster_members->{$method}->{$name};
	}
	else {

		#		print STDERR "Cluster name not found for $method and $name\n";
		#		print Dumper($cluster_members->{$method});
		print STDERR "Cluster name not found for $method and $name\n";
		return;
	}
}

sub get_cluster_members {
	my ( $method, $cluster ) = @_;

	if ( exists $clusters->{$method}->{$cluster} ) {
		return $clusters->{$method}->{$cluster};
	}
	else {
		print STDERR "Members not found for $cluster in $method\n";
		return;
	}
}

sub get_distal_node {
	my ( $genome, $g ) = @_;

	#	print STDERR keys %{$coord};
	$right_node->execute( $genome, $g );
	my $proximal;

	while ( my @row = $right_node->fetchrow_array ) {
		$proximal = $row[0];
	}
	return $proximal ? $proximal : undef;
}

sub get_proximal_node {
	my ( $genome, $g ) = @_;

	#	print STDERR keys %{$coord};
	$left_node->execute( $genome, $g );
	my $proximal;

	while ( my @row = $left_node->fetchrow_array ) {
		$proximal = $row[0];
	}
	return $proximal ? $proximal : undef;
}

sub get_genome_name {
	my $g = shift;

	unless ( exists $genomes{$g} ) {
		print STDERR "Genome name for $g not found\n";
		return;
	}
	return $genomes{$g};
}

sub read_data {
	my $file   = shift;
	my $method = shift;
	open( FILE, $file ) || die "Can't open $file\n";
	my $number_cluster = 0;
	my $total_number   = 0;

	while ( my $line = <FILE> ) {
		chomp $line;

		#		$number_cluster++;
		my ( $cluster, @f ) = split( /\t/, $line );

		#		if ( scalar(@f) != 4 ) { next; }
		my %num_genomes;

		my %members;

		foreach my $i (@f) {

			die "$i not found in $file" unless exists $genomes{$i};

			#			print $genomes{$i}, "\n";
			$num_genomes{ $genomes{$i} }      = 1;
			$cluster_members->{$method}->{$i} = $cluster;
			$members{$i}                      = 1;

		}

		$clusters->{$method}->{$cluster} = \%members;

		if ( keys(%num_genomes) == 4 ) {
			$clusters_only_4->{$method}->{$cluster} = \%members;
			$total_number += ( scalar(@f) );
			$number_cluster++;
		}
	}
	close(FILE);
	return $number_cluster, $total_number;
}
