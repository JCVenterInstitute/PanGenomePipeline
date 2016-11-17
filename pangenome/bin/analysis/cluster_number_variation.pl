#!/usr/bin/env perl
# $Id: cluster_number_variation.pl 32957 2010-08-05 16:22:09Z mbasu $

##---------------------------------------------------------------------------##
##  File: cluster_number_variation.pl
##
##  Author:
##        Malay <malay@bioinformatics.org>
##
##  Description:
##
#******************************************************************************
#* Copyright (C) 2010 Malay K Basu <malay@bioinformatics.org>
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

cluster_number_variation.pl - One line description.

=head1 SYNOPSIS

cluster_number_variation.pl [options] -o <option>


=head1 DESCRIPTION

Write a description of your prgram. 


=head1 ARGUMENTS 

=over 4

=item B<--option|-o>

First option.



=back

=head1 OPTIONS

Something here.


=head1 SEE ALSO

=head1 COPYRIGHT

Copyright (c) 2010 Malay K Basu <malay@bioinformatics.org>

=head1 AUTHORS

Malay K Basu <malay@bioinformatics.org>

=cut

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

