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

use warnings;
use strict;
$|++;

=head1 NAME

get_dna_for_loci.pl - retrieve genomic DNA for specified loci.

=head1 SYNOPSIS

    USAGE: get_dna_for_loci.pl -a <att_file> -f <genome_fasta> -l <loci_list_file> [ -o <output_file> ]

=head1 OPTIONS

=head2 Required Options

B<--att_file, -a>       :   Path to gene attribute file.

B<--genome_fasta, -f>   :   Path to genome fasta file.

B<--loci_list, -l>      :   Path to file containing one-column list of locus IDs.

=head2 Optional Options

B<--output_file, -o>    :   To send output to a file.  Otherwise, it will go to STDOUT.

B<--region, -r>         :   Excise a single region from the collective lowest coordinate to the collective highest coordinate amongst all given loci.

=head1 DESCRIPTION

This script takes in a list of loci and an att file, and creates an instruction file to send to cutFasta on an also provided genome_fasta file.  (See cutFasta --help if you'd like to learn more about that)

If B<--output_file> is supplied, the resulting fasta will be written to that file, otherwise it will dump to STDOUT.

If B<--region> is used, instead of providing one sequence per locus, the lowest coordinate amongst all the locus coordinates will be used for the start of a single region returned, ending at the highest coordinate found amongst all of the locus coordinates.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut


use Getopt::Long  qw(:config no_ignore_case no_auto_abbrev);
use Path::Tiny;
use Pod::Usage;

my $CUTFASTA_EXEC = '/usr/local/common/cutFasta';
my $instructions = './instructions';

my %opts;
GetOptions( \%opts,
    'att_file|a=s',
    'genome_fasta|f=s',
    'loci_list|l=s',
    'output_file|o=s',
    'region|r',
    'help|h',
) || die "Problem getting options: $!\n";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

check_options();

# Pre-load the att_file_hash
my %att;
my $afh = path( $opts{ att_file } )->filehandle( '<' );
while (<$afh>) {

    chomp;
    my ( $molecule_id, $locus, $end5, $end3, $description, $genome ) = split( '\t', $_ );
    @{$att{ $locus }}{ qw/molecule_id end5 end3/ } = ($molecule_id, $end5, $end3);

#use Data::Dumper; print Dumper (\%att); exit;
}

# get the coords for the loci from the gene attribute file.
# add a line to the instructions file.
my $lfh = path( $opts{ loci_list } )->filehandle( '<' );
my $ifh = path( $instructions )->filehandle( '>' );

if ( $opts{ region } ) {

    my ( $min, $max ) = ( 1000000000, 0 );
    my ( $first, $last ) = ( 'NA', 'NA' );
    my $molecule = undef;

    while ( <$lfh> ) {

        chomp;
        my $locus = $_;
        my ( $low, $high ) = sort { $a <=> $b } ( $att{$locus}->{end5}, $att{$locus}->{end3} );

        $min = (sort { $a <=> $b } ( $low, $min ))[0];
        $max = (sort { $a <=> $b } ( $high, $max ))[1];

        $first = $locus if ( $min == $low );
        $last  = $locus if ( $max == $high );

        if ( $first eq $locus || $last eq $locus ) {

            $molecule = $molecule // $att{$locus}->{molecule_id};
            die "Loci must all be on same molecule with --region\n" unless $molecule eq $att{$locus}->{molecule_id};

        }

    }
    print $ifh "$first..$last\t$min\t$max\t$molecule\n";

} else {

    while (<$lfh>) {

        chomp;
        my $locus = $_;
        my $line = "$locus\t$att{$locus}->{end5}\t$att{$locus}->{end3}\t$att{$locus}->{molecule_id}\n";
        print $ifh $line;

    }

}

# Feed the instruction file to cutFasta and we're done!
my @cmd = ( $CUTFASTA_EXEC, '-i', $instructions );
push( @cmd, '-o', $opts{output_file} ) if $opts{ output_file };
push( @cmd, $opts{genome_fasta} );
system( @cmd ) && die "Problem running cutFasta:\n\t" . join( ' ', @cmd ) . "\n";

exit(0);


sub check_options {

    my $errors = '';

    if ( $opts{ att_file } ) {
        $errors .= "--att_file $opts{ att_file } does not exist or is empty\n" 
            unless ( -s $opts{ att_file } );
    } else {
        $errors .= "Please provide a --att_file\n";
    }

    if ( $opts{ genome_fasta } ) { 
        $errors .= "--genome_fasta $opts{ genome_fasta } does not exist or is empty\n"
            unless ( -s $opts{ genome_fasta } );
    } else {
        $errors .= "Please provide a --genome_fasta file\n";
    }

    if ( $opts{ loci_list } ) {
        $errors .= "--loci_list $opts{ loci_list } does not exist or is empty\n" 
            unless ( -s $opts{ loci_list } );
    } else {
        $errors .= "Please provide a --loci_list\n";
    }

    die $errors if $errors;

}
