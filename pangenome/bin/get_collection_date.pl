#!/usr/bin/env perl
use warnings;
use strict;
$|++;

=head1 NAME

get_collection_date.pl

=head1 SYNOPSIS

    USAGE: ./get_collection_date.pl -l <list of genbank files> [-o output_file]

=head1 OPTIONS

=over

=item B<--list_file, -l>    :   File containing paths to all genbank files to parse

=item B<--output, -o>       :   Path to output file [DEFAULT: ./collection_dates]

=item B<--help, -h>         :   Display this help.

=back

=head1 DESCRIPTION

Given a list of genbank files and, optionally, a path to an output file, this script will
attempt to retrieve each genome's collection date, if included in the source feature tags.

A tab-separated output file is created at the specified location (or, by default, at
./collection_dates) in the format:

    GenomeName<tab>CollectionDate

Some notes:

    * The name in the first column is the filename less the final 
      'extension' for each line in the --list_file. For example, 
      if the path is "/foo/bar.gb" it becomes "bar" in the output file.
    * The values in the second column come from the first source 
      feature with a /collection_date tag encountered in each file.
      Since /collection_date doesn't have a particular format, it 
      might pull things in any number of format depending on the 
      diversity of the set of input files. As of now, this script 
      makes no attempt to validate the contents of that tag or put 
      them all into the same format.
    * When no sequence in a file contains a /collection_date on any 
      source features, the genome will have 'NA' as the value in the 
      second column.

=head1 CONTACT

    Jason Inman

=cut

use Getopt::Long qw( :config no_ignore_case no_auto_abbrev );
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;
use Cwd;
use Pod::Usage;

my $DEFAULT_OUTPUT_NAME = './collection_dates';
my $DEFAULT_COLLECTION_DATE = 'NA';

my %opts;
GetOptions( \%opts,
            'list_file|l=s',
            'output|o=s',
            'help|h',
            ) || die "Can't get options!!\n$!\n";

pod2usage( { -verbose => 2, -exitval => 0 } ) if $opts{ help };
check_params();

# Open file handles.
open( my $lfh, '<', $opts{ list_file } ) || die "Can't open input list file $opts{ list_file }: $!\n";
open( my $ofh, '>', $opts{ output } )    || die "Can't open output file $opts{ output }: $!\n";

while ( my $file = <$lfh>) {

    chomp( $file );
    my ( $filename, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $collection_date = undef;

    # Create a new SeqIO object for the file:
    my $seqio_object = Bio::SeqIO->new( -file => $file, -format => 'genbank' );
    my $seq_object;

    # go through until we have found a collection date.
    while ( ( not defined $collection_date ) and ( $seq_object = $seqio_object->next_seq() ) ) {

        # look for the source feature.  This is where the /collection_date tag resides
        for my $feat_object ( $seq_object->get_SeqFeatures ) {

            if ( $feat_object->primary_tag eq 'source' ) { 
                $collection_date = $DEFAULT_COLLECTION_DATE;
                if ( $feat_object->has_tag( 'collection_date' ) ) {
                    ( $collection_date ) = $feat_object->get_tag_values( 'collection_date' );
                }
                # print what we found.
                print $ofh "$filename\t$collection_date\n";
                last;
            } else {
                next;  #shouldn't have to be here, source should be first but, hey.
            }

        }

    }

    # Print default collection date if nothing is found.
    if ( not defined $collection_date ) {
        print $ofh "$filename\t$DEFAULT_COLLECTION_DATE\n";
    }

}


exit(0);


sub check_params {

    my $errors = '';

    $opts{ output } = $opts{ output } // $DEFAULT_OUTPUT_NAME;

    $errors .= "Please supply a --list_file\n" unless $opts{ list_file };

    die $errors if $errors;

}
