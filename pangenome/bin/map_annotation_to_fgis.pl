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

map_annotation_to_fgis.pl - map annotations from search results to fgis

=head1 SYNOPSIS

    USAGE: map_annotation_to_fgis.pl -i <fgi report index> -c <centroid cluster roles> -o <output file>

=head1 OPTION

=over

B<--fgi_index, -i>          :   path to fgi report index, default: ./results/fGIs/fGI_report.txt.index

B<--annotation_files, -c>   :   path(s) to centroid annotation file(s), default: ./results/centroids.cluster_roles.txt and
                                ./results/plasmid_blast.out

B<--output_file, -o>        :   path to output file, default: ./results/fGI_annotation.txt

=back

=head1 DESCRIPTION

Given a cluster roles file and an fgi index file, determine the annotations that belong to each fGI.

=head1 INPUT

B<--fgi_index> expects a file containing two columns: cluster number and fgi number.  This file is generated
by gene_order.pl normally as part of the pangenome pipeline.

B<--annotation_files> expects a file or comma seperated list of files containing three columns: centroid_id, role description, and Gene Ontology (GO) ID.
Note that GO:ID may be blank, especially for clusters that have "Other" as the description.  It should also be noted
that centroid id is expected to be in the format "centroid_<cluster number>" where "cluster number" corresponds to a
value seen in B<--fgi_index> 

=head1 OUTPUT

B<--output_file> will contain three columns: fgi_number, role description, and Gene Onotology (GO) ID.  
It should be pointed out that the role description and GO ID columns are switched relative to how they appear in the 
B<--annotation_files> and the first column (centroid_id) is replaced with the corresponding fgi ids from B<--fgi_index>
and sorted by the same. 

=head1 CONTACT

    Jason Inman
    jinmanAjcvi.org

=cut

use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use Path::Tiny;
use Pod::Usage;

use Data::Dumper;

my $DEFAULT_ANNOTATION_FILES = './results/centroids.cluster_roles.txt,./results/plasmid_blast.out';
my $DEFAULT_FGI_INDEX        = './results/fGIs/fGI_report.txt.index';

my @annot_files;

my %opts;
GetOptions( \%opts,
            'fgi_index|i=s',
            'annotation_files|c=s',
            'output_file|o=s',
            'help|h',
        ) || die "Can't get options: $!\n";

pod2usage( { -exirtval => 0, -verbose => 2 } ) if $opts{ help };

check_options();

# Load the mapping:
my %seqs;
my %fgis;
my $mfh = path( $opts{ fgi_index } )->filehandle( '<' );
while( <$mfh> ) {

    chomp;
    my ( $seq, $fgi ) = split("\t", $_);

#    $seqs{ $seq }->{ fgi } = $fgi;
    $seqs{ $seq }->{ annotation } = [];
    push( @{$fgis{ $fgi }}, $seq );

}

# Read in the annotations
for my $annot_file ( @annot_files ) {

    my $afh = path( $annot_file )->filehandle( '<' );
    while( <$afh> ) {
        chomp;
        my $line = $_; # Only doing this so we can make a nice warning in a bit.
        my @line = split("\t",$_);
        my $centroid = shift @line;
        my $cluster_id;
        if ( $centroid =~ /centroid_(\d+)/ ) {
            $cluster_id = $1;
        } else {
            warn( "Can't parse centroid id from centroid name '$centroid' in the following line:\n$line\n" ); # see.  Toldja.
            next;
        }
#        print "LINE: $line\n\@LINE: ", join("TAB",@line),"\n";
        push( @{$seqs{ $cluster_id }->{ annotation }}, join( "\t", @line ) );
    }

}


# Print the output
my $ofh = path( $opts{ output_file } )->filehandle( '>' );

for my $fgi ( sort { $a <=> $b } keys %fgis ) {
    my %annotations;
    for my $cluster ( @{$fgis{ $fgi }} ) {
        for my $annotation ( @{$seqs{ $cluster }->{annotation}} ) {
            $annotations{ $annotation }++;
        }
    }

    # This might seem redundant.  But sorting on the keystring ends up sorting on both the description and GO ID at the same
    # time, and then splitting lets us print the GO terms first, which was the desired output.
    for my $keystring ( sort { $a cmp $b } keys %annotations ) {

        my ( $description, $GO ) = split("\t",$keystring);
        $GO = $GO // '';
        #print "$fgi\t$GO\t$annotations{ $keystring }\t$description\n";
        print $ofh "$fgi\t$GO\t$annotations{ $keystring }\t$description\n";

    }

}

exit(0);


sub check_options {

    my $errors = '';

    my $annot_file_string = $opts{ annotation_files } // $DEFAULT_ANNOTATION_FILES;
    # Should check if all user-supplied annotation files are present.  However, if
    # relying on the defaults, for example as part of the pangenome pipeline, it is
    # quite possible there were no hits in some of the annotation files.
    for my $annot_file ( split( ',', $annot_file_string ) ) {
        if ( -s $annot_file ) {
            push( @annot_files, $annot_file );
        } else {
            $errors .= "Can't find annotation file: $annot_file\n" if ( $opts{ annotation_files } );
        }
    }

    $opts{ fgi_index } = $opts{ fgi_index } // $DEFAULT_FGI_INDEX;
    $errors .= "Please provide an --fgi_index\n" unless ( -s $opts{ fgi_index } );

    die $errors if $errors;

}
