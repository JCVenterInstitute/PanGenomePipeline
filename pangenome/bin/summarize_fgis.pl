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

summarize_fgis.pl - create files presenting genome->fgi relationships on a per genome and per fgi basis.

=head1 SYNOPSIS

    USAGE: summarize_fgis.pl -p <pangenome dir> 

        - OR -

    summarize_fgis.pl -g [genomes.list] -i [fGI_report.txt.index] -a [combined.att] -m [matchtable.txt] -o <output dir> 


=head1 OPTIONS

=over 

=item B<--pangenome_dir, -p>    :   Path to a Pan-genome pipeline working directory.

=back

 - Alternatively, provide each of the following:

=over

=item B<--genome_list, -g>      :   Path to genomes.list

=item B<--fgi_index, -f>        :   Path to fGI_report.txt.index

=item B<--matchtable, -m>       :   Path to matchtable.txt

=item B<--att_file,-a>          :   Path to combined.att

=item B<--output_dir, -o>       :   Path to directory where output will be written.

=back

 - Miscellaneous options:

=over

=item B<--help>                 :   Display this help.

=back

=head1 DESCRIPTION

When given certain inputs from a Pan-genome Pipeline run, this script produces two files
that summarize fgi to genome relationships in two ways:

    1. fGI 'weight', or percentage of an fGI's clusters for which a genome has a member.
    2. a list of which fGIs are found in each genome.


=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Getopt::Long qw( :config no_ignore_case no_auto_abbrev );
use Pod::Usage;
use Cwd;

my $pangenome_dir;
my $genome_list;
my $fgi_index;
my $matchtable;
my $att_file;
my $output_dir;

my %opts;
GetOptions( \%opts,
            'pangenome_dir|p=s',
            'genome_list|g=s',
            'fgi_index|f=s',
            'matchtable|m=s',
            'att_file|a=s',
            'output_dir|o=s',
            'help|h',
            ) || die "Can't read options! $!\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{help};

check_params();

# Get genome order for printing
my @genome_order;
open( my $gfh, '<', $genome_list ) || die "Can't open genome list $genome_list: $!\n";
while ( <$gfh> ) {

    chomp;
    push @genome_order, $_;

}

# make gene 2 genome hash from combined.att
my %gene2genome;
open( my $afh, '<', $att_file ) || die "Can't open att_file $att_file: $!\n";
while ( <$afh> ) {

    chomp;
    my ( $gene, $genome ) = (split(/\t/,$_))[1,5];
    $gene2genome{ $gene } = $genome;

}

# make matchtable (genome2cluster) hash
## (hash containing arrays of the clusters per genome)
my %genome2cluster;
open( my $mfh, '<', $matchtable ) || die "Can't open matchtable $matchtable: $!\n";
while ( <$mfh> ) {

    chomp;
    my @line = split( /\t/, $_ );
    my $cluster = shift @line;
    for my $gene ( @line ) {

        next if $gene =~ /---/;

        push @{$genome2cluster{ $gene2genome{ $gene } }}, $cluster;

    }

}

# make fgi (fgis2clusters) hashes
my %fgi2cluster;
my %cluster2fgi;
open( my $fhf, '<', $fgi_index ) || die "Can't open fgi_index $fgi_index: $!\n";
while ( <$fhf> ) {

    chomp;
    my ( $cluster, $fgi ) = split( /\t/, $_ );
    push @{$fgi2cluster{ $fgi }}, $cluster;
    $cluster2fgi{ $cluster } = $fgi;

}

# prepare output files
my $fgi_weight = "$output_dir/fgi_weights";
open( my $wofh, '>', $fgi_weight ) || die "Can't open $fgi_weight for writing: $!\n";
print $wofh "fgi\t" . join("\t", @genome_order ) . "\n";

my $genome2fgi = "$output_dir/genome2fgi";
open( my $gofh, '>', $genome2fgi ) || die "Can't open $genome2fgi for writing: $!\n";


# Create weights file
for my $fgi ( keys %fgi2cluster ) {

    print $wofh "$fgi";
    for my $genome ( @genome_order ) {

        my @genome_clusters = @{$genome2cluster{ $genome }};

        my $fgi_cluster_sum = scalar @{$fgi2cluster{ $fgi }};
        my $genome_cluster_sum = 0;

        for my $cluster ( @{$fgi2cluster{ $fgi }} ) {

            for my $gc ( @genome_clusters ) {

                if ( $gc == $cluster ) {

                    $genome_cluster_sum++;
                    next;

                }

            }

        }

        my $weight = $genome_cluster_sum ? sprintf( "%.2f", ( $genome_cluster_sum / $fgi_cluster_sum * 100 ) ) : 0;
        print $wofh "\t" . $weight;

    }

    print $wofh "\n";

}

# Create file listing fgis per genome
for my $genome ( @genome_order ) {

    my %seen_fgis;
    for my $cluster ( @{$genome2cluster{$genome}} ) {
        $seen_fgis{ $cluster2fgi{ $cluster } }++;
    }

    print $gofh "$genome\t";
    print $gofh join(",", sort { $a <=> $b } keys %seen_fgis );
    print $gofh "\n";
    

}

exit(0);


sub check_params {

    my $errors = '';

    if ( $opts{ pangenome_dir } ) {

        $errors .= "Please provide a --pangenome_dir\n" unless ( -d $opts{ pangenome_dir } );

        $pangenome_dir = $opts{ pangenome_dir };

        $fgi_index      = "$pangenome_dir/results/fGIs/fGI_report.txt.index";
        $matchtable     = "$pangenome_dir/results/matchtable.txt";
        $att_file       = "$pangenome_dir/combined.att";
        $genome_list    = "$pangenome_dir/genomes.list";
        $output_dir     = "$pangenome_dir/results/fGIs";

    } else {

        $fgi_index      = $opts{ fgi_index };
        $matchtable     = $opts{ matchtable };
        $att_file       = $opts{ att_file };
        $genome_list    = $opts{ genome_list };
        $output_dir     = $opts{ output_dir } // getcwd();

        $errors .= "Please supply --fgi_index if not using --pangenome_dir\n"
            unless ( $fgi_index );
        $errors .= "Please supply --matchtable if not using --pangenome_dir\n"
            unless ( $matchtable );
        $errors .= "Please supply --att_file if not using --pangenome_dir\n"
            unless ( $att_file );
        $errors .= "Please supply --genome_list if not using --pangenome_dir\n"
            unless ( $genome_list );

    }

    if ( $fgi_index ) {
        $errors .= "$fgi_index empty/not found\n" unless ( -s $fgi_index );
    }
    if ( $matchtable ) {
        $errors .= "$matchtable empty/not found\n" unless ( -s $matchtable );
    }
    if ( $att_file ) {
        $errors .= "$att_file empty/not found\n" unless ( -s $att_file );
    }
    if ( $genome_list ) {
        $errors .= "$genome_list empty/not found\n" unless ( -s $genome_list );
    }
    if ( $output_dir ) {
        $errors .= "$output_dir not found\n" unless ( -d $output_dir );
    }

    die $errors if $errors;

}
