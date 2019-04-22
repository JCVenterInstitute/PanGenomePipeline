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

use strict;
use warnings;
$|++;

=head1 NAME

    compute_pangenome_ani.pl - compute weighted ani between pairs of genomes in a pangenome based on cluster membership.

=head1 SYNOPSIS

    compute_pangenome_ani.pl [ --working_dir <pangenome working directory> ]
                             [ --genome_list <genomes.list> ]
                             [ --combined_att <combined.att> ]
                             [ --cluster_alignment_dir <directory of .afa files> ]
                             [ --output_file <output file> ]
                             [ --help ]

=head1 OPTIONS

B<--working_dir, -w>    :   Path to a JCVI Pan-genome Pipeline working directory.  Default is '.' (Current directory)

B<--genome_list, -g>    :   Path to a list of genomes.  Default is B<--working_dir>/genomes.list

B<--combined_att, -a>   :   Path to a combined gene attribute file.  Default is B<--working_dir>/combined.att

B<--cluster_alignment_dir, -c>  :   Path to a directory containing ortholog cluster alignments in .afa format.  Default: B<--working_dir>/results/cluster_alignments

B<--output_file, -o>    :   Path to desired output file.  Default: B<--working_dir>/results/ani_by_clusters.pl

B<--help, -h>           :   Display this help

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Bio::Seq;
use Bio::SeqIO;
use Cwd qw( getcwd abs_path );
use Getopt::Long qw( :config no_ignore_case no_auto_abbrev );
use Pod::Usage;

use Data::Dumper;

my $working_dir;
my $combined_att;
my $genome_list;
my $clust_align_dir;
my $output_file;

my %opts;
GetOptions( \%opts,
            'working_dir|w=s',
            'combined_att|a=s',
            'genome_list|g=s',
            'cluster_alignments_dir|c=s',
            'output_file|o=s',
            'help|h',
            ) || die "Problem getting options: $!\n";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{ help };

check_params();

# load genome list
print "Loading Genome List: $genome_list\n";
my @genomes;
open( my $glf, '<', $genome_list ) || die "Can't open $genome_list: $!\n";
while ( <$glf> ) {
    chomp;
    push @genomes, $_;
}


# Load gene->genome mapping from combined.att
print "Loading locus->genome mapping from $combined_att\n";
my %locus2genome;
open( my $alf, '<', $combined_att ) || die "Can't open $combined_att: $!\n";
while ( <$alf> ) {

    chomp;
    my ( $locus, $genome ) = (split(/\t/,$_))[1,5];
    $locus2genome{ $locus } = $genome;

}


# Begin opening alignments.
print "Parsing alignments from $clust_align_dir\n";
my %values;
opendir( my $cad, $clust_align_dir ) || die "Can't opendir $clust_align_dir: $!\n";
while ( defined( my $afa_file = readdir($cad) ) ) {

    next if ( $afa_file =~ /^.{1,2}$/ );
    my $afa_path = "$clust_align_dir/$afa_file";

    print "\tParsing $afa_path\n";

    my @seqs;
    my $in_io = Bio::SeqIO->new( -file => "<$afa_path", -format => 'fasta' );
    while ( my $seqin = $in_io->next_seq ) {
        push @seqs, $seqin;
    }

    while ( @seqs > 1 ) {
        # shift off first sequences.
        my $seqOne = shift @seqs;
        my $seqOneSeq = $seqOne->seq();
        my $seqOneID  = $seqOne->primary_id();
        my $seqOneGenome = $locus2genome{ $seqOneID };

        # Find identity between pairs.
        for my $seqTwo ( @seqs ) {
            my $seqTwoSeq       = $seqTwo->seq();
            my $seqTwoID        = $seqTwo->primary_id();
            my $seqTwoGenome    = $locus2genome{ $seqTwoID };

            # Store Identity and length of alignment in both parts of the genome hash.
            my ( $ali_len, $ali_identity ) = compute_identity( $seqOneSeq, $seqTwoSeq );

            $values{ $seqOneGenome }->{$seqTwoGenome}->{ ali_len } += $ali_len;
            $values{ $seqOneGenome }->{$seqTwoGenome}->{ ali_identity } += $ali_identity;
            $values{ $seqTwoGenome }->{$seqOneGenome}->{ ali_len } += $ali_len;
            $values{ $seqTwoGenome }->{$seqOneGenome}->{ ali_identity } += $ali_identity;

        }
    }

}

# Uncomment to see the values hash.
#print Dumper (\%values); exit;

print "Printing results to $output_file\n";
open( my $ofh, '>', $output_file ) || die "Can't write to: $output_file: $!\n"; 

# Print header:
my $header = join( "\t", ('',@genomes[1 .. $#genomes]));
print $ofh "$header\n";
my $orig_size = $#genomes;
# Per genome in genome list (While genomes > 1)
while ( @genomes > 1 ) {

    # First genome starts a pair
    my $genomeOne = shift @genomes;
    my $padding = $orig_size - scalar(@genomes);
    my @line = ($genomeOne);
    while ( $padding-- > 0 ) {
        push @line, '';
    }

    # Next genome forms a pair
    for my $genomeTwo ( @genomes ) {
        # From the first genome, find weighted ANI for the pair
        my $weighted_ani = sprintf( "%.2f", $values{ $genomeOne }->{ $genomeTwo }->{ ali_identity } / $values{ $genomeOne }->{ $genomeTwo }->{ ali_len } * 100 );
        push @line, $weighted_ani;
    }

    print $ofh join("\t",@line),"\n";

}

print "Finished!\n";
exit(0);


sub compute_identity {
#   Given two seqs from an alignment (afa) file, compute the alignment length
#   and identity between the seqs.
#   For now, skip ALL gaps, whether they be on both seqs or just one.

    my ( $s1, $s2 ) = @_;

    my $ali_len     = 0;
    my $identity    = 0;

    for my $i ( 0 .. (length($s1) - 1) ) {

        my $s1c = substr($s1, $i, 1);
        my $s2c = substr($s2, $i, 1);
        #next if ( ($s1c eq '-') && ($s2c eq '-') ); # Only skips gaps from BOTH seqs
        next if ( ($s1c eq '-') || ($s2c eq '-') );  # Skips gaps from EITHER seq
        $ali_len++;
        $identity++ if ( $s1c eq $s2c );

    }

    return( $ali_len, $identity );

}


sub check_params {

    my $errors = '';

    $working_dir = $opts{ working_dir } // getcwd();
    $errors .= "Can't find --working_dir $working_dir\n" unless ( -d $working_dir );
    chop $working_dir if $working_dir =~ m(/$); # This *MIGHT* break symlinks to a working dir

    $combined_att = $opts{ combined_att } // "$working_dir/combined.att";
    if ( ! -f $combined_att ) {
        $errors .=  "Can't find --combined_att $combined_att\n";
    } elsif ( ! -s $combined_att ) {
        $errors .= "--combined_att $combined_att appears to be empty.\n";
    }

    $genome_list = $opts{ genome_list } // "$working_dir/genomes.list";
    if ( ! -f $genome_list ) {
        $errors .=  "Can't find --genome_list $genome_list\n";
    } elsif ( ! -s $genome_list ) {
        $errors .= "--genome_list $genome_list appears to be empty.\n";
    }

    $output_file = $opts{ output_file } // "$working_dir/results/ani_by_clusters.txt";
    
    $clust_align_dir = $opts{ cluster_alignment_dir } // "$working_dir/results/cluster_alignments";
    $errors .= "Can't find --cluster_alignment_dir $clust_align_dir\n" unless ( -d $clust_align_dir );

    die $errors if $errors;

}
