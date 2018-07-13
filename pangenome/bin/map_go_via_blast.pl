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

map_go_via_blast.pl - map GO terms not covered by HMMs to centroids.

=head1 SYNOPSIS

    USAGE: map_go_via_blast.pl -P <project code> | --blast_local | --blast_file <blast file name>
                                --input_seqs <input fasta file>
                                --output_file <output file path>
                                [ --blast_db <path to db> ] [ --evalue <evalue> ] 
                                [ --percent_id <percent identity> ] 
                                [ --percent_coverage <percent coverage> ]

=head1 OPTIONS

Choose ONE of the following THREE options:

=over 

B<--project, -P>        :   SGE/UGE grid accounting project code

B<--blast_local>        :   Do NOT run BLAST in parallel on a computing grid, but run it locally.

B<--blast_file, -b>     :   Skip BLAST, run using this results file instead.

=back

B<--blast_db, -s>       :   Path to BLAST db for searching.

B<--mapping_file, -m>   :   Path to file mapping blast headers to go, roles, etc.

B<--output_file, -o>    :   Path to output file   

B<--evalue, -E>         :   e-value required to pass and be reported as a hit. [DEFAULT: 10e-5 ]

B<--percent_id, -I>     :   percent identity required to be reported as a hit. [DEFAULT: 35 ]

B<--percent_cov, -C>    :   min percent coverage required to be reported as a hit. [DEFAULT: 80 ]

B<--input_seqs, -i>     :   Input fasta file.

B<--use_nuc, -n>        :   Input is in nucleotide space, also use blastn.

=head1 DESCRIPTION

This script takes in an input fasta file, and a predetermined blast database with associated GO mapping, and maps the GO terms to the input sequence after examining the blast results.

=head1 INPUT

B<--input_seqs>   - typically centroids.fasta

B<--blast_db>     - blast db of some group of interesting genes.

B<--mapping_file> - maps the sequences in the B<--blast_db> to associated GO terms.

=head1 OUTPUT

B<--output_file> - contains the mappings in the format:

<input_seq_id>\t<hit description>\t<GO Term>\n

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Data::Dumper;

use Capture::Tiny qw{ capture_merged };
use Cwd;
use File::Path qw( mkpath remove_tree );
use FindBin;
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use IO::File;
use Pod::Usage;

use lib "$FindBin::Bin/../lib";
use grid_tasks;

my $BIN_DIR     = $FindBin::Bin;
my $BLASTN_EXEC = '/usr/local/packages/ncbi-blast+/bin/blastn';
my $BLASTP_EXEC = '/usr/local/packages/ncbi-blast+/bin/blastp';
my $BLASTDB_DIR = "$BIN_DIR/BLAST2GO";
my $BLASTN_DB   = "$BLASTDB_DIR/plasmid_finder_rep.seq";
my $BLASTP_DB   = "$BLASTDB_DIR/plasmid_finder_rep.pep";
my $BLAST_MAP   = "$BLASTDB_DIR/plasmidDB2GO_lookup.txt";
my $SPLIT_FASTA = "$BIN_DIR/split_fasta.pl";

my $DEFAULT_EVALUE      = '10e-5';
my $DEFAULT_PERCENT_ID  = '35';
my $DEFAULT_PERCENT_COV = '80';

my %opts;
GetOptions( \%opts,
            'project|P=s',
            'blast_local',
            'blast_file|b=s',
            'input_seqs|i=s',
            'blast_db|s=s',
            'mapping_file|m=s',
            'evalue|E=s',
            'output_file|o=s',
            'percent_id|I=i',
            'percent_cov|C=i',
            'use_nuc|n',
            'working_dir|w=s',
            'log_dir=s',
            'help|h',
        ) || die "Problem getting options.\n";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{ help };

check_options();

# Run blast on the input, if required.
my $blast_file;
if ( $opts{ blast_file } ) {

    $blast_file = $opts{ blast_file };

} else {

    $blast_file = "$opts{ working_dir }/blast_output";

    my $blast_prog = $opts{ use_nuc } ? $BLASTN_EXEC : $BLASTP_EXEC;
    my $blast_db = $opts{ blast_db } // $opts{ use_nuc } ? $BLASTN_DB : $BLASTP_DB;
    my $blast_cmd = "$blast_prog -db $blast_db -evalue $opts{ evalue }"; 
    $blast_cmd .=   " -qcov_hsp_perc $opts{ percent_cov }";
    $blast_cmd .=   " -perc_identity $opts{ percent_id }" if ( $opts{ use_nuc } );
    $blast_cmd .=   " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\"";

    if ( $opts{ blast_local } ) {

        # add in the query
        $blast_cmd .= " -query $opts{ input_seqs } -out $blast_file";

        #   Run the blast
        my $lf = "$opts{ log_dir }/blast.log";
        my $lh = IO::File->new( $lf, "w+" ) || die( "Couldn't open $lf for logging: $!\n" );
 
        capture_merged{

            system( $blast_cmd ) == 0 || die( "Problem running blast.  See $lf\n" );

        } stdout => $lh;

    } else {

        # Going to run blast on the grid.
        my $blast_dir = "$opts{ working_dir }/blast";
        my $split_dir = "$blast_dir/split_fastas";
        mkpath( $split_dir ) unless ( -d $split_dir );
        my $split_cmd = "$SPLIT_FASTA -f $opts{ input_seqs } -n 1000 -o $split_dir";

        unless ( system( $split_cmd ) == 0 ) {

            # something went wrong and we should check the results carefully.
            # Check for errors of the type "Expected: unique FASTA identifier" in split_fasta.log
            my $split_fasta_log = "$opts{ working_dir }/split_fasta.pl.error";
            my $found_dupe = 0;
            if ( -s $split_fasta_log ) {

                open( my $sflh, '<', $split_fasta_log ) || die "Can't open $split_fasta_log to investigate split_fasta.pl failure\n";
                while( <$sflh> ) {

                    $found_dupe++ if (/Expected: unique FASTA identifier/);

                }

            }

            my $dupe_message = ( $found_dupe ) ?
                                "It looks like duplicate locus tags are invovled.\n" :
                                "It doesn't look like duplicate locus tags are invovled.\n";

            die( "Error running split_fasta.  $dupe_message\n");

        }

        # Write shell script
        my @file_list = <$split_dir/split_fasta.*>;
        my $sh_file = write_blast_shell_script( $blast_file, $split_dir, $blast_dir, $blast_cmd );
        print "Running blast on the grid.\n";

        # Launch blast job array, wait for finish
        my @grid_jobs;
        push( @grid_jobs, launch_grid_job( $opts{ project }, $opts{ working_dir }, $sh_file, 'blast.stdout', 'blast.stderr', "", scalar @file_list ) );
        print "Waiting for blast jobs to finish.\n";
        wait_for_grid_jobs_arrays( \@grid_jobs, 1, scalar( @file_list ) ) if ( scalar @grid_jobs );
        print "Blast jobs finished!\n";

        # Cat all blast files together
        open( my $cfh, ">", $blast_file ) || die ( "Couldn't open $blast_file for writing: $!\n" );
        _cat( $cfh, glob( "$blast_dir/blast_output.*" ) );
        close $cfh; # Force the buffer to flush to the output file so it can be seen as non-empty sooner
                    # rather than later.

        if ( -e $blast_file ) {

            print "Removing intermediate blast files.\n";

            # Remove intermediate blast dir.
            remove_tree( $blast_dir, 0, 1 );

        } else {

            die( "Problem getting blast results.\n");

        }

        if ( -z $blast_file ) {

            print "No results from blast\n";
            exit(0);

        }

    }

}

# Parse blast output ( $blast_file, now ).
print "BLAST results at: $blast_file\n";

open( my $bfh, '<', $blast_file ) || die "Can't open $blast_file for reading: $!\n";
open( my $mfh, '<', $opts{ mapping_file } ) || die "Can't open mapping file $opts{ mapping_file }: $!\n";

my %plasmid2go;

while (<$mfh>) {

    chomp;
    my ( $id, $line );
    if ( /^([^\t]+)\t(.*)$/ ) {
        ( $id, $line ) = ( $1, $2 );
        $id =~ tr/ //; # 
        $plasmid2go{ $id } = $line;
    } else {
        warn "There was a problem parsing this line in mapping_file $opts{ mapping_file }:\n$_\n";
    }

}

my @output_lines;

while (<$bfh>) {

    chomp;
    my ( $input_seq, $hit_description, $perc_id, $evalue, $coverage ) = (split(/\t/,$_))[0,1,2,10,12];

    # For debugging this will show how the blast is getting parsed.
#    print "$input_seq, $hit_description, $perc_id, $evalue, $coverage\n";

    $hit_description =~ tr/ //;
    next unless ( $evalue <= $opts{ evalue } );
    push @output_lines, "$input_seq\t$hit_description\t$plasmid2go{ $hit_description }\n";    

}

open( my $ofh, '>', $opts{ output_file } ) || die "Can't open $opts{ output_file } for writing: $!\n";
print $ofh map { $_ } @output_lines;

exit(0);


sub write_blast_shell_script {

    my ( $fasta, $sdir, $bdir, $cmd ) = @_;

    my $cmd_string = $cmd . " -query $sdir" . '/split_fasta.$SGE_TASK_ID ' . "-out $bdir" . '/blast_output.$SGE_TASK_ID';
    my $script_name = "$bdir/grid_blast.sh";

    open( my $gsh, '>', $script_name ) || _die( "Can't open $script_name: $!\n", __LINE__ );

    print $gsh "#!/bin/tcsh\n\n";
    print $gsh "$cmd_string\n";

    chmod 0755, $script_name;

    return $script_name;

}


sub _cat {
# Given a list of file names, concatonate the first through n minus one-th
# onto the nth.
    
    my ( $output_fh, @input ) = ( @_ );
    
    for ( @input ) {
        if(-s $_){
            open ( my $ifh, '<', $_ );
            while ( <$ifh> ) { print $output_fh $_ };
        }
    }
}


sub check_options {

    my $errors = '';

    # Only want one of ( --project --blast_local --blast_file )
    if ( ( $opts{ project } && ( $opts{ blast_local } || $opts{ blast_file } ) ) ||
         ( $opts{ blast_local } && $opts{ blast_file } ) ) {
        $errors .= "Please specifiy only one of --project, --blast_local, or --blast_file\n";
    }

    # But must have one of them.
    unless ( $opts{ project } || $opts{ blast_local } || $opts{ blast_file } ) {
        $errors .= "Please specify one of --project, --blast_local, or --blast_file\n";
    }

    # If blast_file, make sure it exists:
    if ( $opts{ blast_file } ) {
        unless ( -s $opts{ blast_file } ) {
            $errors .= "--blast_file $opts{ blast_file } has no size or doesn't exists.\n";
        }
    }

    # Gotta have an input file or this is pointless.
    if ( $opts{ input_seqs } ) {
        unless( -s $opts{ input_seqs } ) {
            $errors .= "input seq file $opts{ input_seqs } is empty or non-existant.\n";
        }
    } else {
        $errors .= "--input_seqs is necessary.\n";
    }

    # Also need a file mapping blast headers to the data:
    $opts{ mapping_file } = $opts{ mapping_file } // $BLAST_MAP;
    unless ( -s $opts{ mapping_file } ) {
        $errors .= "mapping file $opts{ mapping_file } is empty -r non-existant.\n";
    }

    unless ( $opts{ output_file } ) {
        $errors .= "--output_file is required.\n";
    }


    # These params have defaults:
    $opts{ evalue }      = $opts{ evalue }      // $DEFAULT_EVALUE;
    $opts{ percent_id }  = $opts{ percent_id }  // $DEFAULT_PERCENT_ID;
    $opts{ percent_cov } = $opts{ percent_cov } // $DEFAULT_PERCENT_COV;

    $opts{ working_dir } = $opts{ working_dir } // getcwd();
    $opts{ log_dir }     = $opts{ log_dir }     // getcwd();

    die $errors if $errors;

}
