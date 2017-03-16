#!/usr/bin/env perl
use strict;
use warnings;
$|++;

=head1 NAME

align_clusters.pl - create a multiple alignment file of the members in each cluster

=head1 SYNOPSIS

    USAGE: ./align_clusters.pl --matchtable <matchtable.txt> --input_fasta <combined.fasta>

=head1 OPTIONS

B<--matchtable, -m>     :   Path to the matchtable file that associates cluster members

B<--combined_fasta, -c> :   Path to the multifasta of containing all member sequences

B<--singleton_cluster_list, -s> :   Path to list of singleton clusters to avoid

B<--project_code, -P>   :   SGE project code to use with qsub for running muscle on grid

-- OR --

B<--no_grid>            :   Do NOT invoke muscle on an SGE grid, but run it locally

B<--working_dir, -w>    :   Set working directory.  [DEFAULT: current directory]

B<--log_file, -l>       :   Path to log. [DEFAULT: ./align_clusters.log]

B<--help, -h>           :   Show this help.

=head1 DESCRIPTION

For each cluster in matchtable, create a multifasta of the members and align them using muscle.

The multifasta are created using the script 'panoct_multifasta.pl'.  If B<--singleton_cluster_list> is used to provide a file with the cluster numbers of singletons in the first column (tab separated), those numbers will be excluded from the list of cluster that have multifasta created.  Multifasta are created in the directory ./cluster_multifasta and named <cluster_number> with no extension.

After multifasta are created, muscle is invoked to create alignment files.  The invocation used is:

    muscle -in <cluster multifasta> -out <cluster number>.afa -diags -quiet -verbose

Alignment files are created in ./cluster_alignments.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Capture::Tiny qw( capture capture_merged );
use Cwd;
use FindBin;
use lib "$FindBin::Bin/../lib";
use File::Basename;
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use grid_tasks;
use IO::File;
use Pod::Usage;

my $MULTIFASTA_EXEC = "$FindBin::Bin/panoct_multifasta.pl";
my $MUSCLE_EXEC     = "$FindBin::Bin/muscle";

my $matchtable;
my $combined_fasta;
my $singleton_clusters;
my $project_code;
my $working_dir = getcwd();
my $debug = 0;

my %opts;
GetOptions( \%opts,
            'matchtable|m=s',
            'combined_fasta|c=s',
            'singleton_clusters|s=s',
            'project_code|P=s',
            'no_grid',
            'working_dir|w=s',
            'log_file|l=s',
            'log_level=i',
            'help|h',
            ) || die( "Failure getting options!\n" );
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};
check_params();

# set up log
my $log_file = $opts{ log_file } // "$working_dir/align_clusters.log";
open( my $lfh, '>', $log_file ) || die "Can't open log_file $log_file: $!\n";
my $log_dir = dirname( $log_file );

## For every cluster, create a multifasta of the cluster, then run muscle on it. ##

# First, get a listing of the cluster ids:
my $cluster_id_file = make_cluster_id_file( $matchtable, $singleton_clusters );

# create dir full of multifastas
my $multifasta_dir = create_multifastas( $cluster_id_file );

# create alignments for each cluster now
create_alignments( $cluster_id_file, $multifasta_dir );

_log( "Finished!" );

exit(0);


sub create_alignments {

    my ( $cluster_id_file, $multifasta_dir ) = @_;

    my $alignment_dir = "$working_dir/cluster_alignments";
    mkdir $alignment_dir unless ( -d $alignment_dir );

    if ( $project_code ) {

        _log( "Running muscle on grid with project_code $project_code", 0 );

        my $sh_file = write_muscle_shell_script( $multifasta_dir, $alignment_dir );

        # get the highest cluster id to use as our highest task number.
        # (dev null hack to prevent 'broken pipe' pseudo-error from 'sort' here.)
        my $max_id = `sort -nr $cluster_id_file 2>/dev/null | head -n 1`;
        chomp $max_id;

        my @grid_jobs;
        push( @grid_jobs, launch_grid_job( $project_code, $working_dir, $sh_file, 'muscle.stdout', 'muscle.stderr', '', $max_id ) );

        wait_for_grid_jobs_arrays( \@grid_jobs, 1, $max_id ) if ( scalar @grid_jobs );

    } elsif ( $opts{ no_grid } ) {

        _log( "Running muscle locally", 0 );

        open( my $cfh, '<', $cluster_id_file ) || _die( "Can't open cluster id file: $!", __LINE__ );

        while ( <$cfh> ) {

            chomp;
            my @cmd = ( $MUSCLE_EXEC, '-in', "$multifasta_dir/$_", '-out', "$alignment_dir/$_.afa", '-diags', '-quiet', '-verbose' );

            system( @cmd ) && _die( "Error running muscle on $multifasta_dir/$_", __LINE__ );

        }

    }

    _log( "muscle finished.", 0 );

}


sub write_muscle_shell_script {

    my ( $multifasta_dir, $alignment_dir ) = @_;

    my $script_name = "$working_dir/grid_muscle.sh";
    _log( "Writing muscle shell script at: $script_name", 0 );

    open( my $sfh, '>', $script_name ) || _die( "Can't write grid script $script_name: $!", __LINE__ );

    print $sfh <<"EOF";
#!/bin/sh

infile="$multifasta_dir/\$SGE_TASK_ID"
outfile="$alignment_dir/\${SGE_TASK_ID}.afa"
if [ -s \$infile ]
then
    $MUSCLE_EXEC -in \$infile -out \$outfile -diags -quiet -verbose
else
    echo "Skipping cluster \$SGE_TASK_ID"
fi
EOF

    chmod 0755, $script_name;

    return $script_name;

}


sub create_multifastas {

    my ( $cluster_id_file ) = @_;

    _log( "Running panoct_multifasta.pl to create multifastas", 0 );

    my $multifasta_dir = "$working_dir/cluster_multifastas";
    mkdir $multifasta_dir unless ( -d $multifasta_dir );

    my @cmd = ( $MULTIFASTA_EXEC, '-b', $multifasta_dir, '-P', $opts{ combined_fasta },
                '-M', $matchtable, '-C', $cluster_id_file );

    my $lf = "$log_dir/panoct_multifastas.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Can't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) && _die( "Couldn't run command: " . join( ' ', @cmd ), __LINE__ );

    } stdout => $lh;

    return $multifasta_dir;

}


sub make_cluster_id_file {

    my ( $matchtable, $singleton_clusters ) = @_;
    my $cluster_id_file = "$working_dir/cluster_ids.txt";

    _log( "Making cluster id file: $cluster_id_file", 0 );

    my %singleton_ids;

    # setup black-listed (singleton) cluster list:
    if ( $singleton_clusters ) {

        open( my $sfh, '<', $singleton_clusters ) || _die( "Can't open singleton list $singleton_clusters", __LINE__ );

        while ( <$sfh> ) {
            next if ( /^cluster/ ); # skip first line
            chomp;
            if ( /^(\d+)\t/ ) {
                $singleton_ids{ $1 }++;
            } else {
                _die( "unexpected line in singletons file: $_", __LINE__ );
            }
        }

    }

    open( my $mfh, '<', $matchtable ) || _die( "Can't open $matchtable: $!", __LINE__ );
    open( my $cfh, '>', $cluster_id_file ) || _die( "Can't open $cluster_id_file: $!", __LINE__ );
    while( <$mfh> ) {
        next if ( /^cluster id/ ); #skip first line
        if ( /^(\d+)\t/ ) {
            unless ( exists $singleton_ids{ $1 } ) {
                print $cfh "$1\n";
            }
        } else {
            _die( "unexpected line in matchtable $matchtable: $_", __LINE__ );
        }
    } 
    return $cluster_id_file;

}


sub check_params {

    my $errors = '';

    if ( $opts{ matchtable } ) {
        $matchtable = $opts{ matchtable };
        $errors .= "matchtable $matchtable empty or nonexistant\n"
            unless ( -s $matchtable );
    } else {
        $errors .= "Please provide a --matchtable.\n";
    }

    if ( $opts{ combined_fasta } ) {
        $combined_fasta = $opts{ combined_fasta };
        $errors .= "combined fasta $combined_fasta empty or nonexistant\n"
            unless ( -s $combined_fasta );
    } else {
        $errors .= "Please provide a --combined_fasta.\n";
    }

    if ( $opts{ singleton_clusters } ) {
        $singleton_clusters = $opts{ singleton_clusters };
        $errors .= "Singlton clusters file $singleton_clusters empty or nonexistant\n"
            unless ( -s $singleton_clusters );
    }

    if ( $opts{ project_code } ) {
        $project_code = $opts{ project_code };
    }
    $errors .= "Please use either --project_code or --no_grid\n"
        unless ( $project_code || $opts{ no_grid } ); 

    $working_dir = $opts{ working_dir } // $working_dir;

    $debug = $opts{ log_level } // $debug;

    die $errors if $errors;

}

sub _log {

    my ( $msg, $lvl ) = ( @_ );
    $lvl = $lvl // $debug;

    chomp( $msg );
    print $lfh "$msg\n" if ( $lvl <= $debug );

}

sub _die {

    my ( $msg, $line ) = @_;

    _log( $msg, 0 );
    die "$msg at line $line\n";

}
