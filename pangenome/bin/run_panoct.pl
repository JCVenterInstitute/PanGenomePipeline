#!/usr/local/bin/perl

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

run_panoct.pl - utility for running the pangenome pipeline

=head1 SYNOPSIS

    USAGE: run_panoct.pl  -d /path/to/genome_list_file
                            -u <sybase user> && -p <sybase password> || -D database_file
                           [
                            --use_nuc
                            --gene_att_file /path/to/gene_attribute_file
                            --blast_file /path/to/blast/output/file
                            --panoct_local
                            --blast_local
                            --project_code
                            --no_grid
                            --no_stats
                            --no_trees
                            --working_dir /path/to/working_dir
                            --log_file /path/to/log_file
                            --reference
                            --role_ids
                            --terms
                            --percent_id
                            --strict
                            --hmm_file
                            --role_lookup
                            --no_graphics
                            --lite
                          ]


=head1 OPTIONS

=head2 Required options

B<--genome_list_file,-d>            :   File containing list of genome names and download location

=head2 Optional Options

B<--att_dir, -A>                :   Directory containing genome .att files

B<--gene_att_file. -a>          :   Use this gene_att_file instead of creating one from att_dir 

B<--blast_file,-b>              :   Path to a file of all-vs-all blast results. Must be -m 8 or -m 9 tab delimited output.

B<--no_stats>                   :   Flag to not run statistics program. 

B<--panoct_local>               :   Runs PanOCT locally <JCVI Specific>

B<--blast_local>                :   Runs Blast locally <JCVI Specific>

B<--limiting_genome>            :   Genome(s) limiter used in stat generation

B<--role_ids>                   :   Role id limiter used in stat generation (Surrounded in quotes,separated by comma)

B<--terms>                      :   Protein name limiter used in stat generation (Surrounded in quotes, separated by commas)

B<--strict>                     :   Level of clustering method strictness(none, low, high) (Default: low)

B<--percent_id>                 :   Blast cutoff to use within a clustering method's algorithm

B<--methods>                    :   Options are either panoct or none (skips the clustering step). (Default: panoct)

B<--hmm_file>                   :   A file of HMMs of interest, will produce files with clusters that have these HMMs

B<--project_code, -P>           :   Valid project code for grid-accounting

B<--no_grid>                    :   Run all subtasks locally.

B<--no_trees>                   :   Don't run tree-making programs following the clustering algorithms

B<--use_nuc>                    :   use nucleotide versions of blast and input files.

B<--no_graphics>                :   Don't execute the various scripts that produce graphics

B<--lite>                       :   Don't run anything after panoct.  Intended for use by run_pangnome.pl in hierarchical mode.

=head2 Misc options

B<--log_dir>                    :   Place to write run_panoct.log, along with a small collection of log files

B<--debug>                      :   Integer specifying log message levels.  Lower is more verbose.

B<--help,-h>                    :   Display this help message.

=head1  DESCRIPTION

Runs several scripts that make up the pangenome pipeline:

=head2 Data Pulling 

<External Data: Genbank>
FTP Files can be downloaded from Genbank/RefSeq. These files will
be stored in the "fasta_dir" directory. Each genome will get it's own
sub directory and the resulting combined.fasta will be created.

Sequences are stored in a directory named 'fasta_dir' within the working dir.  This
step can be skipped by supplying the --no_fetch flag, however -d will
still be needed.  

Note that this is done concurrently with building the gene_att file.  If this step is skipped via --no_fetch, and panoct is to be run, the
user must supply one with --gene_att_file.

=head2 BLAST

<External Data>
An all-vs-all blast file of the proteins to be clusters is expected to be passed in with the option -b. 
The all-vs-all blast must be in the NCBI Blast tab delimited -m9/-m8 output.

<JCVI specific>
Runs all-vs-all blast on the fasta sequences.  Uses formatdb on the combined.fasta
file found in the fasta_dir directory, then blasts with the same file against the new db.  This
step can be passed by supplying a --blast_file of concatenated -m9 formatted blast results.

=head2 Clustering

Third, runs any clustering methods supplied by the --methods parameter, containing a comma-
sepsrated list of method names.  Currently implemented names include:

    panoct [Default]
    none [skips this step.]

Results are placed in the working dir under the directory results, named "method_name.result".

=head2 Generate Statistics

Finally, the statistics program can be run at the end of the pipeline unless --no_stats is set. 

=head1  INPUT

=head2 Database File <JCVI specific>

File listing the username/password/server combination that allows for JCVI specific database
logins. The format is username, password and databse server. Seperated by newline. 

Example file:
username
password123
SYBPROD

=head2 Fasta Files

<External data>
A protein fasta file of all proteins to be clustered must be stored in a directory called "fasta_dir". The 
fasta_dir directory needs to be located in your working directory and the fasta file must be called "combined.fasta".

<JCVI specific>
Protein fasta files will be pulled from the databases passed into this program, default. 

If supplying a fasta file provide the option --no_fetch.

If you use the --no_fetch and are running PanOCT you must pass in an associated --gene_att_file.

=head2 Gene attribute file

<External data>
To create the file each gene will have a corresponding line with 
the following information:

Molecule Name/ID<tab>Gene Identifier<tab>End 5<tab>End 3<tab>Protein Name<tab>Database/Genome<tab>Protein length<tab>TIGRFAM Role ID<tab>HMMs

Note: The TIGRFAM role id and HMM are optional.
The TIGRFAM role id can be substituted for any other classification system you'd wish to later refine your output by. If you
choose not to supply this information put a tab in the column.

=head2 Database List

The genome_list_file is a file that contains what genomes you'd like to compare as well as where the data
should be pulled from.

Format should be:
display name<tab>source<tab>asmbl_id(if JCVI)

Display name: Name of genome showed in the output files
Source: JCVI SGD
Asmbl_id(JCVI specific): If none is specified ISCURRENT is default

Note: The display name must match the genome name found in the gene attribute
file if running PanOCT.

JCVI Specific Note: If only the display name is given then that will be assumed to be the
SGD and the name to be shown

If you are not pulling data from SGD (--no_fetch) and are providing the script with your own
fasta file and gene attribute file then genome_list_file is still required, however, you must
only supply it with the genome names found in the gene attribute file.

Examples of genome_list_file:

  Just JCVI DBs
  % cat genome_list.file
  ntab08
  ntab16
  ntab17
  ntab20

  Used with --no_fetch, these names match what's in the gene attribute file
  & cat genome_list.file
  Acinetobacter_baumannii_ATCC_17978_uid17477
  Acinetobacter_baumannii_AYE_uid28921
  Acinetobacter_baumannii_AB0057_uid21111
  Acinetobacter_baumannii_uid13001

=head2 HMM File

File that contains HMM's you find interesting in your analysis. If a HMM file is provided two additional output files
will be generated. A list of clusters that contain your HMMs of interest. In the overview_stats file stats will be printed
of the number/percentage of HMMs found in the clusters as well as HMMs that were not found.

Example HMM File:
TIGR01391                     
TIGR00344                     
PF00750                       
TIGR00459                     
TIGR00152                     
TIGR00435

=head1  OUTPUT

=head2 Log File

A log file is produced containing various messages, by default this is found at
working_dir/logs/run_panoct.log
But can be specified using --log_file.  The level of messages (verbosity) can be
specified using the --debug option, passing in an integer value. 0 results in the
fewest messages, increasing the number results in more detailed explanation.

=head2 Fasta Files

The fasta files created in the first step include multifastas per genome of each
genomes ORFS, and a multifasta that combines all of those into one 'combined.fasta'.
These files are found in working_dir/fasta_dir.  This location is currently non-configurable.

=head2 BLAST Db Files

During the second step, a blast database is created in the working directory.  
The combined.blast file is the output of the blast run with the -m9 flag set, 
which produces tabular output featuring headers between each input sequence.

=head2 Clustering Results

Results for each clustering algorithm are parsed and stored as 'method_name.result'
in the working_dir/results. All files created from a clustering method are also stored
in this directory. 

=head2 Stats Files

If statistics were run these files are found in the working dir under the results
directory. A description of the files can be found in doc directory of your Pangenome 
code repository.

=head2 Graphics Files

Unless B<--no_graphics> is invoked, the pipeline produces various images summarizing
the pangenome.

=head1  CONTACT

    Jason Inman
    jinman@jcvi.org

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long  qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Capture::Tiny qw{ capture capture_merged };
use IO::File;
use Cwd;
use Cwd 'abs_path';
use feature qw(switch);
use TIGR::Foundation;
use TIGR::FASTAreader;
use TIGR::FASTArecord;
use TIGR::FASTAwriter;
use File::Basename;
use File::Slurp;
use File::Temp;
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use SeqToolBox::SeqDB;
#use Pangenome::DB::ProkSybase;
use Pangenome::ConsistencyChecks;
use Pangenome::Download;
use Data::Dumper;
use lib "$FindBin::Bin/../lib";
use File::Path qw(mkpath remove_tree);
use File::Glob qw(glob);
use File::Spec;
use File::Copy; 
use grid_tasks;

# Executables
my $BIN_DIR = $FindBin::Bin;
my $BLAST_BIN_DIR           = '/usr/local/packages/ncbi-blast+/bin'; #JCVI SPECIFIC
my $MAKEBLASTDB_EXEC        = "$BLAST_BIN_DIR/makeblastdb";          #JCVI SPECIFIC
my $BLASTP_EXEC             = "$BLAST_BIN_DIR/blastp";               #JCVI SPECIFIC
my $BLASTN_EXEC             = "$BLAST_BIN_DIR/blastn";               #JCVI SPECIFIC
my $PANOCT_EXEC             = "$BIN_DIR/panoct.pl";
my $STATS_EXEC              = "$BIN_DIR/pangenome_statistics.pl";
my $SPLIT_FASTA             = "$BIN_DIR/split_fasta.pl";
my $ATT_DATA_EXEC           = "$BIN_DIR/convert_att_file_to_hsh.pl";
my $MATCHTABLE_EXEC         = "$BIN_DIR/paralog_matchtable.pl";
my $COMPUTE_PANGENOME_EXEC  = "$BIN_DIR/compute_pangenome.R";
my $NEW_PLOT_PANGENOME_EXEC = "$BIN_DIR/new_plot_pangenome.R";
my $CORE_CLUSTER_HISTO_EXEC = "$BIN_DIR/core_cluster_histogram.R";
my $GENE_ORDER_EXEC         = "$BIN_DIR/gene_order.pl";
my $PAN_CHROMO_EXEC         = "$BIN_DIR/pan_chromosome/make_pan_chromosome_fig.sh";
my $ANNOT_EXEC              = "$BIN_DIR/get_go_annotations.pl";
my $MAPGOVIABLAST_EXEC      = "$BIN_DIR/map_go_via_blast.pl";
my $MAP_FGI_ANNOT_EXEC      = "$BIN_DIR/map_annotation_to_fgis.pl";
my $GENOME_PROP_DIR         = "$BIN_DIR/genome_properties";
my $ALIGN_CLUSTER_EXEC      = "$BIN_DIR/align_clusters.pl";
my $CHECK_COMBINED_EXEC     = "$BIN_DIR/check_combined_files.pl";
my $GP_EXEC                 = "$GENOME_PROP_DIR/genome_properties.pl";
my $EVALUATION_ORDER        = "$GENOME_PROP_DIR/evaluation_order";

$ENV{RUBYLIB}           = "$FindBin::Bin/../lib/ruby";
my $TREE_SCRIPT_EXEC    = "$BIN_DIR/make_BSR_trees.sh";

# Defaults
my $DEFAULT_WORKING_DIR = getcwd;
my $DEFAULT_DEBUG_LVL   = 0;
my $DEFAULT_MAX_TARGETS = 500;

# Various globals
my $genome_list_file    = '';
my @genomes         = ();
my $blast_file_name     = ''; # for example: 'combined.blast'
my $blast_file_dir      = ''; # for example: '/usr/local/scratch'
my $blast_file_path     = ''; # for example: '/usr/local/scratch/combined.blast'
my $method              = "panoct";
my $project_code        = '';
my $working_dir         = '';
my $results_dir         = '';
my $log_dir             = '';
my $log_file            = '';
my $lfh                 = undef;
my $debug               = '';
my $gene_att_file       = '';
my $att_dir             = '';
my $combined_fasta      = '';
my $fasta_dir           = '';
my $role_lookup_file    = '';
my $asmbl_ids           = '';
my $strict              = '';
my $att_suffix          = 'patt';
my $use_nuc             = '';

my %opts;
GetOptions( \%opts,
    'genome_list_file|g=s',
    'blast_file|b=s',
    'panoct_local',
    'blast_local',
    'working_dir|w=s',
    'gene_att_file|a=s',
    'att_dir|A=s',
    'combined_fasta|f=s',
    'fasta_dir|F=s',
    'log_dir=s',
    'no_grid',
    'no_stats',
    'no_new_plot',
    'no_panoct',
    'no_trees',
    'no_annotation',
    'no_graphics',
    'no_align_clusters',
    'reference|r=s',
    'role_ids|i=s',
    'terms|t=s',
    'strict|s=s',
    'percent_id=i',
    'hmm_file=s',
    'role_lookup=s',
    'project_code|P=s',
    'use_nuc',
    'lite',
    'debug=i',
    'help|h',
    ) || die "Error getting options! $!";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

check_params();

#Step 1: Call Blast
if ( $blast_file_path ) {

    _log( "Using $blast_file_path as blast input.  Skipping blast", 0 );
    
} else {
    
    if ( $opts{ no_panoct } && ( -s "$working_dir/combined.blast" ) ) {
        $blast_file_path = "$working_dir/combined.blast";
    } else {
        blast_genomes();
    }    

}

#Step 2: Run panoct
if ( $opts{ no_panoct } ) {
    
    _log( "Skipping panoct intentionally", 0 );
    find_result_files();
    
} else {
    
    #TODO: replace with direct call to panoct
    call_panoct()

}

#Step 3: Move result files to final locations
clean_move_files();

exit(0) if ( $opts{lite} );

#Step 4: Generate stats
unless ( $opts{no_stats} ) {
    
    call_stat_generation_script();
    move( "$working_dir/formatdb.log", "$log_dir/formatdb.log" );

}

#Step 5: Call paralog_matchtable & compute_pangenome.R & new_plot_pangenome.R
if ( (-s "$results_dir/paralogs.txt" ) && ( -s "$results_dir/matchtable_0_1.txt" ) ) {

    # setup the directory structure
    my $plot_dir = "$results_dir/R.plots";
    my $data_dir = "$plot_dir/data";
    mkdir "$plot_dir";
    mkdir "$data_dir";

    call_paralog_matchtable();

    call_compute_pangenome();

    call_new_plot_pangenome() unless $opts{ no_new_plot };

}


#Step 6: (Optionally) call tree-making scripts.
unless ( $opts{ no_trees } ) {

   call_tree_making_script(); 

}



#Step 7: Call graphics scripts.
unless ( $opts{ no_graphics } ) {

    call_graphics_scripts();

}


# Step 8: Annotate the clusters.
unless ( $opts{ no_annotation } ) {

    #call_genome_properties();
    call_annotation_scripts();

}


#Step 9: Call align_clusters.pl
unless ( $opts{ no_align_clusters } ) {

    call_align_clusters( $results_dir );

}

_log( "Pangenome Pipeline Finished.", 0 );

exit(0);


sub create_core_cluster_data {
# generate the data for the core cluster histogram from the overview_stats.txt file
# original process was a oneliner: 
# cat overview_stats.txt | perl -ne 'chomp; if ($p) {print "$_\n";} if (/^Cluster Size Breakdown/) {$p = 1;}' > R.plots/data/core_cluster_histogram_data.txt
    
    my ( $data_dir ) = @_;

    my $overview_stats = "$results_dir/overview_stats.txt";
    my $core_cluster_data = "$data_dir/core_cluster_histogram_data.txt";

    open( my $osfh, '<', $overview_stats )      || _die( "Can't open overview_stats.txt: $!", __LINE__ );
    open( my $ccfh, '>', $core_cluster_data )   || _die( "Can't write to $core_cluster_data: $!", __LINE__ );

    my $print_now = 0;
    while( <$osfh> ) {

        if ( $print_now ) {
            print $ccfh $_;
        } else {
            $print_now++ if ( /^Cluster Size Breakdown/ );
        }
        
    }

    return $core_cluster_data;

}


sub call_core_cluster_histogram {
# Takes a single argument, the path to core_cluster_histogram_data.txt

    my ( $core_cluster_data, $plots_dir ) = @_;

    my @cmd = ( $CORE_CLUSTER_HISTO_EXEC, '-i', $core_cluster_data );

    _log( "Creating core_cluster_histogram.pdf:\n". join( ' ', @cmd ), 0 );

    my $lf = "$log_dir/core_cluster_histogram.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Couldn't run $CORE_CLUSTER_HISTO_EXEC: $? $!\n", __LINE__ );

    } stdout => $lh;

    if ( -s "$working_dir/core_cluster_histogram.pdf" ) {
        move( "$working_dir/core_cluster_histogram.pdf", $plots_dir ) || _log( "WARNING: Problem moving core_cluster_histogram.pdf to $plots_dir: $!", 0 );
    } else {
        _log( "WARNING: Missing core_cluster_histogram.pdf; should be at:\n$working_dir/core_cluster_histogram.pdf", 0 )
    }

}


sub create_fGI_info {

    my ( $result_dir ) = @_;

    my $fgi_dir = "$result_dir/fGIs";
    mkdir $fgi_dir || _die( "Couldn't create $fgi_dir: $!", __LINE__ );
 
    my @cmd =   (   $GENE_ORDER_EXEC, '-T', '2', '-W', "$result_dir/cluster_weights.txt",
                    '-M', "$result_dir/95_core_adjacency_vector.txt",
                    '-m', "$result_dir/0_core_adjacency_vector.txt", '-l', '5',
                    '-t', "$genome_list_file", '-A', "$fgi_dir/Core.att",
                    '-a', "$fgi_dir/fGI.att", '-I', "$fgi_dir/fGI_report.txt",
                    '-C', "$result_dir/centroids.fasta" );
    push( @cmd, '-P' ) unless ( $opts{ nucleotide } );

    _log( "Running gene_order.pl:\n" . join( ' ', @cmd ), 0 );

    my $of = "$fgi_dir/consensus.txt";
    my $oh = IO::File->new( $of, "w+" ) || _die( "Couldn't open $of for writing: $!", __LINE__ );

    my $lf = "$log_dir/gene_order.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture{

        system( @cmd ) == 0 || _die( "Couldn't run $GENE_ORDER_EXEC: $? $!", __LINE__ );

    } stdout => $oh, stderr => $lh;

    # make_pan_chromosome needs shared_clusters.txt
    my ( $source, $target ) = (  "$result_dir/shared_clusters.txt", "$fgi_dir/shared_clusters.txt" );
    if ( eval { symlink("",""); 1 } ) {  # <-- This is a check to see if we can create symlinks on this system.
        symlink $source, $target;
    } else {
       copy( $source, $target ) || _die( "Couldn't copy $source! $!", __LINE__ );
    }

}


sub call_make_pan_chromosome {

    my @cmd = ( $PAN_CHROMO_EXEC, $working_dir );

    _log( "Running make_pan_chromosome:\n" . join( ' ', @cmd ), 0 ); 

    my $lf = "$log_dir/make_pan_chromosome.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Couldn't run $PAN_CHROMO_EXEC: $? $!", __LINE__ );

    } stdout => $lh;

}


sub call_align_clusters {
# Execute the align_clusters.pl script.

    my @cmd = ( $ALIGN_CLUSTER_EXEC, '-m', "$results_dir/matchtable.txt",
                '-c', $combined_fasta, '-s', "$results_dir/singletons_clusters.txt",
                '-w', $results_dir, 
                '-l', "$log_dir/align_clusters.log" );
    if ( $project_code ) {
        push( @cmd, '-P', $project_code );
    } else {
        push( @cmd, '--no_grid' );
    }

    _log( "Running align_clusters.pl:\n" . join( ' ', @cmd ), 0 );

    system( @cmd ) == 0 || _die( "Couldn't run $ALIGN_CLUSTER_EXEC: $? $!", __LINE__ );

}


sub call_graphics_scripts {
# Execute the following collection of R, shell, and perl scripts:
# 1. core_cluster_histogram.R
# 2. make_pan_chromosome.sh

    my $plots_dir = "$results_dir/R.plots";
    my $data_dir = "$plots_dir/data";

    _log( "Creating various graphical outputs.", 0 );

    # create core_cluster_histogram_data.txt
    my $core_cluster_data = create_core_cluster_data( $data_dir );

    # run core_cluser_histogram.R
    call_core_cluster_histogram( $core_cluster_data, $plots_dir );

    # create fGI data for make_pan_chromosome.sh
    create_fGI_info( $results_dir );

    # run make_pan_chromosome.sh
    call_make_pan_chromosome( $results_dir );


}


sub call_genome_properties {

    # TODO: Make sure centroids.pep exists for nuc runs.

    my @cmd = ( $GP_EXEC, '-all', '-seqs', 'centroids.pep', '-eval_order', $EVALUATION_ORDER, '-name', 'genome_props' );

    _log( "Running Genome Properties script:\n" . join( ' ', @cmd ), 0 );

    my $lf = "$log_dir/genome_properties.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Can't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Problem running Genome Properties script", __LINE__ );

    } stdout => $lh;

}


sub call_annotation_scripts {

    # First up: run hmms and rgi card on centroids
    my @cmd = ( $ANNOT_EXEC, '-i', "$results_dir/centroids.fasta", '-w', $results_dir, '-l', $log_dir );
    if ( $project_code ) {
        push( @cmd, '-P', $project_code );
    } else {
        push( @cmd, '--hmm_local' );
    }
    push( @cmd, '--nuc' ) if ( $opts{ use_nuc } );

    _log( "Running annotation script:\n" . join( ' ', @cmd ), 0 );

    my $lf = "$log_dir/get_go_annotations.pl.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Can't open $lf for logging: $!", __LINE__ );
    capture_merged{

        system( @cmd ) == 0 || _die( "Problem running annotation script", __LINE__ );

    } stdout => $lh;

    # Next, use blast to look for plasmid seqeunces:
    @cmd = ( $MAPGOVIABLAST_EXEC, '-i', "$results_dir/centroids.fasta", '-o', "$results_dir/plasmid_blast.out" );
    push( @cmd, '-n' ) if ( $opts{ use_nuc } );
    if ( $project_code ) {
        push( @cmd, '-P', $project_code );
    } else {
        push( @cmd, '--blast_local' );
    }

    _log( "Running blast annotation script:\n" . join( ' ', @cmd ), 0 );

    my $blf = "$log_dir/map_go_via_blast.pl.log";
    my $blh = IO::File->new( $blf, "w+" ) || _die( "Can't open $blf for logging: $!", __LINE__ );
    capture_merged{

        system( @cmd ) == 0 || _die( "Problem running blast annotation script", __LINE__ );

    } stdout => $blh;

    # Finally, use the previously run HMM/RGI and Blast results for annotating fGIs:
    my $fgi_report_index = "$results_dir/fGIs/fGI_report.txt.index";
    my $cluster_roles    = "$results_dir/centroids.cluster_roles.txt";
    my $plasmid_blast    = "$results_dir/plasmid_blast.out";
    my @input_files;
    push @input_files, $cluster_roles if -s $cluster_roles;
    push @input_files, $plasmid_blast if -s $plasmid_blast;
    if ( scalar @input_files ) {

        @cmd = ( $MAP_FGI_ANNOT_EXEC, '-i', $fgi_report_index, '-c', join(',', @input_files), '-o', "$results_dir/fGI_annotation.txt" );
        _log( "Running fGI annotation mapping:\n" . join( ' ', @cmd ), 0 );

        my $falf = "$log_dir/map_annotation_to_fgis.pl.log";
        my $falh = IO::File->new( $falf, "w+" ) || _die( "Can't open $falf for logging: $!", __LINE__ );
        capture_merged{

            system( @cmd ) == 0 || _die( "Problem running fgi annotation mapping", __LINE__ );

        } stdout => $falh;

    } else {
        _log( "No need to run fgi annotation as neither $cluster_roles nor $plasmid_blast have any data", 0 );
    }

}


sub call_tree_making_script {

    my @cmd = ( $TREE_SCRIPT_EXEC, $working_dir, $BIN_DIR );

    _log( "Running tree building:\n" . join( ' ', @cmd ), 0 );

    my $output_file_base = fileparse($TREE_SCRIPT_EXEC, qr/\.[^.]*/);
    my $lf = "$log_dir/$output_file_base.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Can't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Failure running tree-making script", __LINE__ );

    } stdout => $lh;

}


sub clean_move_files {

    #Move blast output
    move( "$working_dir/blast.stdout", "$log_dir/blast.stdout" );
    move( "$working_dir/blast.stderr", "$log_dir/blast.stderr" );

    #Move any logs
    my @logs = glob( "$working_dir/*.log" );
    for my $log ( @logs ) {
        move( $log, "$log_dir/$log" );
    }

    #move role_id lookup
    move( "$working_dir/role_id_lookup.txt", "$results_dir/role_id_lookup.txt" );

    #move method files
    move( "$working_dir/panoct.stderr",     "$log_dir/panoct.stderr" );
    move( "$working_dir/panoct.stdout",     "$log_dir/panoct.stdout" );
    move( "$working_dir/panoct_script.sh",  "$log_dir/panoct_script.sh" );
    move( "$results_dir/panoct.log",        "$log_dir/panoct.log" );

    #remove duplicate files in result directory
    unlink glob( "$results_dir/combined*" );
    unlink glob( "$results_dir/*pep" );

    my ( $db_name, $path, $fix ) = fileparse( $genome_list_file );
    unlink "$results_dir/$db_name";

}


sub create_att_dat_file {
# given a path to a att_file, create a .dat from it.

    my ( $att_file, $out_dir ) = @_;

    my $att_dat_file = $att_file . ".dat";

    my @cmd = ( $ATT_DATA_EXEC, $att_file, $out_dir );

    _log( "Creating att.dat file:\n" . join( ' ', @cmd ), 1 );

    my $lf = "$log_dir/convert_att_file_to_hsh.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Failure creating att.dat file", __LINE__ );

    } stdout => $lh;

    return $att_dat_file;

}


sub call_stat_generation_script {

    my $frameshift_file = "$results_dir/frameshifts.txt";
    my $centroid_file   = "$results_dir/centroids.fasta";
    my $fragment_file   = "$results_dir/fragments_fusions.txt";

    unless ( -s "$gene_att_file.dat" ) {
        create_att_dat_file( $gene_att_file, $working_dir );
    }

    if ( -s $gene_att_file && -s "$results_dir/panoct.result" ) {

        if ( -s $gene_att_file . ".dat" ) {

            my $params = "-a $gene_att_file.dat";
            $params .= " -m $results_dir/panoct.result";
            $params .= " -l $genome_list_file";
            $params .= " -f $frameshift_file" if (-s $frameshift_file);
            $params .= " -c $centroid_file" if (-s $centroid_file);
            $params .= " --fusion $fragment_file " if (-s $fragment_file);
            $params .= " -s $fasta_dir/combined.fasta" unless $opts{ use_nuc };
            if ( $opts{ use_nuc } ) {
                $params .= " -n $fasta_dir/combined.fasta"
            } elsif ( -s "$fasta_dir/combined.seq" ) {
                $params .= " -n $fasta_dir/combined.seq"
            }
            $params .= " -r \"$opts{ reference }\"" if $opts{ reference };
            $params .= " -i \"$opts{ role_ids }\"" if $opts{ role_ids };
            $params .= " -t \"$opts{ terms }\"" if $opts{ terms };
            $params .= " -o $results_dir";
            $params .= " --hmm_file $opts{ hmm_file }" if $opts{ hmm_file }; 
            $params .= " --role_lookup $role_lookup_file" if $role_lookup_file;

            my $cmd = $STATS_EXEC . " $params";

            _log( "Running stat generation:\n$cmd", 0 );

            my $lf = "$log_dir/pangenome_statistics.log";
            my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

            capture_merged{

                system( $cmd ) == 0 || _die ( "Error running stats: $? $!", __LINE__ );

            } stdout => $lh;

        } else {
            _die( "No $gene_att_file data file found in $working_dir", __LINE__ );
        }

    } else {
        _die( "No $gene_att_file or no result file found in $working_dir", __LINE__ );
    }

}



sub make_input_fasta_string {
# when pointed to a working_dir, retrieves *.pep and concatenates the filenames
# into a single-space-seperated list
    
    my $working_dir = shift;

    my $ext = ( $opts{ use_nuc } ) ? '.nuc' : '.pep';
    
    my $fasta_files = join( " ", glob( "$fasta_dir/*$ext" ) );
    return $fasta_files;
    
}


sub filter_panoct_output {
# Given a panoct.result file, strip away the "_BDM" and "_SYN" appendages from the
# cluster member names.
    
    my $panoct_result_file = shift;
    my $panoct_result_temp = "$working_dir/panoct.result.tmp";

    open( my $panoct, "<", $panoct_result_file ) || _die ( "Can't open $panoct_result_file: $!", __LINE__ );
    open( my $ofh, ">", $panoct_result_temp ) || _die ( "Can't open $panoct_result_temp: $!", __LINE__ );

    while ( <$panoct> ) {

        $_ =~ s/_(BDM|SYN)//g;
        print $ofh $_;

    }
    
    close $ofh;
    close $panoct;

    move( "$working_dir/panoct.result.tmp", $panoct_result_file);
    
}


sub call_paralog_matchtable{

    my $matchtable_file = "$results_dir/matchtable_0_1.txt";
    my $paralogs_file   = "$results_dir/paralogs.txt";
    my $out_file        = "$results_dir/matchtable_paralog.txt";

    my @cmd = ( $MATCHTABLE_EXEC, '-M', $matchtable_file, '-P', $paralogs_file );

    _log( "Calling paralog_matchtable.pl:\n" . join( ' ', @cmd ), 0 );

    my $oh = IO::File->new( $out_file, "w+" ) || _die( "Couldn't open $out_file for output: $!", __LINE__ );
    my $lf = "$log_dir/paralaog_matchtable.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture{

        system( @cmd ) == 0 || _die( "Couldn't run script: $? $!", __LINE__ );

    } stdout => $oh, stderr => $lh;

}


sub call_compute_pangenome {
# compute_pangenome.R takes five arguments: 
# 1. an input file which is the output of paralog_matchtable.pl
# 2. an output file name (or base output file name that it can add suffixes too)
# 3. the percent threshold to be considered core (I would suggest 95%)
# 4. the percent threshold to be considered a new gene (I would suggest 0%)
# 5. the number of combinations to generate when the number of genomes is large 
#   (I would suggest somewhere between 100-500 depending on how fast we want this to run). 

    my $input_file      = "$results_dir/matchtable_paralog.txt";
    my $output_file     = "$results_dir/pangenome_size";
    my $core_threshold  = 95;
    my $new_threshold   = 0;
    my $combinations    = 250;

    my @cmd = ( $COMPUTE_PANGENOME_EXEC, '-i', $input_file, '-o', $output_file );
    push ( @cmd, '-p', $core_threshold, '-q', $new_threshold, '-s', $combinations );

    _log( "Running compute pangenome:\n" . join( ' ', @cmd ), 0 );

    my $lf = "$log_dir/compute_pangenome.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _die( "Couldn't run script: $? $!", __LINE__ );

    } stdout => $lh;

}


sub call_new_plot_pangenome {
# new_plot_pangenome.R just takes two parameters:
# 1. an input file (the output file from compute_pangenome.R)
# 2. a base output file name that it adds suffixes to for the 
#    various outputs such as pdf files that it generates.

    my $input_file = "$results_dir/pangenome_size";
    my $output_file_basename = "$results_dir/R.plots/new_plot";

    my @cmd = ( $NEW_PLOT_PANGENOME_EXEC, '-i', $input_file, '-o', $output_file_basename );

    _log( "Running new plot pangenome script:\n" . join( ' ', @cmd), 0 );

    my $lf = "$log_dir/new_plot_pangenome.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( @cmd ) == 0 || _log( "Couldn't run new plot pangenome script: $? $!", __LINE__ );

    } stdout => $lh;

}


sub call_panoct {

    my $fasta_string = make_input_fasta_string( $working_dir );
    my @params;

    # If the combined.fasta is not in the 'results' directory, it must be
    # Copied (NOT symlinked? :/ ) to this directory.
    my $target = "$results_dir/combined.fasta";
    unless ( abs_path( $combined_fasta ) eq abs_path( $target ) || -e $target ){
        _log( "Copying combined.fasta from $combined_fasta to $target", 0 );
        copy( $combined_fasta, $target ) || _die( "Couldn't copy $combined_fasta to $results_dir! $!", __LINE__ );
    }

    # Likewise with the genome list file
    $target = "$results_dir/genomes.list";
    unless ( abs_path( $genome_list_file ) eq abs_path( $target ) || -e $target ) {
        _log( "Copying genomes.list from $genome_list_file to $target", 0 );
        copy( $genome_list_file, $target ) || _die( "Couldn't copy $genome_list_file to $results_dir!", __LINE__ );
    }

    # Likewise with the combined.att file
    $target = "$results_dir/combined.att";
    unless ( abs_path( $gene_att_file ) eq abs_path( $target ) || -e $target ) {
        _log( "Copying combined.att from $gene_att_file to $target", 0 );
        copy( $gene_att_file, $target ) || _die( "Couldn't copy $gene_att_file to $results_dir!", __LINE__ );
    }

    # Likewise with the combined.blast file
    $target = "$results_dir/combined.blast";
    unless ( abs_path( $blast_file_path ) eq abs_path( $target ) || -e $target ) {
        _log( "Copying combined.blast from $blast_file_path to $target", 0 );
        copy( $blast_file_path, $target ) || _die( "Couldn't copy $blast_file_path to $results_dir!", __LINE__ );
    }

    @params = ('-b', $results_dir, '-t', 'combined.blast', '-f', 'genomes.list',
            '-g', 'combined.att', '-P', 'combined.fasta' ); 
    push( @params, '-i', $opts{ percent_id } ) if $opts{ percent_id };
    push( @params, '-S', $strict ) if $strict;
    push( @params, '-L', '1', '-M', 'Y', '-H', 'Y', '-V', 'Y', '-N', 'Y', '-F', '1.33', '-G', 'y', '-c', '0,50,95,100', '-T' );

    my @grid_jobs;

    unless ( $opts{ panoct_local } || $opts{ no_grid } ) {

        my $panoct_script = write_grid_script( 'panoct', join( ' ', @params ) );
        _log( "Running panoct on grid.  See $panoct_script for invocation.", 0 );
        push ( @grid_jobs, launch_grid_job( $project_code, $working_dir, $panoct_script, 'panoct.stdout', 'panoct.stderr', 'himem' ) );

    } else {

        # Set up capture of panoct's stdout and stderr.  Using Tiny::Capture because simple redirecting wasn't working.
        # Probably more portable, not that it matters with all the other non-portable perl we have going on. #LINUXorBUST
        my $panoct_log = "$log_dir/panoct.log";
        my $mfh = IO::File->new( $panoct_log, "w+" ) || _die( "Can't open $panoct_log: $!", __LINE__);;

        my @cmd = ( $PANOCT_EXEC, @params);
 
        _log( "Running  panoct:\n" . join( ' ', @cmd ), 0 );
        capture_merged {
            system( @cmd ) == 0 || _die("ERROR: Problem running panoct. Check $panoct_log", __LINE__);
        } stdout => $mfh;

        unless ( -s "$log_dir/panoct.log" ) {
            _die("ERROR: Problem running panoct. Check logs/panoct.log", __LINE__);
        }

        chdir $working_dir;

    }

    wait_for_grid_jobs( \@grid_jobs ) if ( scalar @grid_jobs );

    # panoct.result is a filtered copy of matchtable.txt:
    open( my $mfh, '<', "$results_dir/matchtable.txt" ) || _die( "Can't read matchtable.txt: $!", __LINE__ );
    open( my $rfh, '>', "$results_dir/panoct.result" )  || _die( "Can't write to panoct.result: $!\n", __LINE__ );

    while ( <$mfh> ) {

        chomp;
        $_ =~ s/_(BDM|SYN)//g;
        my @f = split( /\t/, $_ );
        my @real_f;

        for my $i ( @f ) {

            next if $i =~ /^\-{2,}/;
            push @real_f, $i;

        }

        print $rfh join("\t", @real_f), "\n";

    }

    close $mfh;
    close $rfh;

}


sub write_grid_script {
# Given a command and param string, creates a shell script for running on the grid.

    my ( $command, $param_string ) = @_;

    my $script_name = "$working_dir/$command".'_script.sh';

    my $cmd_string = "$PANOCT_EXEC $param_string";

    open ( my $gsh, '>', $script_name ) || _die( "Can't open $script_name: $!", __LINE__ );

    select $gsh;

    print "#!/bin/tcsh\n\n";
    print "chdir $working_dir\n\n";
    print "$cmd_string\n";

    close $gsh;
    select STDOUT;

    chmod 0755, $script_name;

    return $script_name;

}


sub find_result_files {
# Finds all .result files in the working directory
    
    my @results = glob( "$results_dir/*.result" );
     
    if ( scalar @results == 1 ) {

        for my $file ( @results ) {

            #Check that file is not of size zero
            if ( ! -s $file ) {
                _die( "Reasult file $file is empty!", __LINE__ );
            }

        }

    } else {

        _die( "Couldn't find any result files in $results_dir", __LINE__ );

    }

}


sub get_max_target_seqs {
# Need to determine how many sequences to return for the blast hits.
# dfouts suggests num_genome * 2 
# Can't hurt to make sure all genomes are unique, too.

    my ( $genome_list ) = @_;

    my %seen;

    open( my $dfh, '<', $genome_list ) || _die( "Can't open $genome_list: $!", __LINE__ );
    while ( <$dfh> ) {
        chomp;
        $seen{$_}++;
    }

    my @dupes = grep { $seen{$_} > 1 } keys %seen;
    if ( scalar @dupes ) {
        _die( "Found duplicate genomes in $genome_list:\n" . join("\n",@dupes) , __LINE__ ) ;
    }

    my $target = 2 * scalar( keys %seen );

    return ($target > $DEFAULT_MAX_TARGETS) ? $target : $DEFAULT_MAX_TARGETS;

}


sub blast_genomes {
   # Assuming working_dir/fasta_dir/combined.fasta exists:
   # Format the blast db based on combined.fasta
   # Run all-vs-all blastp, with -m9 and filtering OFF
    
    _log( "Formatting blastdb", 1 );
    
    my $combined_fasta = "$fasta_dir/combined.fasta";

    my $formatdb_cmd = "$MAKEBLASTDB_EXEC -in $combined_fasta";
    $formatdb_cmd .= ( $opts{ use_nuc } ) ? ' -dbtype nucl ' : ' -dbtype prot ';

    my $lf = "$log_dir/formatdb.log";
    my $lh = IO::File->new( $lf, "w+" ) || _die( "Couldn't open $lf for logging: $!", __LINE__ );

    capture_merged{

        system( $formatdb_cmd ) == 0 || _die ( "Couldn't format blastdb: $? $!", __LINE__ );

    } stdout => $lh;

    $blast_file_name = "combined.blast";
    $blast_file_dir  = "$working_dir";
    $blast_file_path = "$working_dir/combined.blast";

    my $max_target_seqs = get_max_target_seqs( $genome_list_file ) || _die( "Can't determine number of blast hits to get by examining $opts{ genome_list_file }: $!", __LINE__ );

    my $blast_prog = ( $opts{ use_nuc } ) ? "$BLASTN_EXEC -task blastn": "$BLASTP_EXEC -task blastp";
    my $blast_cmd = "$blast_prog -db $combined_fasta -evalue 0.00001";  
    $blast_cmd .= " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" ";
    $blast_cmd .= " -qcov_hsp_perc 60 -max_target_seqs $max_target_seqs";

    if ( $opts{ blast_local } || $opts{ no_grid } ) {

        $blast_cmd .= " -query $combined_fasta -out $blast_file_path";
        
        _log( "Running blast locally:\n$blast_cmd", 0 );

        $lf = "$log_dir/blast.log";
        $lh = IO::File->new( $lf, "w+" ) || _die( "Coudln't open $lf for logging: $!", __LINE__ );

        capture_merged{

            system( $blast_cmd ) == 0 || _die ( "Error running blast: $? $!", __LINE__ );

        } stdout => $lh;

    } else {

        # Run blast on grid

        # Split fasta file first
        my $blast_dir = "$working_dir/blast";
        my $split_dir = "$blast_dir/split_fastas";
        mkpath( $split_dir ) unless ( -d $split_dir );
        my $split_cmd = "$SPLIT_FASTA -f $combined_fasta -n 1000 -o $split_dir";

        unless ( system( $split_cmd ) == 0 ) {

            # Check for errors of the type: "Expected: unique FASTA identifier" in split_fasta.log
            my $split_fasta_log = "$working_dir/split_fasta.pl.error";
            my $found_dupe = 0;
            if ( -s $split_fasta_log ) {

                open( my $sflh, '<', $split_fasta_log ) || _die ("Split fasta failed and I can't open $split_fasta_log to see why. :(\n", __LINE__);

                while( <$sflh> ) {
                    $found_dupe++ if (/Expected: unique FASTA identifier/);
                }

            }

            my $dupe_message = ( $found_dupe ) ? 
                                "It looks like duplicate locus tags are invovled.  Try running check_att_files.pl\n" :
                                "It doesn't look like duplicate locus tags are invovled.  Look closer at $split_fasta_log\n";

            _die ( "Error running split: $? $! '$split_cmd'\n$dupe_message", __LINE__ );

        }

        # Write shell script
        my @file_list = <$split_dir/split_fasta.*>;
        my $sh_file = write_blast_shell_script( $combined_fasta, $split_dir, $blast_dir, $blast_cmd );
        _log( "Running blast on the grid:\n$sh_file", 0 );
        
        # Launch blast job array, wait for finish
        my @grid_jobs;
        push( @grid_jobs, launch_grid_job( $project_code, $working_dir, $sh_file, 'blast.stdout', 'blast.stderr', "", scalar @file_list ) );
        wait_for_grid_jobs_arrays( \@grid_jobs,1,scalar @file_list ) if ( scalar @grid_jobs );

        # Cat all blast files together
        open( my $cfh, ">", "$working_dir/combined.blast" ) || _die ( "Couldn't open combined.blast for writing: $!", __LINE__ );
        _cat( $cfh, glob( "$blast_dir/blast_output.*" ) );

        if ( -s "$working_dir/combined.blast" ) {

            _log( "Removing intermediate blast files.", 1 );

            #Remove intermidiate blast directory
            remove_tree( $blast_dir, 0, 1 );
            $blast_file_name     = 'combined.blast'; 
            $blast_file_dir      = $working_dir;
            $blast_file_path     = "$working_dir/$blast_file_name";

        } else {
            _die( "Error creating blast file, check grid error log:$!", __LINE__ );
        }

    }    

}


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


sub make_combined_fasta {
# Given a genome_list file and a directory, look for each genome's .pep or .nuc file
# and combine them all into combined.fasta file in fasta_dir 

    my ( $fasta_dir ) = @_;
    my $combined_fasta = "$fasta_dir/combined.fasta";

    my ( @fasta_files, @missing );

    my $extension =  ( $opts{ use_nuc } ) ? 'nuc' : 'pep';

    for my $genome ( @genomes ) {

        my $fasta_file = "$fasta_dir/$genome.$extension";

        if ( -s $fasta_file ) {
            push @fasta_files, $fasta_file;
        } else {
            push @missing, $genome;
        }

    }

    if ( scalar( @missing ) ) {
        _die( "Missing (or empty) .$extension files for the following genomes:\n" . join( "\n", @missing ) . "\n", __LINE__ );
    } else {

        open( my $ofh, '>', $combined_fasta ) || _die( "Can't open $combined_fasta: $!\n", __LINE__ );

        for my $fasta_file ( @fasta_files ) {

            open( my $ifh, '<', $fasta_file ) || _die( "Can't open $fasta_file: $!\n", __LINE__ );

            while (<$ifh>) {
                print $ofh $_;
            }

        }

    }

    return $combined_fasta;

}


sub make_combined_att {
# Given a genome_list file and a directory, look for each genome's .att file and
# combine them all into a single combined.att file in the working directory

    my ( $att_dir ) = @_;
    my $gene_att_file = "$working_dir/combined.att";

    my ( @att_files, @missing );

    for my $genome ( @genomes ) {

        my $att_file = "$att_dir/$genome.$att_suffix";

        if ( -s $att_file ) {
            push @att_files, $att_file;
        } else {
            push @missing, $att_file;
        }

    }

    if ( scalar( @missing ) ) {
        _die( "Missing att files for the following genomes:\n" . join( "\n", @missing ) . "\n", __LINE__ );
    } else {

        open( my $ofh, '>', $gene_att_file ) || _die( "Can't oprn $gene_att_file: $!\n", __LINE__ );

        for my $att_file ( @att_files ) {

            open( my $ifh, '<', $att_file ) || _die( "Can't open $att_file: $!\n", __LINE__ );
            while (<$ifh>) {
                print $ofh $_;
            }

        }

    }

    return $gene_att_file;

}


sub check_params {

    my $errors = '';

    # General options
    $debug          = $opts{ debug }        // $DEFAULT_DEBUG_LVL;
    $working_dir    = $opts{ working_dir }  // $DEFAULT_WORKING_DIR;
    $results_dir    = "$working_dir/results";
    $log_dir        = $opts{ log_dir }      // "$working_dir/logs";
    $project_code   = $opts{ project_code };
    if ( $opts{ use_nuc } ) { $att_suffix = 'natt' }

    if ( ! ( $opts{ blast_file } || $opts{ blast_local } ) ) {

        $errors .= "MUST provide a --project_code, --no_grid, or choose between --blast_local and --blast_file\n" unless $project_code || $opts{no_grid};

    } 

    # check genome_list
    $genome_list_file = $opts{ genome_list_file } // "$working_dir/genomes.list";
    if ( -s $genome_list_file ) {
        @genomes = read_file( $genome_list_file );
        for my $genome ( @genomes ) {
            chomp $genome;
        }
    } else {
        $errors .= "Please provide a non-empty --genome_list_file\n";
    }

    #make output directories
    mkpath( "$results_dir" ) unless ( -d "$results_dir" );
    mkpath( "$log_dir" )    unless ( -d "$log_dir" );

    $log_file = "$log_dir/run_panoct.log";
    open( $lfh, '>', $log_file ) || die "Can't open log file $log_file: $!";

    # Options setting up blast or the blast file
    if ( $opts{ blast_file } ) {

        $blast_file_path = abs_path( $opts{ blast_file } );
        $errors .= "$blast_file_path doesn't exist!\n" unless ( -e $blast_file_path );
        ( $blast_file_name, $blast_file_dir ) = fileparse( abs_path( $blast_file_path ) );

    }

    # combined att_file setup
    $gene_att_file  = $opts{ gene_att_file } // '';
    if ( !( $gene_att_file) ) {
        $att_dir = $opts{ att_dir } // "$working_dir/att_dir";
        $gene_att_file = make_combined_att( $att_dir );
    }
    $gene_att_file = abs_path( $gene_att_file ) if $gene_att_file;

    # combined fasta file setup
    $fasta_dir  = $opts{ fasta_dir }  // "$working_dir/fasta_dir";
    $errors .= "The fasta dir $fasta_dir doesn't exist.\n" unless ( -d $fasta_dir );

    $combined_fasta = $opts{ fasta_file } // "$fasta_dir/combined.fasta";

    if ( !( -s $combined_fasta ) ) {
        $combined_fasta = make_combined_fasta( $fasta_dir );
    }
    $combined_fasta = abs_path( $combined_fasta ) if $combined_fasta;

    # Check to make sure the combined.fasta and combined.att
    my @cmd = ( $CHECK_COMBINED_EXEC, '-a', $gene_att_file, '-f', $combined_fasta,
                '-l', "$log_dir/check_combined_files.log" );

    if ( system(@cmd) ) { # non-zero here means an issue was found.
        $errors .= "Found an error with the gene_att_file or combined.fasta.  Please " .
                    "look at logs/check_combined_files.log to see what the issue is.\n";
    }

    # other called-script params
    if ( $opts{ strict } ) {

        my $op = lc( $opts{ strict } );
        $strict = $opts{ strict };
        $errors .= "--strict can be low, high, none\n" unless ( $op =~ /(low|high|none)/ );

    } else {
        $strict = "low";
    }

    if ( $opts{ role_lookup } ) {

        $errors .= "$opts{ role_lookup } does not exist or is size zero\n" unless ( -s $opts{ role_lookup } );
        $role_lookup_file = $opts{ role_lookup }; 

    }

    $errors .= "File is size zero or does not exist: $opts{ hmm_file }\n" if ( $opts{ hmm_file } && !( -s $opts{ hmm_file } ) );
    
    die $errors if $errors;
 
}


sub _log {

    my ( $msg, $lvl ) = ( @_ );

    $lvl = $lvl // $debug;

    chomp( $msg );
    print $lfh "\n$msg\n" if ( $lvl <= $debug );

}


sub _die {

    my ( $msg,$line ) = @_;

    clean_move_files();
    _log( $msg, 0 );
    die "$msg at line $line\n";

}
