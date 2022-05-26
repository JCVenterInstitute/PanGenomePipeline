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

get_go_annotations.pl - assign GO ids to fasta sequences

=head1 SYNOPSIS

    USAGE: get_go_annotations.pl -i input.fasta [ --nuc ] --core_att

=head1 OPTIONS

=over 

=item B<--input_fasta, -i> :   Path to input fasta, usually centroids.fasta

=item B<--nuc, -n>  :   Flag to indicate if this is a nucleotide fasta file

=item B<--core_att, -c> :   [Optional] Path to Core.att

=item B<--hmm_db, -d>   :   [Optional] Specify "pfam","tigrfams","both" [Default: "both"]

=item B<--no_go>        :   Don't look for GO terms

=item B<--no_aro>       :   Don't look for ARO terms

=item B<--no_cleanup>   :   Don't remove hmmer hits files.  [Default behavior is to delete them]

=item B<--roles2go, -r> :   Document linking labels (roles) to GO IDs.

=item B<--working_dir, -w>  :   Working directory [Default: cwd]

=item B<--project_code, -P> :   SGE Project code if this should be run on the grid.

or

=item B<--hmm_local>    : Do NOT run on the grid.  No --project_code is needed when this is used.

=back

=head1 DESCRIPTION

This script will, when given an input fasta file, run hmmer searches against
Pfam/TIGRFAM and assign GO ids to each sequence, where possible.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Getopt::Long qw( :config no_ignore_case no_auto_abbrev );
use Pod::Usage;
use Capture::Tiny qw{ capture capture_merged };
use IO::File;
use Cwd;
use Cwd 'abs_path';
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use File::Basename;
use File::Copy;
use File::Path qw(mkpath remove_tree);
use File::Glob qw(glob);
use grid_tasks;

my $BIN_DIR = $FindBin::Bin;
my $FIX_HEADERS_EXEC = "$BIN_DIR/clean_multifasta.pl";
my $TRANSLATE_ORFS_EXEC = 'transeq';
my $SPLIT_FASTA_EXEC    = "$BIN_DIR/split_fasta.pl";
my $HMMER2GO_DIR        = "$BIN_DIR/HMMER2GO";
my $HMMER2GO_EXEC       = "$HMMER2GO_DIR/bin/hmmer2go";
my $HMMER2GO_LIBDIR     = "$HMMER2GO_DIR/lib";
my $HMMER2GO_DATA       = "$HMMER2GO_DIR/data";
my $RGI_DIR             = "$BIN_DIR/rgi";
my $ORIGINAL_CARD_JSON  = "$RGI_DIR/_data/card.json";
my $RGI_EXEC            = "$RGI_DIR/jcvi.rgi.py";

my $DEFAULT_HMMDB       = 'both';
my $DEFAULT_ROLES2GO    = "$HMMER2GO_DATA/Main_roles2GO";

my $working_dir;
my $project_code;
my $log_dir;

my %opts;
GetOptions( \%opts,
            'input_fasta|i=s',
            'nuc|n',
            'core_att|c=s',
            'hmm_db|d=s',
            'project_code|P=s',
            'roles2go|r=s',
            'working_dir|w=s',
            'log_dir|l=s',
            'hmm_local',
            'no_go',
            'no_aro',
            'no_cleanup',
            'skip_searches|X',
            'help|h',
            ) || die "Can't read options! $!\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{help};

check_params();

# fix headers
my $input_fasta = fix_headers( $opts{ input_fasta }, $opts{ nuc } );

# translate orfs if --nuc
if ( $opts{ nuc } ) {
    translate_orfs( $input_fasta ) unless ( $opts{ skip_searches } );;
}

unless ( $opts{ no_go } ) {

    # run the searches
    run_hmmer_searches( $input_fasta ) unless ( $opts{ skip_searches } );

    # map GO terms
    my $final_go_map = map_go_terms();

    # create slimmed mapping
    create_slimmed_mapping( $final_go_map );

}

run_aro_searches( $input_fasta ) unless ( $opts{ no_aro } );

exit(0);


sub fix_headers {

    my ( $input_fasta, $nuc ) = @_;

    my $output_file = $input_fasta . '.cleaned';

    if ( $nuc ) {

        my @cmd = ($FIX_HEADERS_EXEC, '-i', $input_fasta, '-o', $output_file);

        my $base = (fileparse($FIX_HEADERS_EXEC, qr/\.[^.]*/ ))[0];
        my $lf = "$log_dir/$base.log"; 
        my $lh = IO::File->new( $lf, "w+" ) || die "Can't open $lf: $!\n";

        capture_merged {

            system( @cmd ) == 0 || die ( "Failed running header fixer: ", join( ' ', @cmd ), "\nCheck $lf\n" );

        } stdout => $lh;

    } else {

        # replace header lines with the > and first 'word' in the header file.
        open( my $ifh, '<', $input_fasta ) || die "Can't open $input_fasta: $!\n";
        open( my $ofh, '>', $output_file ) || die "Can't open $output_file: $!\n";

        while( <$ifh> ) {

            if ( /^(>\S+)\s/ ) {
                $_ = "$1\n";
            }

            print $ofh $_;

        }

    }

    return $output_file;
    
}


sub translate_orfs {

    my ( $input_fasta ) = @_;

    my $output_file = "$input_fasta.tmp";

    my $translation_table = 11;

    my @cmd = ( $TRANSLATE_ORFS_EXEC, '-sequence', $input_fasta, '-outseq', $output_file, '-table', $translation_table );

    my $base = (fileparse($TRANSLATE_ORFS_EXEC, qr/\.[^.]*/ ))[0];
    my $log_file = "$log_dir/$base.log";
    my $err_file = "$log_dir/$base.err";
    my $lh = IO::File->new( $log_file, "w+" ) || die "Can't open log $log_file: $!\n";
    my $eh = IO::File->new( $err_file, "w+" ) || die "Can't open error file $err_file: $!\n";

    capture{

        system( @cmd ) == 0 || die( "Failed running orf translater: ", join( ' ', @cmd ), "\n" );

    } stdout => $lh, stderr => $eh;

    post_process_orf_headers( $output_file );

    move( $output_file, $input_fasta );

}


sub post_process_orf_headers {
# Remove the '_1' that transeq puts after each orf.

    my ( $input_fasta ) = @_;

    my $output_file = "$input_fasta.tmp";

    open( my $ifh, '<', $input_fasta ) || die "Can't open $input_fasta: $!\n";
    open( my $ofh, '>', $output_file ) || die "Can't open $output_file: $!\n";

    while( <$ifh> ) {

        if ( /^(.*_\d+)_1/ ) {
            $_ = "$1\n";
        }

        print $ofh $_;

    }

    close $ofh;
    close $ifh;

    move( $output_file, $input_fasta );

}


sub copy_aro_files {

    my ( $orig_dir, $new_dir ) = @_;

    my $orig_card_json = $ORIGINAL_CARD_JSON;
    my $new_card_json  = "$new_dir/card.json";

    copy( $ORIGINAL_CARD_JSON, $new_card_json ) || die "Didn't copy card.json: $!\n";

    my @to_copy = ('proteindb.fsa','protein.db.pin','protein.db.phr','protein.db.psq','protein.db.dmnd');
    
    for my $file ( @to_copy ) {

        my $old_file = "$orig_dir/$file";
        my $new_file = "$new_dir/$file";

        copy( $old_file, $new_file ) || die "Couldn't copy $old_file to $new_file: $!\n";

    }

}


sub run_aro_searches {

    my ( $input_fasta ) = @_;

    my $aro_dir = "$working_dir/aro_searches";
    my $abs_aro_dir = abs_path($aro_dir); 
    mkdir $aro_dir || die "Can't make $aro_dir: $!\n";
    my $old_dir = getcwd();

    copy( $input_fasta, $aro_dir ) || die "Can't copy $input_fasta: $!\n";

    chdir $aro_dir || die "Can't chdir into $aro_dir: $!\n";
    $input_fasta = (fileparse($input_fasta))[0];

    my $output_json = $input_fasta . '.aro'; # Will actually be $input_fasta.aro.json
    my $output_file = 'dataSummary'; # Will actually be datasummary.txt
    my $new_card_json  = "$abs_aro_dir/card.json";

    # copy card.json & other files needed for the run. Every time.  sigh.
    copy_aro_files( "$RGI_DIR/_db", $aro_dir );

    # Run rgi.
    my @cmd = ( 'python2', $RGI_EXEC, '-t', 'protein', '-i', $input_fasta, '-o', $output_json );
    my $base = (fileparse($RGI_EXEC, qr/\.[^.]*/ ))[0];
    my $lf = "$log_dir/$base.log";
    my $ef = "$log_dir/$base.err";
    my $lh = IO::File->new( $lf, "w+" ) || die "Can't open log_file $lf: $!\n";
    my $eh = IO::File->new( $ef, "w+" ) || die "Can't open error_file $ef: $!\n";
    capture{

        system( @cmd ) == 0 || die "Error running rgi: ", join( ' ', @cmd ), "\n";

    } stdout => $lh, stderr => $eh;

    chdir $old_dir;

}


sub run_hmmer_searches {
# Execute the hmm searches

    my ( $input_fasta ) = @_;

    my $database;
    my ( @grid_jobs, @file_list, $split_dir );

    unless ( $opts{ hmm_local } ) {

        # Need to split the input fasta.
        $split_dir = "$working_dir/split_fastas";
        mkpath( $split_dir ) unless ( -d $split_dir );
        my @split_cmd = ( $SPLIT_FASTA_EXEC, '-f', $input_fasta, '-n', '1000', '-o', $split_dir, '-e', '.fasta' );
        my $base = (fileparse($SPLIT_FASTA_EXEC, qr/\.[^.]*/ ))[0];
        my $lf = "$log_dir/$base.log";
        my $lh = IO::File->new( $lf, "w+" ) || die "Can't open logfile: $lf: $!\n";
        capture_merged {
            system( @split_cmd ) == 0 || die( "Error running splitfasta: ", join( ' ', @split_cmd ), "\n" );
        } stdout => $lh;
        @file_list = <$split_dir/split_fasta.*>;

    }

    my @hmmer2go_cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'run' );

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'pfam' ) {

        $database = "$HMMER2GO_DATA/Pfam-A.hmm";
        my @cmd = ( @hmmer2go_cmd, '-d', $database );

        if ( $opts{ hmm_local } ) {
            # Run locally.
            @cmd = ( @cmd, '-i', $input_fasta );
            my $base = (fileparse($HMMER2GO_EXEC, qr/\.[^.]*/ ))[0];
            my $lf = "$log_dir/$base.run.pfam.log";
            my $lh = IO::File->new( $lf, "w+" ) || die "Can't open log file: $lf: $!\n";
            capture_merged{
                system( @cmd ) == 0 || die( "Problem running hmmer2go 'run': ", join( ' ', @cmd ), "\n" );
            } stdout => $lh;

        } else {
            # Farming to grid
            my $hmm_dir = "$working_dir/pfam_hmms";
            mkpath( $hmm_dir ) unless ( -d $hmm_dir );

            # write shell script
            my $sh_file = &write_hmm_shell_script( $input_fasta, $split_dir, $hmm_dir, \@cmd, 'pfam' );

            push( @grid_jobs, launch_grid_job( $project_code, $hmm_dir, $sh_file, "$log_dir/grid_logs/pfam_hmmer2go/", "$log_dir/grid_logs/pfam_hmmer2go/" , "", scalar @file_list ) );

        }

    }

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'tigrfams' ) {

        $database = "$HMMER2GO_DATA/TIGRFAMs_15.0_HMM.LIB";

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'run','-d', $database );

        if ( $opts{ hmm_local } ) {
            # Run Locally.
            @cmd = ( @cmd, '-i', $input_fasta );
            my $base = (fileparse($HMMER2GO_EXEC, qr/\.[^.]*/ ))[0];
            my $lf =  "$log_dir/$base.run.base.log";
            my $lh = IO::File->new( $lf, "w+" ) || die "Can't open log file: $lf: $!\n";
            capture_merged{
                system( @cmd ) == 0 || die( "Problem running hmmer2go 'run': ", join( ' ', @cmd ), "\n" );
            } stdout => $lh;

        } else {
            # Farming to grid.
            my $hmm_dir = "$working_dir/tigrfams_hmms";
            mkpath( $hmm_dir ) unless ( -d $hmm_dir );

            # write shell script
            my $sh_file = &write_hmm_shell_script( $input_fasta, $split_dir, $hmm_dir, \@cmd, 'tigrfams' );

            push( @grid_jobs, launch_grid_job( $project_code, $hmm_dir, $sh_file, "$log_dir/grid_logs/tigrfam_hmmer2go/", "$log_dir/grid_logs/tigrfam_hmmer2go/", "", scalar @file_list ) );

        }

    }

    unless ( $opts{ hmm_local } ) {

        # bide our time until the grid is done with our jobs.
        wait_for_grid_jobs_arrays( \@grid_jobs, 1, scalar @file_list ) if ( scalar @grid_jobs );

        # cat the files together by db:
        if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'pfam' ) {
            my $tblout_file = $opts{ input_fasta }. "_Pfam-A.tblout";
            open( my $pfh, '>', $tblout_file ) || die "Can't open $tblout_file: $!\n"; 
            &_cat( $pfh, glob( "$working_dir/pfam_hmms/*.tblout" ) );            

            # can delete interim files
            if ( -s $tblout_file ) {
               remove_tree( "$working_dir/pfam_hmms", 0, 1 ) unless $opts{ no_cleanup };
            } else {
                die "Problem with collecting the pfam tblout files.\n";
            }
        }
         
        if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'tigrfams' ) {
            my $tblout_file = $opts{ input_fasta }. "_TIGRFAMs_15.0_HMM.tblout";
            open( my $tfh, '>', $tblout_file ) || die "Can't open $tblout_file: $!\n";
            &_cat( $tfh, glob( "$working_dir/tigrfams_hmms/*.tblout" ) );

            # can delete interim files
            if ( -s $tblout_file ) {
               remove_tree( "$working_dir/tigrfams_hmms", 0, 1 ) unless $opts{ no_cleanup };
            } else {
                die "Problem with collecting the tigrfams tblout files.\n";
            }

        }

    }

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


sub write_hmm_shell_script {
# Create a shell script suitable for grid submission:

    my ( $input_fasta, $split_dir, $hmm_dir, $hmmer2go_cmd_ref, $db ) = @_;

    my @cmd = ( @$hmmer2go_cmd_ref, '-i', "$split_dir".'/split_fasta.$SGE_TASK_ID.fasta' );
    my $script_name = "$hmm_dir/grid_$db".'_hmmer2go.sh';

    open( my $gsh, '>', $script_name ) || die "Can't write grid script: $script_name: $!\n";

    print $gsh "#!/bin/tcsh\n\n";
    print $gsh join( ' ', @cmd ), "\n";

    close $gsh;
    chmod 0755, $script_name;

    return $script_name;

}


sub map_go_terms {
# map the go terms to the hmm hits

    my ( $tblout_file, $map_file, $final_go_map );

    my $input_fasta = $opts{ input_fasta }; 
    my $input_fasta_base = ( fileparse( $opts{ input_fasta }, qr/\.[^.]*/) )[0];

    # All runs based on the same input will have the same GOterm_mapping.tsv file
    # So let's cat the results once 
    ( my $input_base  = $opts{ input_fasta } ) =~ s/([^\.+])\..*/$1/;
    my $term_mapping_tsv = $input_base . '_GOterm_mapping.tsv';

    # So we will concatonate the results to this file for all runs in case we're setting both up.
    $final_go_map = $input_fasta_base . '.mapterm.txt';
    unlink $final_go_map if ( -f $final_go_map );

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'pfam' ) {

        $tblout_file = $input_fasta . "_Pfam-A.tblout";
        $map_file    = $input_fasta . '_pfam_GO.txt';

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'mapterms', '-i', $tblout_file, '-p', "$HMMER2GO_DATA/pfam2go.map", '-o', $map_file, '--map' );

        my $base = (fileparse($HMMER2GO_EXEC, qr/\.[^.]*/ ))[0];
        my $lf = "$log_dir/$base.mapterms.pfam.log";
        my $lh = IO::File->new( $lf, "w+" ) || die "Can't open log file: $lf: $!\n";
        capture_merged{
            system( @cmd ) == 0 || die ( "Problem running hmmer2go 'mapterms': ", join( ' ', @cmd ), "\n");
        } stdout => $lh;

        if ( -s $term_mapping_tsv ) {
            open( my $ifh, '<', $term_mapping_tsv ) || die "Can't open $term_mapping_tsv: $!";
            open( my $ofh, '>>', $final_go_map )     || die "Can't open $final_go_map: $!";

            while (<$ifh>) {
                print $ofh $_;
            }

        } else {

            die "No $term_mapping_tsv !?\n";

        }

        unless ( -s $final_go_map ) {

            die "Something happened... no $final_go_map\n";
        }

    }

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'tigrfams' ) {

        $tblout_file = $input_fasta . "_TIGRFAMs_15.0_HMM.tblout";
        $map_file    = $input_fasta . '_TIGRFAMs_15.0_HMM.GO.txt';

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'mapterms', '-i', $tblout_file, '-p', "$HMMER2GO_DATA/tigrfams2go.map", '-o', $map_file, '--map' );
        my $base = (fileparse($HMMER2GO_EXEC, qr/\.[^.]*/ ))[0];
        my $lf = "$log_dir/$base.mapterms.tigrfam.log";
        my $lh = IO::File->new( $lf, "w+" ) || die "Can't open log file: $lf: $!\n";
        capture_merged{
            system( @cmd ) == 0 || die ( "Problem running hmmer2go 'mapterms': ", join( ' ', @cmd ), "\n");
        } stdout => $lh;

        if ( -s $term_mapping_tsv ) {
            open( my $ifh, '<', $term_mapping_tsv ) || die "Can't open $term_mapping_tsv: $!";
            open( my $ofh, '>>', $final_go_map )     || die "Can't open $final_go_map: $!";

            while (<$ifh>) {
                print $ofh $_;
            }

        }
    }

    return $final_go_map;

}


sub create_slimmed_mapping {

    my ( $final_go_map ) = @_;

    open( my $rfh, '<', $opts{roles2go} ) || die "Can't open $opts{roles2go}: $!\n";

    # Cast aside header
    readline($rfh);

    my %roles;
    my %go2roles;

    # build mapping hashes for go2roles catalogue
    while ( <$rfh> ) {

        chomp;
        my ( $role, $line_go_ids ) = split( "\t", $_ );

        for my $go_term ( split( ',', $line_go_ids ) ) {

            $roles{ $role }->{ $go_term }++;
            if ( exists $go2roles{ go_term } ) {
                warn "WARNING: $go_term was already seen in $go2roles{ $go_term }\n";
            } else {
                $go2roles{ $go_term } = $role;
            }

        }

    }

    # build hash of arrays of go_terms per input id.
    open( my $fgfh, '<', $final_go_map ) || die "Can't open $final_go_map: $!\n";
    my %input_ids;

    # This way is for working off the final .mapping.txt file
    while (<$fgfh>) {

        chomp;
        my ( $id, $go_ids ) = split("\t",$_);
        my @go_id_list = split(',',$go_ids);

        push @{$input_ids{ $id }}, @go_id_list;

    }

#    # This way is for working directly off the _GO.txt files
#    while (<$fgfh>) {
#
#        chomp;
#        my ( $id, $go_id ) = ( split( "\t", $_ ) )[ 0, 4 ];
#
#        push @{$input_ids{ $id }}, $go_id;
#
#    }

    ( my $input_base  = $opts{ input_fasta } ) =~ s/([^\.+])\..*/$1/;
    my $output_file = "$input_base.cluster_roles.txt";
    open( my $ofh, '>', $output_file ) || die "Can't open $output_file: $!\n";

    # convert go_terms to the category
    for my $id ( keys %input_ids ) {

        my $category = 'Other';
        my $cat_go = '';

        for my $go_id ( @{$input_ids{ $id }} ) {

            next unless defined $go_id;

            if ( exists $go2roles{ $go_id } ) {

                $category = $go2roles{ $go_id };
                $cat_go = $go_id;

            }

        }

        print $ofh "$id\t$category\t$cat_go\n";

    }

}


sub check_params {

    my $errors = '';

    if ( $opts{ input_fasta } ) {
        $errors .= "Can't find --input_fasta file $opts{ input_fasta } (or it is empty)\n"
             unless ( -s $opts{ input_fasta } );
    } else {
        $errors .= "Please provide --input_fasta\n";
    }

    if ( $opts{ core_att } and not ( -s $opts{ core_att } ) ) {
        $errors .= "Could not find --core_att $opts{ core_att } (or it is empty)\n";
    }

    unless ( $opts{ no_go } ) {

        # Need to know which hmms to run.
        $opts{ hmm_db } = $opts{ hmm_db } // 'both';
        $opts{ hmm_db } = lc $opts{ hmm_db };
        $opts{ hmm_db } = 'tigrfams' if ( $opts{ hmm_db } eq 'tigrfam' ); # Cuz people gonna do it.

        # Need a default roles2go mappign
        $opts{ roles2go } = $opts{ roles2go } // $DEFAULT_ROLES2GO;

        if ( $opts{ project_code } ) {

            $project_code = $opts{ project_code };

        } else {

            $errors .= "Please supply a --project_code or use --hmm_local\n" unless $opts{ hmm_local };

        }

    }

    if ( $opts{ no_go } && $opts{ no_aro } ) {

        $errors .= "Don't use --no_go and --no_aro!  I'll have nothing to do!\n";

    }

    $working_dir = $opts{ working_dir } // getcwd();
    chdir $working_dir;
    $log_dir     = $opts{ log_dir } // "$working_dir/logs";
    unless ( -d $log_dir ) {
        mkdir($log_dir);
        $errors .= "Can't make log_dir $log_dir: $!\n" unless ( -d $log_dir );
    }

    die $errors if $errors;

}
