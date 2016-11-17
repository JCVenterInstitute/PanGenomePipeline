#!/usr/bin/env perl
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
use Cwd;
use Cwd 'abs_path';
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use File::Basename;
use File::Copy;
use File::Path qw(mkpath remove_tree);
use File::Glob qw(glob);


my $BIN_DIR = $FindBin::Bin;
my $FIX_HEADERS_EXEC = "$BIN_DIR/clean_multifasta.pl";
my $TRANSLATE_ORFS_EXEC = 'transeq';
my $SPLIT_FASTA_EXEC    = "$BIN_DIR/split_fasta.pl";
my $HMMER2GO_DIR        = "$BIN_DIR/HMMER2GO";
my $HMMER2GO_EXEC       = "$HMMER2GO_DIR/bin/hmmer2go";
my $HMMER2GO_LIBDIR     = "$HMMER2GO_DIR/lib";
my $HMMER2GO_DATA       = "$HMMER2GO_DIR/data";
my $RGI_DIR             = "$BIN_DIR/rgi";
my $ORIGINAL_CARD_JSON  = "$RGI_DIR/card.json.orig";
my $RGI_EXEC            = "$RGI_DIR/rgi.py";
my $RGI_CONVERT_EXEC    = "$RGI_DIR/convertJsonToTSV.py";
my $RGI_CLEAN_EXEC      = "$RGI_DIR/clean.py";

my $DEFAULT_HMMDB   = 'both';

my $working_dir;
my $project_code;

my %opts;
GetOptions( \%opts,
            'input_fasta|i=s',
            'nuc|n',
            'core_att|c=s',
            'hmm_db|d=s',
            'project_code|P=s',
            'roles2go|r=s',
            'working_dir|w=s',
            'hmm_local',
            'no_go',
            'no_aro',
            'help|h',
            ) || die "Can't read options! $!\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{help};

check_params();

# fix headers
my $input_fasta = fix_headers( $opts{ input_fasta }, $opts{ nuc } );

# translate orfs if --nuc
if ( $opts{ nuc } ) {
    translate_orfs( $input_fasta );
}

unless ( $opts{ no_go } ) {

    # run the searches
    run_hmmer_searches( $input_fasta );

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

        system( @cmd ) == 0 || die ( "Failed running header fixer: ", join( ' ', @cmd ), "\n" );

    } else {

        # replace header lines with the > and first 'word' in the header file.
        open( my $ifh, '<', $input_fasta ) || die "Can't open $input_fasta: $!\n";
        open( my $ofh, '>', $output_file ) || die "Can't open $output_file: $!\n";

        while( <$ifh> ) {

            if ( /^(>\S+)\s/ ) {
                $_ = "$1\n";;
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

    system( @cmd ) == 0 || die( "Failed running orf translater: ", join( ' ', @cmd ), "\n" );

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


sub run_aro_searches {

    my ( $input_fasta ) = @_;

    my $aro_input_fasta = (fileparse($input_fasta))[0];  # Necessary because of rgi.py pretending to know better than us where our stuff is.

    my $output_json = $aro_input_fasta . '.aro'; # Will actually be $input_fasta.aro.json
    my $output_file = 'dataSummary'; # Will actually be $input_fasta.aro.txt
    my $new_card_json  = "$working_dir/card.json";

    # copy card.json. Every time.  sigh.
    copy( $ORIGINAL_CARD_JSON, $new_card_json );

    # Go to the working_dir because rgi.py doesn't accept a working_dir as a parameter.
    my $old_dir = getcwd;
    chdir $working_dir;

    # Run rgi.
    my @cmd = ( '/usr/bin/env', 'python', $RGI_EXEC, '-t', 'protein', '-i', $input_fasta, '-o', $output_json );
    system( @cmd ) == 0 || die "Error running rgi: ", join( ' ', @cmd ), "\n";

    # Run conversion from json to tab-delimitted:
    @cmd = ( '/usr/bin/env', 'python', $RGI_CONVERT_EXEC, '-i', "$output_json.json", '-o', $output_file );
    system( @cmd ) == 0 || die "Error converting rgi json into tabbed-text: ", join( ' ', @cmd ), "\n";

    # Back to the old dir.
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
        my @split_cmd = ( $SPLIT_FASTA_EXEC, '-f', $input_fasta, '-n', '1000', '-o', $split_dir );
        system( @split_cmd ) == 0 || die( "Error running splitfasta: ", join( ' ', @split_cmd ), "\n" );
        rename_split_files( $split_dir );
        @file_list = <$split_dir/split_fasta.*>;

    }

    my @hmmer2go_cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'run' );

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'pfam' ) {

        $database = "$HMMER2GO_DATA/Pfam-A.hmm";
        my @cmd = ( @hmmer2go_cmd, '-d', $database );

        if ( $opts{ hmm_local } ) {
            # Run locally.
            @cmd = ( @cmd, '-i', $input_fasta );

            system( @cmd ) == 0 || die( "Problem running hmmer2go 'run': ", join( ' ', @cmd ), "\n" );

        } else {
            # Farming to grid
            my $hmm_dir = "$working_dir/pfam_hmms";
            mkpath( $hmm_dir ) unless ( -d $hmm_dir );

            # write shell script
            my $sh_file = &write_hmm_shell_script( $input_fasta, $split_dir, $hmm_dir, \@cmd, 'pfam' );

            push( @grid_jobs, launch_grid_job( $sh_file, 'pfam.hmmer2go.stdout', 'pfam.hmmer2go.stderr', "", scalar @file_list, $hmm_dir ) );

        }

    }

    if ( $opts{ hmm_db } eq 'both' || $opts{ hmm_db } eq 'tigrfams' ) {

        $database = "$HMMER2GO_DATA/TIGRFAMs_15.0_HMM.LIB";

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'run','-d', $database );

        if ( $opts{ hmm_local } ) {
            # Run Locally.
            @cmd = ( @cmd, '-i', $input_fasta );

            system( @cmd ) == 0 || die( "Problem running hmmer2go 'run': ", join( ' ', @cmd ), "\n" );

        } else {
            # Farming to grid.
            my $hmm_dir = "$working_dir/tigrfams_hmms";
            mkpath( $hmm_dir ) unless ( -d $hmm_dir );

            # write shell script
            my $sh_file = &write_hmm_shell_script( $input_fasta, $split_dir, $hmm_dir, \@cmd, 'tigrfams' );

            push( @grid_jobs, launch_grid_job( $sh_file, 'tigrfams.hmmer2go.stdout', 'tigrfams.hmmer2go.stderr', "", scalar @file_list, $hmm_dir ) );

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
               remove_tree( "$working_dir/pfam_hmms", 0, 1 );
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
               remove_tree( "$working_dir/tigrfams_hmms", 0, 1 );
            } else {
                die "Problem with collecting the tigrfams tblout files.\n";
            }

        }

    }

}


sub rename_split_files {
# Hack to get around hmmer2go's convention of stripping away the last 'extension' on output files.

    my ( $split_dir ) = @_;

    my @files = glob("$split_dir/split_fasta.*");

    for my $orig_file ( @files ) {

        my $extension = ( fileparse( $orig_file, qr/\.[^.]*/) )[2];
        my $new_name = $orig_file.$extension;
        move( $orig_file, $new_name );

    }

}


sub launch_grid_job {
# Given a shell script, launch it via qsub.

    my ( $shell_script, $outfile, $errfile, $queue, $job_array_max, $grid_work_dir ) = @_;

    my $std_error  = "$working_dir/$errfile";
    my $std_out    = "$working_dir/$outfile";

    my $qsub_command = "qsub -P $project_code -e $std_error -o $std_out";
    $qsub_command .= " -l $queue" if $queue;
    $qsub_command .= " -t 1-$job_array_max" if $job_array_max;
    $qsub_command .= " -wd $grid_work_dir";

    $qsub_command .= " $shell_script";

    my $response = `$qsub_command`;
    my $job_id;

    if ($response =~ (/Your job (\d+) \(.*\) has been submitted/) || $response =~ (/Your job-array (\d+)\./)) {

        $job_id = $1;

    } else {
        die "Problem submitting the job!: $response";
    }

    return $job_id;

}
 

sub wait_for_grid_jobs_arrays {
# given an array of job_ids, wait until they are all done.

    my ($job_ids,$min,$max) = @_;

    my $lch = build_task_hash_arrays( $job_ids,$min, $max );
    my $stats_hash = build_task_hash_arrays($job_ids, $min, $max);

    while ( keys %{$lch} ) {

        for my $job_id ( keys %{$lch} ) {

            my $response = `qacct -j $job_id 2>&1`;
            parse_response_arrays( $response, $lch,$stats_hash );
            sleep 1;

        }
    }
}


sub build_task_hash_arrays {

    my ($job_ids, $min_id, $max_id) = @_;

    my $lch;
    for my $job_id ( @{$job_ids} ) {

        for my $task_id ( $min_id .. $max_id ) {

            $lch->{ $job_id }->{ $task_id } = 0;

        }

    }

    return $lch;

}


sub parse_response_arrays {
# given a qacct response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $lch,$stats_hash ) = @_;
    return if ( $response =~ /error: job id \d+ not found/ );  # hasn't hit the grid yet.

    my @qacct_array = split ( /=+\n/, $response );
    @qacct_array = grep { /\S/ } @qacct_array; # get rid of empty record at beginning.

    for my $record ( @qacct_array ) {

        next if $record =~ /error: ignoring invalid entry in line/;

        chomp $record;

        my @rec_array = split ( "\n", $record );

        my %rec_hash;
        for my $line (@rec_array) {

            $line =~ s/(.*\S)\s+$/$1/;
            my ( $key, $value ) = split ( /\s+/, $line, 2 );
            $rec_hash{ $key } = $value;

        }

        if ( defined $rec_hash{taskid} && defined $rec_hash{jobnumber} ) {

            my ($task_id, $job_id) = @rec_hash{'taskid','jobnumber'};

            unless ( $stats_hash->{ $job_id }->{ $task_id } ) {

                $stats_hash->{ $job_id }->{ $task_id } = \%rec_hash;

                # clear the task from the lch
                delete $lch->{ $job_id }->{ $task_id };

                # clear the job if all tasks are accounted for
                delete ( $lch->{ $job_id } ) unless ( keys %{ $lch->{ $job_id } } );

                print "Found task $task_id from job $job_id\n" if ($opts{debug});

            }

        } else {

            print "Problem with one of the jobs' qacct info.\n";

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

    # duplicate $SGE_TASK_ID string is intentional and necessary, and unintentionally hilarious
    my @cmd = ( @$hmmer2go_cmd_ref, '-i', "$split_dir".'/split_fasta.$SGE_TASK_ID.$SGE_TASK_ID' );
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

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'mapterms', '-i', $tblout_file, '-p', 'pfam2go', '-o', $map_file, '--map' );

        system( @cmd ) == 0 || die ( "Problem running hmmer2go 'mapterms': ", join( ' ', @cmd ), "\n");

        until ( -s $term_mapping_tsv ) { print "Waiting for $term_mapping_tsv to magically finish getting written.\n"; sleep 10; }

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

        my @cmd = ( '/usr/bin/env', 'perl', '-I', $HMMER2GO_LIBDIR, $HMMER2GO_EXEC, 'mapterms', '-i', $tblout_file, '-p', "$HMMER2GO_DATA/tigrfams2go", '-o', $map_file, '--map' );

        system( @cmd ) == 0 || die ( "Problem running hmmer2go 'mapterms': ", join( ' ', @cmd ), "\n");

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

#        chomp;
#        my ( $id, $go_id ) = ( split( "\t", $_ ) )[ 0, 4 ];

#        push @{$input_ids{ $id }}, $go_id;
    
#    }

    ( my $input_base  = $opts{ input_fasta } ) =~ s/([^\.+])\..*/$1/;
    $input_base = (fileparse($input_base))[0];
    
    my $output_file = "$working_dir/$input_base.cluster_roles.txt";
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

        $opts{ hmm_db } = $opts{ hmm_db } // 'both';
        $opts{ hmm_db } = lc $opts{ hmm_db };
        $opts{ hmm_db } = 'tigrfams' if ( $opts{ hmm_db } eq 'tigrfam' ); # Cuz people gonna do it.

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

    die $errors if $errors;

}
