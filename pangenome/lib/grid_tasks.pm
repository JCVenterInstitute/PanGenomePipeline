#!/usr/bin/env perl

## TODO: Create pod documentation for this module.

## quick how-to guide:
# 1. create a shell script for the job.
# 2. call launch_grid_job
# 3. if invoked with $job_array_max defined, call wait_for_grid_jobs_arrays,
#    otherwise, call wait_for_grid_jobs
# control is returned upon grid job completion.

sub launch_grid_job {
# Given a shell script, launch it via qsub.

    my ( $project_code, $working_dir, $shell_script, $outfile, $errfile, $queue, $job_array_max ) = @_;

    my $std_error  = "$working_dir/$errfile";
    my $std_out    = "$working_dir/$outfile";

    my $qsub_command = "qsub -P $project_code -e $std_error -o $std_out";
    $qsub_command .= " -l $queue" if $queue;
    $qsub_command .= " -t 1-$job_array_max" if $job_array_max;

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


sub wait_for_grid_jobs {
# given an array of job_ids, wait until they are all done.

    my $job_ids = shift;

    my $lch = build_task_hash( $job_ids );

    while ( keys %{$lch} ) {

        for my $job_id ( keys %{$lch} ) {

            my $response = `qacct -j $job_id 2>&1`;
            parse_response( $response, $lch );
            sleep 1;

        }
    }
}


sub parse_response {
# given a qacct response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $lch ) = @_;
    return if ( $response =~ /error: job id \d+ not found/ );  # hasn't hit the grid yet.

    my @qacct_array = split ( /=+\n/, $response );
    @qacct_array = grep { /\S/ } @qacct_array; # get rid of empty record at beginning.

    for my $record ( @qacct_array ) {

        next if ( $record =~ /error: ignoring invalid entry in line/ );

        chomp $record;

        my @rec_array = split ( "\n", $record );

        my %rec_hash;
        for my $line (@rec_array) {

            $line =~ s/(.*\S)\s+$/$1/;
            my ( $key, $value ) = split ( /\s+/, $line, 2 );
            $rec_hash{ $key } = $value;

        }

        if ( defined $rec_hash{ jobnumber } ) {

            my $job_id = $rec_hash{ 'jobnumber' };

            delete ( $lch->{ $job_id } )

        } else {

            print "Problem with one of the jobs' qacct info.\n";
        }

    }

}


sub wait_for_grid_jobs_arrays {
# given an array of job_ids, wait until they are all done.

    my ( $job_ids, $min, $max ) = @_;

    my $lch = build_task_hash_arrays( $job_ids, $min, $max );
    my $stats_hash = build_task_hash_arrays( $job_ids, $min, $max );

    while ( keys %{$lch} ) {

        for my $job_id ( keys %{$lch} ) {

            my $response = `qacct -j $job_id 2>&1`;
            parse_response_arrays( $response, $lch,$stats_hash );
            sleep 1;

        }
    }
}

sub parse_response_arrays {
# given a qacct response, delete a job id from the loop-control-hash when
# a satisfactory state is seen.

    my ( $response, $lch, $stats_hash, $debug ) = @_;
    return if ( $response =~ /error: job id \d+ not found/ );  # hasn't hit the grid yet.

    my @qacct_array = split ( /=+\n/, $response );
    @qacct_array = grep { /\S/ } @qacct_array; # get rid of empty record at beginning.

    for my $record ( @qacct_array ) {

        next if ( $record =~ /error: ignoring invalid entry in line/ );

        chomp $record;

        my @rec_array = split ( "\n", $record );

        my %rec_hash;
        for my $line ( @rec_array ) {

            $line =~ s/(.*\S)\s+$/$1/;
            my ( $key, $value ) = split ( /\s+/, $line, 2 );
            $rec_hash{ $key } = $value;

        }

        if ( defined $rec_hash{ taskid } && defined $rec_hash{ jobnumber } ) {

            my ( $task_id, $job_id ) = @rec_hash{ 'taskid', 'jobnumber' };

            unless ( $stats_hash->{ $job_id }->{ $task_id } ) {

                $stats_hash->{ $job_id }->{ $task_id } = \%rec_hash;

                # clear the task from the lch
                delete $lch->{ $job_id }->{ $task_id };

                # clear the job if all tasks are accounted for
                delete ( $lch->{ $job_id } ) unless ( keys %{ $lch->{ $job_id } } );

                print "Found task $task_id from job $job_id\n" if ( $debug );

            }

        } else {

            print "Problem with one of the jobs' qacct info.\n";

        }

    }

}


sub build_task_hash {
# create a hash representing a single task.

    my $job_ids = shift;

    my $lch;
    for my $job_id ( @{$job_ids} ) {

        $lch->{ $job_id } = 0;

    }

    return $lch;

}


sub build_task_hash_arrays {
# create an array of hashes representing the tasks in a job

    my ( $job_ids, $min_id, $max_id ) = @_;

    my $lch;
    for my $job_id ( @{$job_ids} ) {

        for my $task_id ( $min_id .. $max_id ) {

            $lch->{ $job_id }->{ $task_id } = 0;

        }

    }

    return $lch;

}


1;
