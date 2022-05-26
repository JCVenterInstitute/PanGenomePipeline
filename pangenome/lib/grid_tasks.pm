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

=head1 NAME

grid_tasks.pm - not-so-fancy module for launching/monitoring grid jobs with a strong slant towards the JCVI grid configuration.

=head1 SYNOPSIS

    uae grid_tasks;

    my @job_ids;
    push @job_ids, launch_grid_job( $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue, $job_array_max );
    # ... repeat for more jobs if you want.

    wait_for_grid_jobs( \@job_ids );
    # Or if you launched as an array-job:
    wait_for_grid_jobs_arrays( \@job_ids );

=head1 DESCRIPTION

Quick start guide:

 1. Create a shell script for the job.
 2. Call launch_grid_job, store the returned job_id in an array
 3. Repeat if more jobs are to be launched.
 4. If invoked with $job_array_max defined, call wait_for_grid_jobs_arrays,
    otherwise, call wait_for_grid_jobs

 Control is returned upon grid job completion.

There are numerous use cases unsupported by this module at the moment.  Notably:

 When using jobs that wait for other jobs.
 Using task ids that increment with an integer other than 1.
 When resubmitting job_arrays that have errored-tasks reset by using qmod -cj
 
This module will likely grow and change to support fancier submissions.

=head1 Subprocedures

=over 4

=item C<launch_grid_job>

    For array jobs:
    my $job_id = launch_grid_job( $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue, $job_array_max );

    For single task jobs:
    my $job_id = launch_grid_job( $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue );

=cut

sub launch_grid_job {
# Given a shell script, launch it via qsub.

    my ( $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue, $job_array_max ) = @_;

    my $qsub_command = "qsub -l centos7 -V -P $project_code -o $stdoutdir -e $stderrdir -wd $working_dir";
    #$qsub_command .= " -l $queue" if $queue;
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

=item C<wait_for_grid_jobs>

    Given an array ref of job_ids, wait until they are all done.

    wait_for_grid_jobs( \@job_ids );

=cut

sub wait_for_grid_jobs {

    my $job_ids = shift;

    my $lch = build_task_hash( $job_ids );

    while ( keys %{$lch} ) {

        for my $job_id ( keys %{$lch} ) {

            my $response = `qacct -j $job_id 2>&1`;
            parse_response( $response, $lch, $job_id );
            sleep 1;

        }
    }

}


sub parse_response {
# NOT INTENDED TO BE CALLED DIRECTLY.
# Given a qacct response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $lch, $job_id ) = @_;
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


=item C<wait_for_grid_jobs_arrays>

    Given an array ref of array-job job_ids (Wow!), wait until they are all finished.

    wait_for_grid_jobs_arrays( \@job_ids, $min_task_id, $max_task_id, $debug);

    Note that $min_task_id and $max_task_id must be the same for all job_ids.
    (This should change.  Hopefully sooner rather than later.)

=back

=cut

sub wait_for_grid_jobs_arrays {
# given an array of job_ids, wait until they are all done.

    my ( $job_ids, $min, $max, $debug ) = @_;

    my $lch = build_task_hash_arrays( $job_ids, $min, $max );
    my $stats_hash = build_task_hash_arrays( $job_ids, $min, $max );

    while ( keys %{$lch} ) {

        for my $job_id ( keys %{$lch} ) {

            my $response = `qacct -j $job_id 2>&1`;
            parse_response_arrays( $response, $lch, $stats_hash, $job_id, $debug );
            sleep 1;

        }
    }

    #check_for_grid_errors( $job_id );  See note in subprocedure.

}


sub check_for_grid_errors {
# Given a job_id of a supposedly finished job_array, check that
# none of the qacct job sections have reported an error.

# NOTE This currently over-reports errors: Quite often when jobs appear on qacct
# they are filled with tasks that have error codes.  By the time someone can type
# qacct -j <jobid> on the command line, they've cleared theymselves up, but this
# script has a habit of catching them and throwing an error.
# Besides which, the original problem this was trying to catch is no longer an issue.
# Should the time arise in which this looks like a great feature to re-instate, I
# recoommend introducing a sleep statement, of oh, 30-60 seconds to give the grid time
# to sort itself.
#   JMI 20181023

    my ( $job_id ) = @_;

    my $response = `qacct -j $job_id 2>&1`;

    my @qacct_array = split( /=+\n/, $response );
    @qacct_array = grep { /\S/ } @qacct_array; # get rid of empty recoed at beginning;

    my %failures;

    for my $record ( @qacct_array ) {

        next if ( $record =~ /error: ignoring invalid entry in line/ );

        my $task_id;
        if ( $record =~ /taskid +(\d+)/ ) {
            $task_id = $1;
        } elsif ( $record =~ /taskid +undefined/ ) {
            $task_id = "$job_id (not a job-array)";
        }

        if ( $record =~ /failed +([^\n]+)\n/ ) {
            my $excuse = $1;
            $failures{ $task_id } = "'$excuse'" unless $excuse =~ /^0/;
        } 

    }

    if ( scalar keys %failures ) {

        my $message = "The following grid tasks (for job_id $job_id ) have failed.  The grid says the following reasons:\n";
        $message .= join("\n", map {"$_\t$failures{ $_ }"} keys %failures);

        die "$message\n";

    }

}


sub parse_response_arrays {
# NOT INTENDED TO BE CALLED DIRECTLY.
# Given a qacct response, delete a job id from the loop-control-hash when
# a satisfactory state is seen.

    my ( $response, $lch, $stats_hash, $job_id, $debug ) = @_;
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
# NOT INTENDED TO BE CALLED DIRECTLY.
# create a hash representing a single task.

    my $job_ids = shift;

    my $lch;
    for my $job_id ( @{$job_ids} ) {

        $lch->{ $job_id } = 0;

    }

    return $lch;

}


sub build_task_hash_arrays {
# NOT INTENDED TO BE CALLED DIRECTLY.
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
