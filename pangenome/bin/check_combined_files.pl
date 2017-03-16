#!/usr/bin/env perl
use warnings;
use strict;
$|++;

=head1 NAME

check_combined_files.pl - utility for validating combined.fasta and combined.att

=head1 SYNOPSIS

    USAGE: check_combined_files.pl --att_file <combined.att> --fasta_file <combined.fasta>

=head1 OPTIONS

B<--att_file, -a>   :   Path to the combined.att file

B<--fasta_file, -f> :   Path to the combined.fasta file

B<--log, -l>        :   Path to log file [DEFAULT: ./check_combined_files.log]

B<--help, -h>       :   Print help info

=head1 DESCRIPTION

Compare corresponding combined.fasta and combined.att files for a pangenome (panoct) run.

Checks include:
duplicate ids in both files
1-to-1 relationship between feature ids in the .att and .fasta file

If any are found, this script will exit with a non-zero value.  Specifically, the total number of
ids associated with any kind of error.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my $DEFAULT_ATT     = './combined.att';
my $DEFAULT_FASTA   = './combined.fasta';
my $DEFAULT_LOGFILE = './check_combined_files.log';

my %opts;
GetOptions( \%opts,
            'att_file|a=s',
            'fasta_file|f=s',
            'log|l=s',
            'help|h',
            ) || die "Can't retrieve options: $!\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{ help };

#set files.
my $att_file    = $opts{ att_file }   // $DEFAULT_ATT;
my $fasta_file  = $opts{ fasta_file } // $DEFAULT_FASTA;
my $log_file    = $opts{ log }   // $DEFAULT_LOGFILE;

open( my $lfh, '>', $log_file ) || die "Can't open log file $log_file: $!\n";

my %fasta_ids;
my %att_ids;

my $fasta_duplicate_id_count = 0;
my @fasta_duplicate_ids = ();

my $att_duplicate_ids_count = 0;
my @att_duplicate_ids = ();

# Step ONE: Read in the fasta file headers.  Locus should be first string of non-space
# characters in the header.

_log("Parsing fasta file $fasta_file");

open( my $ffh, '<', $fasta_file ) || _die( "Can't open fasta file $fasta_file: $!");

while ( <$ffh> ) {

    next unless /^>/;

    if ( /^>(\S+)/ ) {
        $fasta_ids{ $1 }++;
    } else {
        _log( "Can't get id from line:\n$_" );
    }

}

# Should NOT have any duplicates!
@fasta_duplicate_ids = grep { $fasta_ids{ $_ } > 1 } keys %fasta_ids;
$fasta_duplicate_id_count = scalar @fasta_duplicate_ids;
if ( $fasta_duplicate_id_count ) {

    _log("The following $fasta_duplicate_id_count ids were seen more than once in $fasta_file");
    _log( join( "\n", @fasta_duplicate_ids ) );

}

_log("Parsing att file $att_file");

# Step TWO: Read in ids from the att file.
open( my $afh, '<', $att_file ) || _die( "Can't open $att_file: $!" );

while( <$afh> ) {

    my $locus_tag = (split(/\t/))[1];
    if ( $locus_tag ) {
        $att_ids{ $locus_tag }++
    } else {
        _log( "Can't get id from line:\n$_" );
    }

}

# Shouldn't have duplicates:
@att_duplicate_ids = grep { $att_ids{ $_ } > 1 } keys %att_ids;
$att_duplicate_ids_count = scalar @att_duplicate_ids;
if ( $att_duplicate_ids_count ) {

    _log("The following $att_duplicate_ids_count ids were seen more than once in $att_file:");
    _log( join( "\n", @att_duplicate_ids ) );

}

# Step THREE: Compare the lists. They should be identical.
# First, find ids from att not in fasta:
my @in_att_not_in_fasta;
my $in_att_not_in_fasta_count = 0;
for my $id ( sort keys %att_ids ) {

    unless ( exists $fasta_ids{ $id } ) {
        push @in_att_not_in_fasta, $id;
    }

}
$in_att_not_in_fasta_count = scalar @in_att_not_in_fasta;
if ( $in_att_not_in_fasta_count ) {

    _log("The following $in_att_not_in_fasta_count ids were seen in $att_file but not in $fasta_file:");
    _log( join( "\n", @in_att_not_in_fasta ) );

}

my @in_fasta_not_in_att;
my $in_fasta_not_in_att_count = 0;
for my $id ( sort keys %fasta_ids ) {

    unless ( exists $att_ids{ $id } ) {
        push @in_fasta_not_in_att, $id;
    }

}
$in_fasta_not_in_att_count = scalar @in_fasta_not_in_att;
if ( $in_fasta_not_in_att_count ) {

    _log("The following $in_fasta_not_in_att_count ids were seen in $fasta_file but not in $att_file:");
    _log( join( "\n", @in_fasta_not_in_att ) );

}

# Print a summary.
_log("\nSUMMARY:");
_log("Duplicates in att:\t\t$att_duplicate_ids_count");
_log("Duplicates in fasta:\t\t$fasta_duplicate_id_count");
_log("In att but not in fasta:\t$in_att_not_in_fasta_count");
_log("In fasta but not in att:\t$in_fasta_not_in_att_count");

my $exit_value = $att_duplicate_ids_count + $fasta_duplicate_id_count +
                 $in_att_not_in_fasta_count + $in_fasta_not_in_att_count;

exit( $exit_value );

sub _log {

    my ( $msg ) = @_;
    chomp $msg;
    print $lfh "$msg\n";

}

sub _die {

    my ( $msg ) = @_;
    chomp $msg;
    _log( $msg );
    die "$msg\n";

}
