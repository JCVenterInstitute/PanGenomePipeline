#!/usr/local/bin/perl
use warnings;
use strict;
$|++;

=head1 SYNOPSIS

check_att_files.pl - pre-processor to ensure reliability of att_files for pangenome runs

=head1 USAGE

    check_att_files.pl --db_list <file path> --att_dir <directory> [ --fasta_dir <directory> --resolve ]

=head1 OPTIONS

=over

=item B<--db_list, -d>  :   List of genomes whose files will be examined

=item B<--att_dir, -a>  :   Path to a directory containing .att files for the given genomes

=item B<--att_suffix>   :   Provide an alternate suffix for the .att file [default is ".att"]

=item B<--resolve>      :   Attempt to make identical thigns different in att files and fasta files [ Optional, requires B<--fasta_dir> ]

=item B<--fasta_dir, -f>    :   Path to a directory containing .pep or .nuc files for the given genomes

=item B<--nuc>          :   Look for .nuc files instead of .pep files in --fasta_dir when using --resolve

=back

=head1 DESCRIPTION

Examines .att files and .pep (and/or .nuc) files in the given directories, building a list of seen loci.  When a locus prefix duplication is found, creates a report of the conflicting loci.  When B<--resolve> is used, will also attempt to make those conflicting loci different, appending lowercase letters in sequence from a, b, c... etc. to the end of the duplicated prefix in the affected files.  

This script can also determine when a locus id appears more than once within a particular genome.  When B<--resolve> is used, the script will remove all features with this locus except for the longest feature that appears first.  (For example, if four features share locus FOO12345 and have lengths in the order: 5000, 2500, 350, 5000, only the first feature 5000 long will be kept, the others will be removed.)

=head1 OUTPUT

Ouptut goes to stdout, so if the user wishes to preserve the report of conflicts, they should redirect it into a file as per their shell.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Bio::SeqIO;
use File::Copy qw( mv );
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use Pod::Usage;

my %opts;
GetOptions( \%opts,
            'db_list|d=s',
            'att_dir|a=s',
            'att_suffix=s',
            'fasta_dir|f=s',
            'resolve',
            'nuc',
            'help|h',
           ) || die "Can't get options! $!\n";
pod2usage( { verbose => 2, exitval => 0 } )  if ( $opts{ help } );

my $att_suffix = '.att';

check_params( \%opts );

open( my $dlh, '<', $opts{ db_list } ) || die "Can't open dblist: $!";

my $loci = {};  # Store per-genome dupes
my $locus_prefixes = {}; # look for duplicate locus tag prefixes between genomes
my $file_count = 0;
while ( <$dlh> ) {

    $file_count++;

    chomp;

    my $att_file_path = $opts{ att_dir } . '/' . $_ . $att_suffix;
    unless ( -s $att_file_path ) {

        print "Skipping this att_file that doesn't exist or has zero-size: $att_file_path\n";
        next;

    }

    open( my $afh, '<', $att_file_path ) || die "Can't open $att_file_path: $!";

    print "Reading $file_count:  $att_file_path\n";

    while( <$afh> ) {

        chomp;
        my ( $molecule, $locus_tag, $start, $end, $description, $genome ) = split( "\t" );
        # Nothing scary, just assiging into a hash slice:
        my %locus;
        @locus{ qw(molecule locus_tag start end description genome) } = 
                ( $molecule, $locus_tag, $start, $end, $description, $genome );

        my $index = $genome . '_' . ( $start < $end ) ? ( $start . '_' . $end ) : ( $end . '_' . $start );

        my $locus_prefix = get_prefix( $locus_tag );
        $locus_prefixes->{ $locus_prefix }->{ $genome }++;

        $loci->{ $genome }->{ $locus_tag }++;

    }

}

for my $genome ( keys %$loci ) {

    my @dupes = grep { $loci->{ $genome }->{ $_ } > 1 } keys %{$loci->{ $genome }};

    if ( scalar @dupes ) {

        print "The following loci were seen more than once in $genome:\n";
        print map { "$_ : $loci->{ $genome }->{ $_ }\n" } @dupes;

        if ( $opts{ resolve } ) {
            print "Resolving duplicate loci in $genome\n";
            resolve_dupes( $genome, \@dupes );
        }

    }

}

my $dupe_prefixes;
for my $prefix ( keys %$locus_prefixes ) {

    my $genome_count = scalar keys %{$locus_prefixes->{ $prefix }};

    if ( $genome_count > 1 ) {

        $dupe_prefixes->{ $prefix } = $locus_prefixes->{ $prefix };

        print "$prefix is used by $genome_count genomes: ";
        print join(",", sort { $a cmp $b }  keys %{$locus_prefixes->{ $prefix }});
        print "\n";

    }

}

resolve_prefixes( $dupe_prefixes )  if ( $opts{ resolve } );

exit(0);


sub resolve_dupes {
# Will only keep the longest of all loci using the same locus id!

    my ( $genome, $dupes ) = @_;
    my $att_file = $opts{ att_dir } . '/' . $genome . $att_suffix;

    # find the lines in the att file that match this.
    for my $locus ( @$dupes ) {

        my @lines = `grep $locus $att_file`;
        map { chomp } @lines;

        # get the max length of these things
        my $max_len = 0;
        for my $line ( @lines ) {

            my ( $start, $end ) = ( split( "\t",$line ) )[2,3];
            my $len = ( $end > $start ) ? $end - $start + 1 : $start - $end + 1;
            $max_len = ( $len > $max_len ) ? $len : $max_len;

        }

        # we want to keep the first entry in both the att and fasta files that have the max_len
        # and get rid of all the rest.
        strip_from_fasta_file( $genome, $locus, $max_len );
        strip_from_att_file( $genome, $locus, $max_len );

    }


}


sub strip_from_att_file {
# Remove all but the longest features from the given genome's att_file
# Leave only the first feature that has the same length as the max length

    my ( $genome, $locus, $max_len ) = @_;

    my $att_file = $opts{ att_dir } . '/' . $genome . $att_suffix;
    my $tmp_file = $att_file . '.tmp';

    open( my $afh, '<', $att_file ) || die "Can't open $att_file for reading: $!\n";
    open( my $tfh, '>', $tmp_file ) || die "Can't open $tmp_file for writing: $!\n";

    my $have_kept = 0;

    while ( <$afh> ) {

        my $line = $_;

        if ( /$locus/ ) {
            next if $have_kept;
            my ( $start, $end ) = ( split( "\t", $line ) )[2,3];
            my $len = ( $end > $start ) ? $end - $start + 1 : $start - $end + 1;
            if ( $len == $max_len ) {
                $have_kept = 1;
            } else {
                next;
            }
        }

        print $tfh $line;

    }

    mv $tmp_file, $att_file;
    
}


sub strip_from_fasta_file {
# remove all sequences that share locus ids except
# for the first sequence that is max_len

    my ( $genome, $locus, $max_len ) = @_;

    my $fasta_file = $opts{ fasta_dir } . '/' . $genome;
    $fasta_file .= ( $opts{ nuc } ) ? '.nuc' : '.pep';
    my $tmp_file = $fasta_file . '.tmp';

    my $seq_out = Bio::SeqIO->new( -file => ">$tmp_file", '-format' => 'Fasta');
    my $seq_in = Bio::SeqIO->new( -file => $fasta_file, '-format' => 'Fasta');

    my $have_kept = 0;

    while ( my $seq = $seq_in->next_seq ) {

        my ( $id, $sequence, $len );
        $id = $seq->id();
        $sequence = $seq->seq;
        $len = $seq->length;

        if ( $id =~ /$locus/ ) {

            next if $have_kept;

            if ( $len == $max_len ) {
                $have_kept = 1;
            } else {
                next;
            }

        }

        $seq_out->write_seq( $seq );

    }
 
    mv $tmp_file, $fasta_file;

}


sub run_system {

    my ( $cmd ) = @_;

#    print "Running $cmd\n";

    system( $cmd );

    if ($? == -1) {
        print "failed to execute: $!\n";
    } elsif ($? & 127) {
        printf "child died with signal %d, %s coredump\n",
        ($? & 127),  ($? & 128) ? 'with' : 'without';
    } else {
        my $ret_val = $? >> 8;
        unless ( $ret_val == 0 ) {
            printf "child exited with value %d\n", $ret_val;
        }
    }

}


sub resolve_prefixes { 

    my ( $dupes ) = @_;

    if ( ! defined $dupes ) {
        print "--resolve used, but no duplicate prefixes to fix.\n";
        return undef;
    }

    for my $prefix ( keys %$dupes ) {

        my $prefix_suffix = 'a';

        for my $genome ( sort keys %{$dupes->{ $prefix }} ) {

            my $new_prefix = $prefix . $prefix_suffix;
            $prefix_suffix++;
            my $att_file = "$opts{ att_dir }/$genome".$att_suffix;
            my $fasta_file = "$opts{ fasta_dir }/$genome.";
            $fasta_file .= ( defined $opts{ nuc } ) ? 'nuc' : 'pep';

            print "Fixing $genome, changing $prefix to $new_prefix\n";

            for my $file ( $att_file, $fasta_file ) {
                run_system( "perl -pi.bak -e  's/$prefix/$new_prefix/' $file" );
            }

        }

    }


}


sub get_prefix {

    my ( $locus ) = @_;

    my $prefix = undef;

    if ( $locus =~ /([^_]+)_/ ) {

        $prefix = $1;

    } elsif ( $locus =~ /^(\D+)\d+/ ) {

        $prefix = $1;

    } else {

        die "Yo, handle me: $locus\n";

    }

    return $prefix;

}


sub check_params {

    my ( $opts ) = @_;

    my $errors = '';

    if ( ! $opts->{ db_list } ) {
        $errors .= "Please supply --db_list\n";
    } else {
        $errors .= "Please supply a non-empty --db_list\n" unless ( -s $opts->{ db_list } );
    }

    $att_suffix = $opts->{ att_suffix } // $att_suffix;

    $errors .= "Please supply --att_dir\n" unless ( $opts->{ att_dir } );

    $errors .= "Please supply --fasta_dir with --resolve\n" if ( $opts->{resolve} && ! $opts->{fasta_dir} );

    die $errors if $errors;

}
