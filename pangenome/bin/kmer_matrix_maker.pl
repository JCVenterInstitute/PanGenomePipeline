#!/usr/local/bin/perl
use warnings;
use strict;
$|++;

=head1 NAME

kmer_matrix_maker.pl - create kmer count matrix given a list/directory of genome fasta files.

=head1 SYNOPSIS

    USAGE: kmer_matrix_maker.pl [ --fasta_dir | --fasta_list ] <path to input>
                                --output <path to output file>

=head1 OPTIONS

B<--fasta_dir,-d>   :   Path to a directory containing the set of input genome fasta files  
                        ALL files within the directory will be considered input genome fasta files.

B<--fasta_list,-l>  :   Path to a file containing paths to each input genome fasta file.  Meryl files
                        will be created alongside them, and removed once meryl is done with them.

B<--output,-o>      :   Path to desired output file location.  REQUIRED.

=head1 DESCRIPTION

Will run meryl on each input genome fasta file and build a matrix file containing genomes
as rows and kmers on columns.  Values represent the number of times each kmer (specifically,
6-mers) are found in a genome.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use Pod::Usage;
use File::Basename;
use File::Path;
use File::Temp qw( tempfile );;

my $MERYL_BIN = '/usr/local/bin/meryl';

my %opts;
GetOptions( \%opts,
            'fasta_dir|d=s',
            'fasta_list|l=s',
            'output|o=s',
            'help|h',
        ) || die "Problem getting options.\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{help};


check_options( \%opts );

# Build list of inputs
my $input_files = build_input_list( \%opts );

my %kmer_hash;
my %kmer_names;
my @filenames;

my $count = 0;
my $total = scalar( @$input_files );

for my $file ( @{$input_files} ) {

    $count++;
    print STDOUT "Working on $file; $count of $total\n";

    my $small_name = basename( $file );
    push @filenames, $small_name;

    # run meryl and collect results
    ## First, build the table. 
    my $temp_file = '';
    ( undef, $temp_file ) = tempfile(); 
    my $cmd = "$MERYL_BIN -B -C -m 6 -s $file -o $temp_file";
    system( $cmd ) && die "Trouble running meryl.\n";

    ## Read these kmers into our hash.
    parse_kmer_file( \%kmer_hash, $small_name, $temp_file );

    ## Cleanup meryl's leftovers or stuff will fill up fast-like.
    my @trash = ( "$file.fastaidx", "$temp_file.mcdat", "$temp_file.mcidx", $temp_file );
    unlink for @trash;

}

# Print the matrix.
open( my $ofh, '>', $opts{output} ) || die "Bummer... all that work and can't open $opts{output} to save it because: $!\n";
select $ofh;

@filenames = sort {$a cmp $b} @filenames;
## Print the header (kmer names)
my @kmer_names = sort {$a cmp $b} keys %kmer_names;
print "\t", join( "\t", @kmer_names ), "\n";
## Print the rest of the rows
for my $file ( @filenames ) {

    my @row; 
    print "$file\t"; # row "header"

    # Wanted to do this with a hash slice, but have to accept the possibility (rare as it may be)
    # that a particular 6mer may actually not be present in a given input sequence.  So we'll loop:
    for my $kmer ( @kmer_names ) {
        my $bit = ( exists $kmer_hash{$file}->{$kmer} ) ? $kmer_hash{$file}->{$kmer} : 0;
        push @row, $bit;
    }
    print join( "\t", @row ), "\n";
   

}

exit(0);


sub parse_kmer_file {

    my ( $kmer_hash, $input_file, $kmer_file ) = @_;

    # Read a dump of the mers
    my $cmd = "$MERYL_BIN -Dt -s $kmer_file";

    my $count = '';
    for (`$cmd`) {
        chomp;
        if ( /^>(\d+)$/ ) {
            $count = $1;
        } else {
            $kmer_hash->{ $input_file }->{ $_ } = $count;
            $kmer_names{ $_ } = ();
        }
    }

}


sub build_input_list {

    my ( $opts ) = @_;

    my @input_files;

    if ( $opts->{fasta_dir} ) {

        chdir( $opts->{fasta_dir} ) || die "Can't chdir to $opts->{fasta_dir}\n";
        opendir( FASTADIR, $opts->{fasta_dir} ) || die "Can't opendir $opts->{fasta_dir}: $!\n";
        @input_files = grep !/^\.\.?$/, readdir( FASTADIR );
        closedir FASTADIR;

    } else {

        open( my $ifh, '<', $opts->{fasta_list} ) || die "Can't open $opts->{fasta_list}: $!\n";
        while ( <$ifh> ) {
            chomp;
            push @input_files, $_;
        }

    }

    return \@input_files;

}


sub check_options {

    my ( $opts ) = @_;

    my $errors = '';

    if ( $opts->{fasta_dir} && $opts->{fasta_list} ) {
        $errors .= "Please specify EITHER --fasta_dir OR --fasta_list\n";
    }

    if ( ( ! $opts->{fasta_dir} ) && ( ! $opts->{fasta_list} ) ) {
        $errors .= "Please specify --fasta_dir OR --fasta_list\n";
    }

    if ( ! $opts->{output} ) {
        $errors .= "Please specify --output\n";
    }

    die $errors if $errors;

}
