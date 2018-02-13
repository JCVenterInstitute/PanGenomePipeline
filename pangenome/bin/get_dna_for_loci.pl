#!/usr/bin/env perl
use warnings;
use strict;
$|++;

=head1 NAME

get_dna_for_loci.pl - retrieve genomic DNA for specified loci.

=head1 SYNOPSIS

    USAGE: get_dna_for_loci.pl -a <att_file> -f <genome_fasta> -l <loci_list_file> [ -o <output_file> ]

=head1 OPTIONS

=head2 Required Options

B<--att_file, -a>       :   Path to gene attribute file.

B<--genome_fasta, -f>   :   Path to genome fasta file.

B<--loci_list, -l>      :   Path to file containing one-column list of locus IDs.

=head2 Optional Options

B<--output_file, -o>    :   To send output to a file.  Otherwise, it will go to STDOUT.

=head1 DESCRIPTION

This script takes in a list of loci and an att file, and creates an instruction file to send to cutFasta on an also provided genome_fasta file.  (See cutFasta --help if you'd like to learn more about that)

If B<--output_file> is supplied, the resulting fasta will be written to that file, otherwise it will dump to STDOUT.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut


use Getopt::Long  qw(:config no_ignore_case no_auto_abbrev);
use Path::Tiny;
use Pod::Usage;

my $CUTFASTA_EXEC = '/usr/local/common/cutFasta';
my $instructions = './instructions';

my %opts;
GetOptions( \%opts,
    'att_file|a=s',
    'genome_fasta|f=s',
    'loci_list|l=s',
    'output_file|o=s',
    'help|h',
) || die "Problem getting options: $!\n";
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

check_options();

# Pre-load the att_file_hash
my %att;
my $afh = path( $opts{ att_file } )->filehandle( '<' );
while (<$afh>) {

    chomp;
    my ( $molecule_id, $locus, $end5, $end3, $description, $genome ) = split( '\t', $_ );
    @{$att{ $locus }}{ qw/molecule_id end5 end3/ } = ($molecule_id, $end5, $end3);

#use Data::Dumper; print Dumper (\%att); exit;
}

# get the coords for the loci from the gene attribute file.
# add a line to the instructions file.
my $lfh = path( $opts{ loci_list } )->filehandle( '<' );
my $ifh = path( $instructions )->filehandle( '>' );
while (<$lfh>) {

    chomp;
    my $locus = $_;
    my $line = "$locus\t$att{$locus}->{end5}\t$att{$locus}->{end3}\t$att{$locus}->{molecule_id}\n";
    print $ifh $line;

}

# Feed the instruction file to cutFasta and we're done!
my @cmd = ( $CUTFASTA_EXEC, '-i', $instructions );
push( @cmd, '-o', $opts{output_file} ) if $opts{ output_file };
push( @cmd, $opts{genome_fasta} );
system( @cmd ) && die "Problem running cutFasta:\n\t" . join( ' ', @cmd ) . "\n";

exit(0);


sub check_options {

    my $errors = '';

    if ( $opts{ att_file } ) {
        $errors .= "--att_file $opts{ att_file } does not exist or is empty\n" 
            unless ( -s $opts{ att_file } );
    } else {
        $errors .= "Please provide a --att_file\n";
    }

    if ( $opts{ genome_fasta } ) { 
        $errors .= "--genome_fasta $opts{ genome_fasta } does not exist or is empty\n"
            unless ( -s $opts{ genome_fasta } );
    } else {
        $errors .= "Please provide a --genome_fasta file\n";
    }

    if ( $opts{ loci_list } ) {
        $errors .= "--loci_list $opts{ loci_list } does not exist or is empty\n" 
            unless ( -s $opts{ loci_list } );
    } else {
        $errors .= "Please provide a --loci_list\n";
    }

    die $errors if $errors;

}
