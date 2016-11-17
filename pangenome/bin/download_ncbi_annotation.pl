#!/usr/local/bin/perl -w

use warnings;
use strict;
$|++;

=head1 NAME

download_ncbi_annotation.pl

=head1 SYNOPSIS

    USAGE: ./download_ncbi_annotation.pl -a <file of accessions> [-o output/dir]

=head1 OPTIONS

B<--accessions, a>    :   File of accessions to download

B<--output, o>        :   Output directory. [Default: Current working dir]

B<--help, -h>         :   Displays Help

=head1 DESCRIPTION

This script downloads from NCBI through their Eutils interface genbank files for a list
of given accessions.

=head1 INPUTS

This script takes a file listing the accessions that should be downloaded.
It can accept WGS & GenBank Accessions or BioSample IDs (and a combination of all). 
Along with an accession please enter a 'name' you'd like to use for that Genome.

If a genome has multiple GenBank accessions associated with it, seperate them by commas (as seen below in the example)

File Example:
ANGE00000000          ANGE00000000
SAMN02141705          SAMN02141705
Klebsiella pneumoniae 30684/NJST258_2          CP006918.1,CP006921.1,CP006922.1,CP006919.1          

**Note: You can just list an genus name and it will pull all associated completed chromosomes/plasmids as well as WGS 
annotation

File Example:
genome   Enterobacter
          
=head1 OUTPUTS

This script outputs a fasta file labeled <name>.gb . The name being what was supplied in the accessions file passed in.
It also creates a combined list file of all gb files(genbank_acc.list), this file can then be passed into the parsed_genbank_file script to 
create files that can be fed into the Pangenome.

=head1 CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use Cwd;
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use Pangenome::NCBIDownload;

my $MAX_RETRIES = 3;

my %opts;
GetOptions( \%opts, 
        'accessions|a=s',
        'output|o=s',
        'help|h') || die "Error getting options! $!";
pod2usage( {-exitval => 0, -verbose => 2} ) if $opts{ help };

my $output = check_params();

my $dobj = Pangenome::NCBIDownload->new($output);

my @accessions = read_file($opts{accessions});

open(my $list_fh, ">", "$output/genbank_acc.list");
open(my $failed_fh, ">", "$output/look_closer.list");
open(my $incomplete_fh, '>', "$output/incomplete_CDS_counts.list");

# For each passed in accession download Fasta file from NCBI
foreach my $line (@accessions){

    $line =~ s/\s+$//;

    my($name,$acc) = split(/\t/,$line);

    unless($name && $acc){
        die("ERROR: accessions file must have the a genome name and accession. ie: <name><tab><accession>\n");
    }

    $dobj->set_accession_name($name);
    $dobj->set_accession($acc);
    
    my $tries = 1;

    if ($name eq 'genome') {

        my $accessions = $dobj->download_gb_genomes();

        foreach my $name (keys %$accessions){
            
            #only care about genomes that match given term
            if ($name =~ /$acc/) {

                print "$name\n";
                
                my @accs = keys %{$accessions->{$name}};
                my $string = join(",",@accs);
                
                print "$string\n";
                
                $dobj->set_accession_name($name);
                $dobj->set_accession($string);
                
                do {

                    $dobj->download_files("gb");
                    $tries++;

                } until ( ! is_corrupt( $dobj->get_accession_name . '.gb' ) || $tries == $MAX_RETRIES );

                if ( -s $dobj->get_accession_name .'.gb' ) {
                    if ( is_corrupt( $dobj->get_accession_name.'.gb' ) ) {
                        print $failed_fh "$name\t$acc\n";
                    } else {
                        print $list_fh "$name.gb\n";
                    }  
                } else {   
                    print $failed_fh "$name\t$acc\n";
                }
            
            }

        }
    
    } else {

        do {    

            $dobj->download_files("gb");
            $tries++;
            sleep 1;

        } until ( ! is_corrupt( $output . '/' .$dobj->get_accession_name . '.gb' ) || $tries == $MAX_RETRIES );

        my $path = $output . '/' . $dobj->get_accession_name . '.gb';

        if ( -s $path ) {
            if ( is_corrupt( $path ) ) {
                print $failed_fh "$name\t$acc\n";
            } elsif ( is_incomplete( $path ) ) {
                print $incomplete_fh "$name\t$acc\n";
            } else {
                print $list_fh "$name.gb\n";
            }
        } else {
            print $failed_fh "$name\t$acc\n";
        }

    }

}

close $failed_fh;
close $list_fh;
close $incomplete_fh;

exit(0);


sub is_incomplete {

    my ( $file ) = @_;

    my $ret_value = 0;

    # Check for completeness.  First grab expected number of CDS
    my $expected_cds = get_expected_cds( $file );

    # now look inside for the exact counts:
    my $found_CDS = count_cds( $file );
    if ( $found_CDS ) {

        # If we don't have an expected cds, print what we have and move on.
        if ( $expected_cds == 0 ) {
            warn "Don't know how many CDS to look for in $file, but found $found_CDS\n";
        } elsif ( $found_CDS != $expected_cds ) {
            $ret_value++;
            warn "Found $found_CDS CDS but expected $expected_cds CDS from $file!\n" 
        }
    } else {
        $ret_value++;
        warn "Found ZERO CDS features in $file!\n";
    } 

    return $ret_value;

}


sub is_corrupt {

    my ( $file ) = @_;
    my $ret_value = undef;

    if ( -s $file ) {

        # check for '//' on second to last line of genbank file
        my $sysout = `tail -n 2 $file | head -n 1`;

        if ( $sysout !~ m{//} ) {
            $ret_value++;
            warn "$file is missing the trailing // for last record.\n";
        }

        # check for yucky strings like 'Resource temporarily unavailable' in the file
        if ( `grep -m 1 -n 'Resource temporarily unavailable' $file` ) {
            $ret_value++;
            warn "$file contains the string 'Resource temporarily unavailable'\n";
        }


    } else {

        warn "$file does not exist/is empty!\n";
        $ret_value++;

    }

    return $ret_value;

}


sub get_expected_cds {

    my ( $file ) = @_;
    my $expected_cds_count = 0;

    open( my $cds_fh, '<', $file ) || die "Error opening $file!\n$!\n";

    while ( <$cds_fh> ) {
#        if ( /CDS\s+([\d\,]+)/ ) {
        if ( / {12}CDS\s+[:]{0,2}\s+([\d\,]+)/ ) {
            ( $expected_cds_count = $1 ) =~ s/,//g; # don't forget to strip out the commas.
        } elsif ( /Pseudo Genes\s+::\s+([\d\,]+)/ ) { 
            ( my $pseudos = $1 ) =~ s/,//g;
            $expected_cds_count += $pseudos;
            last;
        }
    }

    return $expected_cds_count;

}


sub count_cds {

    my ( $file ) = @_;

    my $cds_count = 0;

    open( my $cds_fh, '<', $file ) || die "Error opening $file!\n$!\n";

    while ( <$cds_fh> ) {
        $cds_count++ if /^ {5}CDS/;
    }

    return $cds_count;

}


sub check_params{

    my $output = $opts{output} // cwd;
    my $errors;

    mkdir($output) unless (-s $output);
   
    if ($opts{accessions}) {
        $errors .= "File does not exist or is size zero: $opts{accessions}\n" unless(-s $opts{accessions});
    } else {
        $errors .= "Must provide parameter: --accessions\n" unless($opts{accessions});
    }

    $errors .= "\nUsage: ./download_nbci_annotation.pl -a accessions.txt [-o output/dir]\n" if $errors;

    die ($errors) if $errors;
 
    return($output);

}

