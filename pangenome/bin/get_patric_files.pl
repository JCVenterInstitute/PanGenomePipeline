#!/usr/local/bin/perl -w

use warnings;
use strict;
$|++;

=head1 NAME

get_patric_files.pl

=head1 SYNOPSIS

    USAGE: get_patric_files.pl  --db_list <list of genome names>
                                --file_source <source of files, Patric or RefSeq>
                                --file_type <type of file to download: cds, faa or all>
         
=head1 OPTIONS

B<--db_list, d>       :  File of genome directory names to download from Patric FTP site (ie.Acetobacter_aceti/)
                         http://brcdownloads.patricbrc.org/patric2/genomes/

B<--file_source, s>   :  Specifies to download the Patric or RefSeq annotation (options: Patric or RefSeq, not case sensitive)

B<--file_type, t>     :  Specifies what type of file to download, cds or faa (all for both). Not case sensitive, default all. 

B<--help,-h>          :  Display this help message.

=cut 
 
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Cwd;
use Data::Dumper;
use File::Slurp;
use File::Path;

my %opts;
my($patric_file,$file_source,$file_type);

GetOptions(\%opts, 'db_list|d=s',
	   'file_source|s=s',
	   'file_type|t=s',
	   'help|h');
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

&check_params;

my $patric_ftp_url = 'http://brcdownloads.patricbrc.org/patric2/genomes/';

my @file_lines = read_file($patric_file);

foreach my $line(@file_lines){
    $line =~ s/\s+//;
    $line =~ s/\/$//;
    my $path = $patric_ftp_url . "$line/$line";
    my $source = ($file_source eq 'PATRIC') ? "PATRIC" : "RefSeq";

    system("wget $path" .".$source.cds.tab") if $file_type =~ /(CDS|ALL)/;
    system("wget $path" .".$source.faa") if $file_type =~ /(FAA|ALL)/;
}

exit(0);

sub check_params{
    my $errors;

    if($opts{db_list}){
	$errors .= "$opts{db_list} does not exist or is size zero\n" unless (-s $opts{db_list});
	$patric_file = $opts{db_list};
    }else{
	$errors .= "Must provide list of genomes, --db_list\n";
    }

    if($opts{file_source}){
	$file_source = uc($opts{file_source});
	$errors = "$opts{file_source} must be either patric or refseq (not case sensative)\n" unless ($file_source =~ /(PATRIC|REFSEQ)/);
    }else{
	$errors .= "Must provide --file_source: patric or refseq (not case sensitive)\n";
    }

    if($opts{file_type}){
	$file_type = uc($opts{file_type});
	$errors = "$opts{file_type} must be cds,faa or all. (not case sensitive) \n" unless ($file_type =~ /(FAA|CDS|ALL)/);
    }else{
	$file_type = "ALL";
    }

    die($errors) if $errors;
    
}
