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

get_patric_files.pl

=head1 SYNOPSIS

    USAGE: get_patric_files.pl  --genomes_list <list of genome names>
                                --file_source <source of files, Patric or RefSeq>
                                --file_type <type of file to download: cds, faa or all>
         
=head1 OPTIONS

B<--genomes_list, g>       :  File of genome directory names to download from Patric FTP site (ie.Acetobacter_aceti/)
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

GetOptions(\%opts, 'genomes_list|g=s',
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

    if($opts{genomes_list}){
	$errors .= "$opts{genomes_list} does not exist or is size zero\n" unless (-s $opts{genomes_list});
	$patric_file = $opts{genomes_list};
    }else{
	$errors .= "Must provide list of genomes, --genomes_list\n";
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
