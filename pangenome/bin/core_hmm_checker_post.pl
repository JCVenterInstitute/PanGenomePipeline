#!/usr/bin/env perl
#Copy (C) 2016 The J. Craig Venter Institute (JCVI).  All rights reserved

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

use warnings;
use strict;
$|++;

=head1 NAME

core_hmm_checker_post.pl - Modifies Mapping and Download files based on HMM cutoffss

=head1 SYNOPSIS

  USAGE: core_hmm_checker_post.pl --mapping_file
                                  --download_file
                                  --matched_hmm
                                  --output <dir> [Optional]
                                  --help [Optional]

=head1 OPTIONS

B<--mapping_file|m>  : Mapping file created from FTP downloader

B<--download_file|d> : Download file created from FTP downloader

B<--matched_hmm|m>   : List of genomes that matched HMM

B<--output|o>      : Output directory [Default: Current]

B<--help|h>        : Prints help

=head1 OUTPUT

 Files created in hmm directory within the output directory

 hmm_mapping - Mapping file only containing genomes that matched the hmm cutoff
 hmm_download - Download file only containing genomes that matched the hmm cutoff
 hmm_genome_list - List of genomes names

=cut
    
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use File::Path;
use File::Copy;
use File::Touch;
use Cwd;

my %opts;

#GLOBAL VARS
my $OUTPUT;

#Arguments/Options
GetOptions( \%opts, 'mapping_file|m=s',
	    'download_file|d=s',
	    'matched_hmm|t=s',
	    'output|o=s',
	    'help|h');
	    
$OUTPUT = &check_params;

my $matched_hmm = parse_hmm_file($opts{matched_hmm});
create_new_mapping($matched_hmm,$opts{mapping_file});
create_new_download($matched_hmm,$opts{download_file});

sub check_params{

    my $errors;
    my $output;
    
    unless($opts{mapping_file} || $opts{download_file} || $opts{matched_hmm}){

	$errors .= "Usage: ./core_hmm_check_post.pl --mapping_file <file> --download_file <file> --matched_hmm <file> \n";

    }

    if($opts{OUTPUT}){
	$output = $OUTPUT;
    }else{
	$output = cwd;
	mkpath("$output/hmm") unless (-d "$output/hmm");
	$output .= "/hmm";
    }

    return $output;
}
sub create_new_download{
    my ($matches, $download) = @_;

    open(my $nfh, ">", "$OUTPUT/hmm_download.txt");
    open(my $fh, "<", $download);

    while(<$fh>){
	my @values = split(/\t/,$_);

	# Only print lines that have passed HMM cutoff
	# Adds % of HMM matched to col 7
	print $nfh $_ if (exists $matches->{$values[0]});
    }

    close $nfh;
    close $fh;
}

sub create_new_mapping{
    my ($matches, $mapping) = @_;

    open(my $nfh, ">", "$OUTPUT/hmm_mapping.txt");
    open(my $fh, "<", $mapping);

    while(<$fh>){
	my $line = $_;
	$line =~ s/\s+$//;
	
	my @values = split(/\t/,$line);
	
	push(@values,$matches->{$values[0]}) if (exists $matches->{$values[0]});

	print $nfh join("\t", @values);
    }

    close $nfh;
    close $fh;
}
sub parse_hmm_file{

    my $file = shift;
    open (my $fh, "<" , $file);

    my $hsh;
    
    while(<$fh>){
	my $line = $_;
	
	unless($line =~ /^#/){
	    my ($genome, $perc) = split(/\t/,$line);
	    $hsh->{$genome} = $perc;
	}       
    } 

    return $hsh;
}
