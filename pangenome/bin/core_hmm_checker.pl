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

core_hmm_checker.pl-  Script for running hmm3search against a list of fasta genomes
 
=head1 SYNOPSIS

  USAGE: core_hmm_checker.pl --hmm <hmm model file>
                             --fasta_list <list of peptide fasta file locations>
                             --cutoff <1-100>
                             --output <dir> [Optional]
                             --help [Optional]

=head1 OPTIONS

B<--hmm>           : Hmm model file containing HMMs of interest

B<--fasta_list|f>  : List of fasta file locations(pep files)

B<--cutoff|c>     : Percentage of HMMs from model that a genome must match [1-100]

B<--output|o>      : Output directory [Default: Current]

B<--no_cleanup>    : Keeps hmm files instead of removing

B<--help|h>        : Prints help

=head1 OUTPUT

hmm_matched.txt - Lists the files that matched the cutoff of interested HMMs
hmm_missing.txt - Lists the files that did not match the cutoff of interested HMMs

=cut
    
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Basename;
use File::Path;
use File::Copy;
use File::Touch;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin";


#Arguments/Options
my %opts;

GetOptions( \%opts, 'hmm=s',
	    'fasta_list|f=s',
	    'cutoff|c=i',
	    'no_cleanup',
	    'output|o=s',
	    'help|h') || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my $HMM_EXE = "/usr/local/bin/hmm3search";

#Check input params and set output directory
my ($CUTOFF,$OUTPUT) = &check_params;

my $core_accs = find_core_acc($opts{hmm});
my ($pass_hmm,$fail_hmm) = run_hmm_search($opts{hmm},$opts{fasta_list},$core_accs);
print_results($pass_hmm,$fail_hmm);

#Clean up tmp directory
unless($opts{no_cleanup}){
    rmtree("$OUTPUT/tmp/") or die "$!: Could not remove directory\n";
}

exit(0);

sub print_results{
    my ($pass,$fail) = @_;

    open(my $mfh, ">" , "$OUTPUT/hmm_matched.txt");
    open(my $fhf, ">" , "$OUTPUT/hmm_missed.txt");
      
    print $mfh "#genome\t% matched\n";
    print $fhf "#genome\t% matched\n";

    map{print $mfh "$_\n"} @$pass;
    map{print $fhf "$_\n"} @$fail;

  
    close $mfh;
    close $fhf;
}
sub find_core_acc{

    my $file = shift;
    my $count = `grep \"NAME\" $file`;
    my @core_hmms;
    
    my @values = split(/\s+/,$count);
    foreach my $value (@values){
	push(@core_hmms,$value) if($value ne "NAME");
    }

    return \@core_hmms;

}

sub run_hmm_search{

    my ($hmm,$fasta_list,$accs) = @_;
    my $tmp_dir = "$OUTPUT/tmp";
    mkdir($tmp_dir);

    my (@pass_hmm,@fail_hmm);
    
    open(my $fh , "<", $fasta_list);
    open(my $gfh, ">", "$OUTPUT/hmm_genome.list");
    
    while(<$fh>){
	my $line = $_;
	my ($genome,$location) = split(/\t/,$line);

	$genome =~ s/\s+$//;
	$location =~ s/\s+$//;
	
	my $tmp_file = "$tmp_dir" . "/$genome.tbl";
	
	my $cmd = $HMM_EXE . " --tblout $tmp_file --cut_tc $opts{hmm} $location";
	
	system($cmd) == 0 || die("ERROR: $cmd failed");

	my $hmm_count = check_hmm($tmp_file);
	
	my $perc = sprintf("%.3f",($hmm_count / scalar(@$accs))) * 100;
	
	$perc >= $CUTOFF ? push(@pass_hmm,"$genome\t$perc") : push(@fail_hmm,"$genome\t$perc");

	print $gfh "$genome\n" if($perc >= $CUTOFF);
		    
    }

    close $gfh;

    return(\@pass_hmm,\@fail_hmm);
}

sub check_hmm{

    my $file = shift;

    open(my $fh, "<", $file);
    my $hmms;
    
    while(<$fh>){

	my $line = $_;
	
	if($line !~ /^#/){
	    my @values = split(/\s+/,$line);
	    $hmms->{$values[2]}++;
	}
    }

    return (scalar keys %$hmms);
}

sub check_params{

    my ($errors,$cutoff,$output);
    
    unless($opts{hmm} || $opts{fasta_list} || $opts{hmm_acc}){
	die("Usage: ./core_hmm_checker.pl --hmm <core_hmm_library> --cutoff 95  --fasta_list <list_of_fasta_files\n");
    }

    $errors .= "$opts{hmm} : File does not exist or is size zero\n" unless (-s $opts{hmm});
    $errors .= "$opts{fasta_list} : File does not exist or is size zero\n" unless (-s $opts{fasta_list});

    $opts{output} ? $output = $opts{output} : $output = cwd;
    
    mkpath("$output/hmm") unless(-d "$output/hmm");
    $output .= "/hmm";
    
    if($opts{cutoff}){
	$cutoff = $opts{cutoff};
    }else{
	$errors .= "Must provide --cutoff [1-100]\n";
    }
    
    if($errors){
	print $errors;
	exit;
    }

    return($cutoff,$output);
}
