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
#!/usr/bin/env perl

use warnings;
use strict;
$|++;

=head1 NAME

core_hmm_prep.pl-  Creates input files for core_hmm_checker.pl

=head1 SYNOPSIS

  USAGE: core_hmm_checker_prep.pl --gb_list <file>
                             --output <dir> [Optional]
                             --help [Optional]

=head1 OPTIONS

B<--gb_list|g>     : List of genome genbank files (<genome><location of gb file>)

B<--output|o>      : Output directory [Default: Current]

B<--help|h>        : Prints help

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use File::Path;
use File::Copy;
use File::Touch;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin";

my %opts;

#Arguments/Options
GetOptions( \%opts,
	    'gb_list|g=s',
	    'output|o=s',
	    'help|h');

my $OUTPUT = &check_params;

run_parse_gb($opts{gb_list});

sub run_parse_gb{

    my $gb_list = shift;
    my $parse_dir = "$OUTPUT/pep";

    mkdir($parse_dir) unless (-d $parse_dir);
    
    open(my $fh, "<", $gb_list);

    foreach (<$fh>){

	my ($genome,$location) = split(/\t/,$_);
        my $parse_file = "$parse_dir/parse.list";
	my $cmd = "cut -f 2 $opts{gb_list} > $parse_file";
	system($cmd) == 0 || die("Failed:$!");

	$cmd = "$Bin/parse_genbank_files.pl -l $parse_file -o $parse_dir --no_dos2unix";
	system($cmd) == 0 || die("Failed:$!");
    }

    make_genome_list_file($parse_dir);
}

sub make_genome_list_file{

    my $dir= shift;

    if(-s "$dir/genomes.list"){
	open(my $fh, "<" , "$dir/genomes.list");
	open(my $ofh, ">" , "$dir/pep.list");

	foreach (<$fh>){
	    my $genome = $_;
	    $genome =~ s/\s+$//;
	    print $ofh "$genome\t$dir/$genome" . ".pep\n";
	}

	close $fh;
	close $ofh;
    }else{
	die("ERROR: $dir/genomes.list file does not exist\n");
    }
}

sub check_params{

    my $output;
    my $errors;

    if($opts{gb_list}){
	$errors .= "File does not exist or is size zero\n" unless(-s $opts{gb_list});
    }else{
	$errors .= "Must supply --gb_list\n" unless($opts{gb_list});
    }
    
    if($opts{output}){
	mkdir($opts{output}) unless(-d $opts{output});
	$output = $opts{output};
    }else{
	$output = cwd;
    }

    die($errors) if $errors;
    
    return $output;
}
	    
