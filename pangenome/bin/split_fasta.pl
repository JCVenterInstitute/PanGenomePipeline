#!/usr/local/bin/perl -w

use warnings;
use strict;
$|++;

use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Cwd;
use Data::Dumper;
use File::Slurp;
use File::Path;
use FindBin;
use lib "$FindBin::Bin/../lib";
use TIGR::FASTArecord;
use TIGR::FASTAreader;
use TIGR::FASTAwriter;
use TIGR::Foundation;

my %opts;
GetOptions(\%opts,'fasta_file|f=s',
	   'max_num|n=i',
	   'output|o=s',
           'help|h');

&check_params;

my $tf_object = new TIGR::Foundation;
my @errors;

my $seq_count = 0;
my $file_count = 1;
my $output;

my $fr = new TIGR::FASTAreader ($tf_object,\@errors, $opts{fasta_file});
my $fw = new TIGR::FASTAwriter($tf_object, \@errors);
$fw->open("$output/split_fasta.$file_count");

while ( $fr->hasNext() ) {
    my $seq_obj = $fr->next();

    if($seq_count < $opts{max_num}){
	$fw->write($seq_obj);
	$seq_count++;
    }else{
	#Close and open new file
	$file_count++;
	$seq_count = 0;
	$fw->open("$output/split_fasta.$file_count");
	$fw->write($seq_obj);
	$seq_count++;
    }
}

exit(0);

sub check_params{
    my $usage = "\nUsage: ./split_fasta.pl --fasta_file <fasta file to split> --max_num <max number of sequences per file> --output <location>\n\n";
    
    if (!($opts{fasta_file}) || !($opts{max_num}) || $opts{help}){
	print STDERR "$usage";
	exit;
    }

    $output = $opts{output} // cwd;
    
    mkpath($output) unless(-d $output);
}
