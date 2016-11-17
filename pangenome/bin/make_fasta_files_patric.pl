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
use TIGR::FASTArecord;
use TIGR::FASTAreader;
use TIGR::FASTAwriter;
use TIGR::Foundation;

my %opts;
GetOptions(\%opts,'fasta_dir|f=s',
           'help|h');

my $dir;
&check_params;

my $fasta = "$dir/initial_combined.fasta";
my $cmd = "cat $dir/*faa > $fasta";
system($cmd) == 0 || die("ERROR: $cmd failed");

$cmd = "/usr/local/common/cleanFasta $fasta";
system($cmd) == 0 || die("ERROR: $cmd failed");

my $tf_object = new TIGR::Foundation;
my @errors;

my $fr = new TIGR::FASTAreader ($tf_object,\@errors, $fasta );
my $fw = new TIGR::FASTAwriter($tf_object, \@errors, "$dir/combined.fasta");

while ( $fr->hasNext() ) {
    my $seq_obj = $fr->next();
    my $header= $seq_obj->{header};
    my @values = split(/\|/,$header);
    my $fasta_record = new TIGR::FASTArecord(">$values[3]",$seq_obj->{data_rec});
    $fw->write($fasta_record);
}

unlink(glob("$dir/initial_combined.fasta*"));

exit(0);

sub check_params{
    my $usage = "\nUsage: ./make_fasta_files_patric.pl --fasta_dir <directory with fasta files, default current dir> \n\n";
    
    if($opts{help}){
	print $usage;
	exit(0);
    }

    if($opts{fasta_dir}){
	$dir = $opts{fasta_dir};
    }else{
	$dir = cwd();
    }
    
}
