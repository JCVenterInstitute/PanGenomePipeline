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
use File::Basename;

my %opts;
my ($DIR,$GENOMES);

GetOptions(\%opts,'cds_dir|c=s',
	   'db_list|d=s',
           'help|h');

&check_params;

my @cds_files = glob("$DIR/*cds.tab");
my @errors;

#Check names first
foreach my $file(@cds_files){
    my($name,$path,$suffix) = fileparse($file);
    my @values = split(/\./,$name);
    
    unless(exists $GENOMES->{$values[0]}){
	push(@errors, $values[0]);
    }
}

if(scalar @errors > 0){
    my $msg = "The following genome names parsed from the cds.tab files do not match a name in list.txt\n";
    $msg .= "Please model the name in list.txt to what cds.tab file is named\n";

    print "$msg\n";

    foreach my $e (@errors){
	print "$e\n";
    }
    
    exit(1);
}

open(my $ofh, ">","$DIR/combined.att_file");

my $asmbl_id_hsh;

foreach my $file(@cds_files){

    open(my $fh, "<", $file);

    my($name,$path,$suffix) = fileparse($file);
    my @file_names = split(/\./,$name);
    my $current_molecule = "";
    my $current_mol_id = 1;

    while(<$fh>){

	my @values = split(/\t/,$_);

	unless(/GENOME_NAME/ ~~ @values){
	    if($values[0] =~ /^\w/){
		
		if($current_molecule){
		    if($current_molecule eq $values[1]){
			print $ofh $asmbl_id_hsh->{$file}->{$values[1]} . "\t";
		    }else{
			$current_mol_id++;
			$asmbl_id_hsh->{$file}->{$values[1]} = $current_mol_id;
			print $ofh $asmbl_id_hsh->{$file}->{$values[1]} . "\t";
		    }
		}else{
		    $asmbl_id_hsh->{$file}->{$values[1]} = $current_mol_id;
		    print $ofh $asmbl_id_hsh->{$file}->{$values[1]} . "\t";
		}

		$current_molecule = $values[1];

		print $ofh $values[5] . "\t"; #Locus

		#Coords
		if($values[8] eq '-'){
		    print $ofh $values[7] . "\t";
		    print $ofh $values[6] . "\t";
		}else{
		    print $ofh $values[6] . "\t";
		    print $ofh $values[7] . "\t";
		}

		$values[11] =~ s/\s+$//;
		print $ofh $values[11] . "\t"; #Protein Name
		print $ofh $file_names[0] . "\n"; #Genome
	    }
	}	
    }
}
close $ofh;

open(my $afh,">","$DIR/asmbl_id_lookup.txt");
foreach my $file(keys %$asmbl_id_hsh){
    my ($name,$path,$post) = fileparse($file);
    my @values = split(/\./,$name);
    
    print $afh "$values[0]\n";
    foreach my $molecule(sort keys %{$asmbl_id_hsh->{$file}}){
	print $afh "\t$molecule\t$asmbl_id_hsh->{$file}->{$molecule}\n";
    }
}
close $afh;

###

sub check_params{
    my $usage = "./make_combined_att_file_patric.pl --db_list <db list> --cds_dir <location of cds files>\n";
    
    if($opts{cds_dir}){
	$DIR = $opts{cds_dir};
    }else{
	$DIR = cwd();
    }

    if($opts{db_list}){
	my @genomes = read_file($opts{db_list});
	foreach my $g(@genomes){
	    $g =~ s/\s+$//;
	    $g =~ s/\/$//;
	    $GENOMES->{$g} = 1;
	}
    }else{
	print $usage;
	exit (0);
    }

    if($opts{help}){
	print $usage;
	exit (0);
    }
}
