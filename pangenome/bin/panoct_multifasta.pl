#!/usr/bin/env perl
#Copy (C) 2013  The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Granger Sutton, Ph.D.

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

#Revision notes
my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use warnings;
use Getopt::Std;
getopts ('Dhb:M:P:C:');
our ($opt_D,$opt_h,$opt_b,$opt_M,$opt_P,$opt_C);

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver3_17GGS";
my $basedir;
my $pep_file;
my $matchtable_file;
my $cluster_file;
my $DEBUG;
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_b) {$basedir = $opt_b;} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if (($opt_P) && (-s "$opt_P")) {$pep_file = $opt_P;} else { print STDERR "Error with -P $opt_P\n"; &option_help; } # if no value for option P (multifasta input file), quit with help menu
if (($opt_M) && (-s "$opt_M")) {$matchtable_file = $opt_M;} else { print STDERR "Error with -M $opt_M\n"; &option_help; } # if no value for option M (matchtable input file), quit with help menu
if (($opt_C) && (-s "$opt_C")) {$cluster_file = $opt_C;} else { print STDERR "Error with -C $opt_C\n"; &option_help; } # if no value for option  (cluster number input file), quit with help menu

my %feat_hash = ();         # Key = feat_name Key2 = struct members with their values
my %cluster_ids = ();       # Key = cluster id Value = 1 just to define it

sub get_cluster_ids {  # obtain list of genomes to compareread cluster ids into a hash
   
    open (my $infile, "<", "$cluster_file") || die ("ERROR: can't open file $cluster_file\n");
    while (<$infile>)  {
	chomp;
	my $cluster_id = $_;
	$cluster_id =~ s/\r$//; #strip ^M for Windows files
	$cluster_id =~ s/\s+$//; #strip trailing white space
	$cluster_id =~ s/\r$//; #do it again in case
	if (defined $cluster_ids{$cluster_id})  {
               die ("ERROR:  You have more than one occurance of $cluster_id in $cluster_file!\n");
	}
	else  {
	    $cluster_ids{$cluster_id} = 1; # used to be genome_hash
	}
    }  
    close($infile);
    return;
}

sub get_protein_info { #read in a multifasta sequence file storing the fasta header and the sequence

  my @line = ();
  my $id;
  my $title = "";
  my $sequence = "";
  my $length = "";

  unless (open (PEPFILE, "<$pep_file") )  {
    die ("can't open file $pep_file.\n");
  }
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<PEPFILE>) {
    ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool)
    @line = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @fasta
    $id = $line[0]; # unique orf identifier is in column 0, com_name is in rest
    $id =~ s/>//;
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
    $length = length($sequence);
    $feat_hash{$id}->{'header'} = join(' ', @line[1..$#line]); # put the identifier into the hash as the "key" and the header as value "header" (joining all columns after first space/tab)
    $feat_hash{$id}->{'length'}= $length;
    $feat_hash{$id}->{'sequence'}= $sequence;
    #print STDERR "$id ($feat_hash{$id}->{'header'}) = $feat_hash{$id}->{'length'}\n";
    $title = ""; # clear the title for the next round.
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (PEPFILE);
  return;
}

sub process_matchtable {

    unless (open (TABLEFILE, "<$matchtable_file") )  {
	die ("ERROR: can not open file $matchtable_file.\n");
    }
    while (<TABLEFILE>) {
	my @feat_names = ();
	chomp;
	@feat_names = split(/\t/, $_);  # split the scalar $line on tab
	my $cluster_id = shift @feat_names;
	if (!defined $cluster_ids{$cluster_id}) { #skip clusters not in the cluster_file
	    next;
	}
	unless (open (OUTFILE, ">$basedir/$cluster_id") )  {
	    die ("ERROR: can not open file $basedir/$cluster_id.\n");
	}
	foreach my $feat_name (@feat_names) {
	    if ($feat_name eq "----------") { #this is a placeholder and can be skipped
		next;
	    }
	    if (!defined $feat_hash{$feat_name}) { # should not happen
		die ("ERROR: gene identifier $feat_name in $matchtable_file is not in $pep_file!\n");
	    }
	    print STDERR "$feat_name $feat_hash{$feat_name}->{'header'} feat_hash{$feat_name}->{'length'}\n" if ($DEBUG);
	    print OUTFILE ">$feat_name $feat_hash{$feat_name}->{'header'}\n";
	    my $seq_len = $feat_hash{$feat_name}->{'length'};
	    my $sequence = $feat_hash{$feat_name}->{'sequence'};
	    for ( my $pos = 0 ; $pos < $seq_len ; $pos += 60 ) {
		print OUTFILE substr($sequence, $pos, 60), "\n";
	    }
	}
	close (OUTFILE);
    }
    close (TABLEFILE);
    return;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - Pan-genome Ortholog Clustering Tool, multifasta cluster files from cluster ids
Copy (C) 2013  The J. Craig Venter Institute (JCVI).  All rights reserved

License:   This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

Citation:  Derrick E. Fouts, Lauren Brinkac, Erin Beck, Jason Inman, and Granger Sutton (2011) "PanOCT: Automated Clustering of Orthologs using 
           Conserved Gene Neighborhood for Pan-Genomic Analysis of Bacterial Strains and Closely Related Species" Nucleic Acids Res. 2012 Dec;40(22).

  Usage: $prog <options>
Example: panoct_multifasta.pl -b output_dir -P multifasta_input_file -M matchtable_file -C cluster_ids_files
Version: $version
 Option:
     -h: print this help page
     -b: base directory path for where to put output multifasta files[DEFAULT = PWD]
     -P: multifasta input file of all genes input to PanOCT (can be peptide or nucleotide)
     -M: matchtable input file
     -C: cluster identifiers input file (this should usually be numbers)
     -D: DEBUG MODE (DEFAULT = off)
 Output: All stored within a directory specified using -b
          1) cluster_id.fasta:  a multifasta file containing the sequences in the specified cluster

 Authors: Granger Sutton, Ph.D.
 Date: July 29, 2013; last revised July 29, 2013
panoct_multifasta.pl requires several input files to specify the gene sequences, cluster identifiers, and clusters.

The input files are:

A multifasta file containing all of the gene sequences input to PanOCT (peptide or nucleotide) using the same gene identifiers as
were used for PanOCT for instance the peptide file specified to PanOCT using -P. The file is specified by the -P option.

A file of cluster identifiers. The file is specified by the -C option. The cluster identifiers must be placed one per line in the
file with nothing else on the line. The cluster identifiers must be the same as in the first column of the matchtable - default for
PanOCT is to number the clusters from 1-N.

A matchtable produced by PanOCT with the first column being a cluster identifier and subsequent columns being gene identifiers or a
placeholder indicating no gene for that genome

The output files are multifasta files, one per cluster identifier specified, containing the gene sequences for the cluster.
_EOB_
    exit;
}

########################################  M A I N  #####################################################
print STDERR "Reading cluster identifiers from $cluster_file\n";
&get_cluster_ids;
print STDERR "Gathering sequence information from $pep_file\n";
&get_protein_info;
print STDERR "Reading matchtable from $matchtable_file and outputting cluster multifasta files to $basedir\n";
&process_matchtable;
exit(0);
