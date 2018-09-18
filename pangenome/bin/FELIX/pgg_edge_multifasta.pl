#!/usr/local/bin/perl -w
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
getopts ('Dhb:m:p:a:g:');
our ($opt_D,$opt_h,$opt_b,$opt_m,$opt_p,$opt_a,$opt_g);

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver1.0";
my $basedir;
my $matchtable_file;
my $att_file;
my $pgg_file;
my $DEBUG;
my $genomes_file_name;
my @genome_array = ();
my $genome_number;
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_b) {$basedir = $opt_b;} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if (($opt_p) && (-s "$opt_p")) {$pgg_file = $opt_p;} else { print STDERR "Error with -p $opt_p\n"; &option_help; } # if no value for option p (pan-genome graph input file), quit with help menu
if (($opt_m) && (-s "$opt_m")) {$matchtable_file = $opt_m;} else { print STDERR "Error with -m $opt_m\n"; &option_help; } # if no value for option m (matchtable input file), quit with help menu
if (($opt_a) && (-s "$opt_a")) {$att_file = $opt_a;} else { print STDERR "Error with -a $opt_a\n"; &option_help; } # if no value for option a (attribute input file), quit with help menu
if (($opt_g) && (-s "$opt_g")) {$genomes_file_name = $opt_g;} else { print STDERR "Error with -g\n"; &option_help; } # if no value for option g (genome tags and contig file names input file), quit with help menu

my %feat_hash = ();         # Key1 = feat_name Key2 = struct members with their values (5p,3p,anno,gtag)
my %cluster_to_feat_hash = ();       # Key1 = genome tag Key2 = cluster_id Value = feat_name
my %genseq_hash = ();       # Key1 = genome tag Key2 = contig_name Value = contig sequence

sub get_genomes {  # obtain list of genomes - must be in the same order as the matchtable columns - and the mulitfasta contigs file for the genomes
   
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: can't open file $genomes_file_name\n");
    print "Order of genomes in $genomes_file_name with array index\n";
    while (<$infile>)  {
	chomp;
	(my $name, my $contig_file) = split(/\t/, $_);  # split the scalar $line on tab

	if (defined $genseq_hash{$name})  {
	    die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
	} else  {
	    push (@genome_array, $name); # populate the genome_array in the order of the genome file
	    print "$name\t$genome_number\n";
	    $genome_number++;
	}
	my @line = ();
	my $id;
	my $title = "";
	my $sequence = "";

	unless (open (CONTIGFILE, "<$contig_file") )  {
	    die ("cannot open file $contig_file.\n");
	}
	my ($save_input_separator) = $/;
	$/="\n>";
	while (<CONTIGFILE>) {
	    ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool)
	    @line = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @fasta
	    $id = $line[0]; # unique orf identifier is in column 0, com_name is in rest
	    $id =~ s/>//; # remove leading >
	    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
	    $genseq_hash{$name}->{$id} = $sequence;
	    $title = ""; # clear the title for the next round.
	    $sequence = ""; #clear out the sequence for the next round.
	}
	$/ = $save_input_separator; # restore the input separator
	close (CONTIGFILE);
    }
    close($infile);
    print "$genome_number genomes\n\n";
}

sub get_attributes {

    my $tag = "";
    my $end5 = "";
    my $end3 = "";
    my $asmbl_id = "";
    my $feat_name = "";
    my $anno = "";
    my $failed = 0;

    unless (open (ATTFILE, "<$att_file") )  {
	die ("ERROR: can not open file $att_file.\n");
    }
    while (<ATTFILE>) {
	my @att_line = ();
	chomp;
	@att_line = split(/\t/, $_);  # split the scalar $line on tab
	$asmbl_id = $att_line[0];
	if ($asmbl_id eq "") {
	    print STDERR "ERROR: assembly id/contig id must not be empty/null in the gene attribute file\n$_\n";
	    $failed = 1;
	}
	$feat_name = $att_line[1];
	if (defined $feat_hash{$feat_name}) {
	    print STDERR "ERROR: $feat_name appears more than once in the gene attribute file $att_file!\n";
	    $failed = 1;
	}
	$end5 = $att_line[2];
	$end3 = $att_line[3];
	$anno = $att_line[4];
	$tag = $att_line[5];
	if (!defined $genseq_hash{$tag}) {
	    print STDERR "ERROR: $tag is a genome tag in the gene attribute file $att_file but not in the genome tag file $genomes_file_name!\n";
	    $failed = 1;
	}
	$feat_hash{$feat_name}->{'5p'} = $end5;
	$feat_hash{$feat_name}->{'3p'} = $end3;
	$feat_hash{$feat_name}->{'anno'} = $anno;
	$feat_hash{$feat_name}->{'gtag'} = $tag;
	$feat_hash{$feat_name}->{'contig'} = $asmbl_id;
	print STDERR "$feat_name $feat_hash{$feat_name}->{'anno'} $feat_hash{$feat_name}->{'5p'} $feat_hash{$feat_name}->{'3p'} $feat_hash{$feat_name}->{'gtag'} $feat_hash{$feat_name}->{'contig'}\n" if ($DEBUG);
    }
    close (ATTFILE);

    if ($failed) {
	die ("ERROR: problems detected in attribute file $att_file!\n");
    }
    return;
}

sub process_matchtable {

    unless (open (TABLEFILE, "<$matchtable_file") )  {
	die ("ERROR: cannot open file $matchtable_file.\n");
    }
    while (<TABLEFILE>) {
	chomp;
	my @feat_names = split(/\t/, $_);  # split the scalar $line on tab
	my $cluster_id = shift @feat_names;
	unless (open (OUTFILE, ">$basedir/cluster_$cluster_id.fasta") )  {
	    die ("ERROR: can not open file $basedir/cluster_$cluster_id.fasta\n");
	}
	my @tmp_array = @genome_array;
	my $genome_tag;
	foreach my $feat_name (@feat_names) {
	    $genome_tag = shift @tmp_array;
	    if (($feat_name eq "----------") || ($feat_name eq "")) { #this is a placeholder and can be skipped
		next;
	    }
	    if (!defined $feat_hash{$feat_name}) { # should not happen
		die ("ERROR: gene identifier $feat_name in $matchtable_file is not in $att_file!\n");
	    }
	    if ($genome_tag ne $feat_hash{$feat_name}->{'gtag'}) {
		die ("ERROR: genome tag in $att_file (feat_hash{$feat_name}->{'gtag'}) not the same as expected column in $matchtable_file ($genome_tag)");
	    }
	    if (!defined $feat_hash{$feat_name}->{'contig'}) { # should not happen
		die ("ERROR: contig identifier was not assigned for $feat_name in $matchtable_file should have come from $att_file!\n");
	    }
	    if (!defined $genseq_hash{$genome_tag}) { # should not happen
		die ("ERROR: genome tag identifier was not assigned for $feat_name in $matchtable_file should have come from $att_file!\n");
	    }
	    if (!defined $genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}) { # should not happen
		die ("ERROR: contig sequence was not assigned for $feat_name in $matchtable_file should have come from $att_file!\n");
	    }
	    $cluster_to_feat_hash{$genome_tag}->{$cluster_id} = $feat_name;
	    my $fivep = $feat_hash{$feat_name}->{'5p'};
	    my $threep = $feat_hash{$feat_name}->{'3p'};
	    my $seq_len;
	    my $sequence;
	    print STDERR "$feat_name $feat_hash{$feat_name}->{'anno'} $feat_hash{$feat_name}->{'gtag'} $genome_tag\n" if ($DEBUG);
	    print OUTFILE ">$genome_tag";
	    print OUTFILE "_$feat_name $fivep $threep $feat_hash{$feat_name}->{'anno'}\n";
	    if ($fivep <= $threep) {
		$seq_len = ($threep - $fivep) + 1;
		$sequence = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - 1), $seq_len);
	    } else {
		$seq_len = ($fivep - $threep) + 1;
		my $tmp_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - 1), $seq_len);
		$sequence = reverse($tmp_seq);
		$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
	    }
	    if ($seq_len < 0) { #should not happen
		die ("ERROR: coordinates on contig sequence reulted in negative seq_len $seq_len for $feat_name!\n");
	    } elsif ($seq_len == 0) {
		    print OUTFILE "\n";
	    } else {
		my $pos;
		my $tmp_seq_len = $seq_len;
		for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
		    print OUTFILE substr($sequence, $pos, 60), "\n";
		    $tmp_seq_len -= 60;
		}
		print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
	    }
	}
	close (OUTFILE);
    }
    close (TABLEFILE);
    return;
}

sub process_pgg {

    unless (open (PGGFILE, "<$pgg_file") )  {
	die ("ERROR: cannot open file $pgg_file.\n");
    }
    while (<PGGFILE>) {
	chomp;
	my $cluster1;
	my $cluster2;
	my $whichend1;
	my $whichend2;
	my @edge_values = split(/\t/, $_);  # split the scalar $line on tab
	my $edge_id = shift @edge_values;
	if ($edge_id =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
	    $edge_id = "edge".$1."_".$2."to".$3."_".$4;
	    $cluster1 = $1;
	    $cluster2 = $3;
	    $whichend1 = $2;
	    $whichend2 = $4;
	} else {
	    die ("ERROR: Bad edge formatting $edge_id in file $pgg_file.\n");
	}
	if ($cluster1 > $cluster2) {
	    next; # we only need to process one orientation of an edge
	}
	    
	unless (open (OUTFILE, ">$basedir/$edge_id.fasta") )  {
	    die ("ERROR: cannot open file $basedir/$edge_id.fasta\n");
	}
	my @tmp_array = @genome_array;
	my $genome_tag;
	foreach my $edge_value (@edge_values) {
	    $genome_tag = shift @tmp_array;
	    if ($edge_value == 0) { #this is a placeholder and can be skipped
		next;
	    }
	    my $feat_name1 = $cluster_to_feat_hash{$genome_tag}->{$cluster1};
	    my $feat_name2 = $cluster_to_feat_hash{$genome_tag}->{$cluster2};
	    my $contig1 = $feat_hash{$feat_name1}->{'contig'};
	    my $contig2 = $feat_hash{$feat_name2}->{'contig'};
	    my $start1;
	    my $start2;
	    my $end1;
	    my $end2;
	    if ((!defined $feat_name1) || (!defined $feat_name2) || (!defined $contig1) || (!defined $contig2)) { # should not happen
		die ("ERROR: cluster to feat_name mapping is in conflict $genome_tag $cluster1 $cluster2!\n");
	    }
	    if ($contig1 ne $contig2) { # should not happen
		die ("ERROR: for edge $edge_value cluster features $feat_name1:$feat_name2 are not on the same contig $contig1:$contig2");
	    }
	    if (($genome_tag ne $feat_hash{$feat_name1}->{'gtag'}) || ($genome_tag ne $feat_hash{$feat_name2}->{'gtag'})) { # should not happen
		die ("ERROR: Inconsistency in genome tag between $att_file, $matchtable_file, and $pgg_file for $feat_name1 and $feat_name2");
	    }
	    if (!defined $genseq_hash{$genome_tag}) { # should not happen
		die ("ERROR: genome tag identifier was not assigned for $pgg_file should have come from $att_file!\n");
	    }
	    if ((!defined $genseq_hash{$genome_tag}->{$contig1}) || (!defined $genseq_hash{$genome_tag}->{$contig2})) { # should not happen
		die ("ERROR: contig sequence was not assigned for $contig1 $feat_name1 or $contig2 $feat_name2 in $pgg_file should have come from $att_file!\n");
	    }
	    if ($whichend1 == 5) {
		$start1 = $feat_hash{$feat_name1}->{'3p'};
		$end1 = $feat_hash{$feat_name1}->{'5p'};
	    } else {
		$start1 = $feat_hash{$feat_name1}->{'5p'};
		$end1 = $feat_hash{$feat_name1}->{'3p'};
	    }
	    if ($whichend2 == 5) {
		$start2 = $feat_hash{$feat_name2}->{'5p'};
		$end2 = $feat_hash{$feat_name2}->{'3p'};
	    } else {
		$start2 = $feat_hash{$feat_name2}->{'3p'};
		$end2 = $feat_hash{$feat_name2}->{'5p'};
	    }
	    print STDERR "$edge_value $feat_name1 $feat_name2 $genome_tag\n" if ($DEBUG);
	    print OUTFILE ">$genome_tag";
	    print OUTFILE "_$edge_value $feat_name1 $start1 $end1 $feat_name2 $start2 $end2\n";
	    my $seq_len;
	    my $sequence;
	    if ((($start1 < $end2) && (($end1 + 1) >= $start2)) || (($start1 >= $end2) && (($end1 - 1) <= $start2))) {
		# cluster/gene features overlap or abut on contig so edge sequence is empty but set length to 1 so that empty line gets output
		$seq_len = 1;
		$sequence = "";
	    } else {
		if ($start1 < $end2) {
		    $seq_len = ($start2 - $end1) - 1;
		    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $seq_len);
		} else {
		    $seq_len = ($end1 - $start2) - 1;
		    my $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $seq_len);
		    $sequence = reverse($tmp_seq);
		    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		}
	    }
	    if ($seq_len < 0) { #should not happen
		die ("ERROR: coordinates on contig sequence reulted in negative seq_len $seq_len for $edge_value $feat_name1 $start1 $end1 $feat_name2 $start2 $end2!\n");
	    } elsif ($seq_len == 0) {
		    print OUTFILE "\n";
	    } else {
		my $pos;
		my $tmp_seq_len = $seq_len;
		for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
		    print OUTFILE substr($sequence, $pos, 60), "\n";
		    $tmp_seq_len -= 60;
		}
		print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
	    }
	}
	close (OUTFILE);
    }
    close (PGGFILE);
    return;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - Pan-genome Graph, multifasta cluster and edge files
Copy (C) 2018  The J. Craig Venter Institute (JCVI).  All rights reserved

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
Example: $prog -b output_dir -g panoct_genome_tag_file -m panoct_matchtable_file -a panoct_attribute_file -p panoct_pan-genome_graph_file
Version: $version
 Option:
     -h: print this help page
     -b: base directory path for where to put output multifasta files[DEFAULT = PWD]
     -g: two colums tab delimited file (input): col1 genome_tag used in panoct in same order as panoct; col2 file name for the fasta contig file for this genome
     -m: panoct matchtable file (input)
     -a: combined panoct attribute file (input)
     -p: a panoct pan_genome_graph_file
     -D: DEBUG MODE (DEFAULT = off)
 Output: All stored within a directory specified using -b
          1) cluster_id.fasta:  a multifasta file containing the sequences in the specified cluster

 Authors: Granger Sutton, Ph.D.
  Date: September 11, 2018; last revised  September 11, 2018
  Input: The input files are primarily files generated by a panoct run and must be consistent with and across that run.
$prog requires several input files to specify the contig sequences, cluster genome tags, attribute file, matchtable, and pan-genome graph.

The input files are:

A two column tab delimited file: column 1 has the genome identifier(tag) used in the panoct run,
column 2 has the file name for a multifasta file of the contig sequences for that genome (full paths recommended).

A file of gene attributes used as input to the panoct run (combined.att).

A file containing a pan-genome graph generated by panoct.

A matchtable produced by PanOCT with the first column being a cluster identifier and subsequent columns being gene identifiers or a
placeholder indicating no gene for that genome

The output files are multifasta files, one for each cluster and edge in the pan-genome graph, containing the sequences for the clusters and edges.
_EOB_
    exit;
}

########################################  M A I N  #####################################################
print "Getting genome names and contig sequences from $genomes_file_name\n";
&get_genomes;
print "Gathering gene coordinates and attribute information from $att_file\n";
&get_attributes;
print "Reading matchtable from $matchtable_file and outputting cluster multifasta files to $basedir\n";
&process_matchtable;
print "Reading pan-genome graph from $pgg_file and outputting edge multifasta files to $basedir\n";
&process_pgg;
exit(0);
