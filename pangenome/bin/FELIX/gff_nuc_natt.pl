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
getopts ('g:');
our ($opt_g);

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver1.0";
my $genome_file_name;
if (($opt_g) && (-s "$opt_g")) {$genome_file_name = $opt_g;} else { print STDERR "Error with -g $opt_g\n"; &option_help; } # if no value for option g (genome list input file), quit with help menu

sub get_genomes {  # obtain list of genomes and ancilliary information
   
    open (my $infile, "<", "$genome_file_name") || die ("ERROR: cannot open file $genome_file_name\n");
    print "Genomes in $genome_file_name\n";
    my %genome_name_hash = (); # used to check for duplicate names
    while (my $line1 = <$infile>)  {
	my %genseq_hash = ();      # key = contig id, value = sequence of contig
	my %genseq_len = ();       # key = contig id, value = length of contig
	my %gene_id_hash = ();     # key = gene id, value = current exon count
	my @exons = ();            # array of hashs (struct) keys = contig (ctg), start coordinate (beg), attribute line (att), and stop coordinate (end)
	my @sorted_exons = ();
	my $exon_count = 0;
	chomp $line1;
	(my $name, my $gff_file_name, my $contig_file_name, my $att_file_name, my $fasta_file_name) = split(/\t/, $line1);  # split the scalar $line on tab

	if (defined $genome_name_hash{$name})  {
	    die ("ERROR: You have more than one occurance of $name in $genome_file_name!\n");
	} else  {
	    $genome_name_hash{$name} = 1;
	    print "$name\t$gff_file_name\t$contig_file_name\t$att_file_name\t$fasta_file_name\n";
	}
	my $contig_file;
	unless (open ($contig_file, "<", $contig_file_name) )  {
	    die ("ERROR: cannot open file $contig_file_name.\n");
	}
	my $save_input_separator = $/;
	$/="\n>";
	while (my $line2 = <$contig_file>) {
	    (my $title, my $sequence) = split(/\n/, $line2, 2); # split the header line and sequence (very cool)
	    my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	    my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	    $id =~ s/>\s*//; # remove leading > and spaces
	    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	    $genseq_hash{$id} = $sequence;
	    $genseq_len{$id} = length($sequence);
	    $title = ""; # clear the title for the next contig
	    $sequence = ""; #clear out the sequence for the next contig
	}
	$/ = $save_input_separator; # restore the input separator
	close ($contig_file);
	my $gff_file;
	unless (open ($gff_file, "<", $gff_file_name) )  {
	    die ("ERROR: cannot open file $gff_file_name.\n");
	}
	while (my $line2 = <$gff_file>) {
	    chomp $line2;
	    if ($line2 =~ /^##/) {# Skip header/pragma line
		next;
	    }
	    (my $seqid, my $source, my $type, my $start, my $end, my $score, my $strand, my $phase, my $attributes) = split(/\t/, $line2); # split on tab
	    (my $first, my $secod, my $third) = split(/:/, $attributes); # split attributes on :
	    if ($type eq "exon") {
		if ($first ne "ID=exon") {
		    die "ERROR: type of $type is not compatible first attribute field of $first expecting ID=exon\n";
		}
	    } else {
		next; # currently only parsing exon attributes - everything else will be on the edges
	    }
	    my @fields = split(/;/, $third);  # split third field of attributes on ;
	    foreach my $field (@fields) {
		(my $tag, my $value) = split(/=/, $field);  # split the pair on =
		if ($tag eq "gene_id") {
		    my $exon_number;
		    my $contig_len = $genseq_len{$seqid};
		    if (($start < 1) || ($end < 1) || ($start > $contig_len) || ($end > $contig_len)) {
			die "ERROR: feature coordinate falls outside of contig boudaries (1:$contig_len) $seqid $start $end\n";
		    }
		    if (defined $gene_id_hash{$value}) {
			$exon_number = $gene_id_hash{$value};
			$exon_number++;
		    } else {
			$exon_number = 1;
		    }
		    $gene_id_hash{$value} = $exon_number;
		    my $feat_name = $value . "_exon_" . $exon_number;
		    if ($value !~ /^$name/) { #append genome name if needed
			$feat_name = $name . "_" . $feat_name;
		    }
		    my $anno = "";
		    my $attribute
		    $attribute = "$seqid\t$feat_name\t";
		    if ($strand eq "+") {
			$attribute .= "$start\t$end\t";
		    } else {
			$attribute .= "$end\t$start\t";
		    }
		    $attribute .= "$anno\t$name\n";
		    $exons[$exon_count++] = { 'ctg' => $seqid,       # contig identifier
					      'beg' => $start,       # start coordinate on contig
					      'att' => $attribute,   # attribute line
					      'end' => $end          # end  coordinate on contig
		    }
		    last;
		}
	    }
	}
	close ($gff_file);
	my $sort_by_contig_start = sub { # sort by contig, then start coordinate, then stop coordinate
	    my $contig_test = $a->{'ctg'} cmp $b->{'ctg'};
	    
	    if ($contig_test) {
		return ($contig_test);
	    } elsif ($a->{'beg'} <=> $b->{'beg'}) {
		return ($a->{'beg'} <=> $b->{'beg'});
	    } else {
		return (($b->{'end'} - $b->{'beg'}) <=> ($a->{'end'} - $a->{'beg'}));
	    }
	};
	
	@sorted_exons = sort $sort_by_contig_start (@exons);
	
	my $att_file;
	unless (open ($att_file, ">", $att_file_name) )  {
	    die ("ERROR: cannot open file $att_file_name.\n");
	}
	my $fasta_file;
	unless (open ($fasta_file, ">", $fasta_file_name) )  {
	    die ("ERROR: cannot open file $fasta_file_name.\n");
	}
	foreach my $i (0 .. $#sorted_exons) {
	    my $ctg1 = $sorted_exons[$i]->{'ctg'};
	    my $beg1 = $sorted_exons[$i]->{'beg'};
	    my $end1 = $sorted_exons[$i]->{'end'};
	    if (!defined $sorted_exons[$i]->{'dup'}) {
		    print $att_file $sorted_exons[$i]->{'att'};
		    my $sequence;
		    my $seq_len;
		    my $anno = "";
		    if ($strand eq "+") {
			$seq_len = ($end - $start) + 1;
			$sequence = substr($genseq_hash{$seqid}, ($start - 1), $seq_len);
		    } else {
			$seq_len = ($end - $start) + 1;
			my $tmp_seq = substr($genseq_hash{$seqid}, ($start - 1), $seq_len);
			$sequence = reverse($tmp_seq);
			$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		    my @fields = split(/\t/, $sorted_exons[$i]->{'att'});  # $fields[1] contains feat_name
		    print $fasta_file ">$fields[1]\n";
		    my $pos;
		    my $tmp_seq_len = $seq_len;
		    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
			print $fasta_file substr($sequence, $pos, 60), "\n";
			$tmp_seq_len -= 60;
		    }
		    print $fasta_file substr($sequence, $pos, $tmp_seq_len), "\n";
	    }
	    foreach my $j (($i + 1) .. $#sorted_exons) {
		my $ctg2 = $sorted_exons[$j]->{'ctg'};
		my $beg2 = $sorted_exons[$j]->{'beg'};
		my $end2 = $sorted_exons[$j]->{'end'};
		if (($ctg1 ne $ctg2) || ($beg2 > $end1)){
		    last;
		}
		if ($end2 <= $end1) {
		    print "WARNING: Ignoring:\n$sorted_exons[$j]->{'att'} which is contained by $sorted_exons[$i]->{'att'}";
		    $sorted_exons[$j]->{'dup'} = 1;
		}
	    }
	}
	close ($att_file);
	close ($fasta_file);
    }
    close($infile);
    return;
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog  - input is list of genomes, gff files, genome files, output is attribute and gene fasta files
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
Example: $prog -g genome list file
Version: $version
 Option:
     -h: print this help page
     -g: genome list file name: five colums tab delimited file (input):
                                col1 genome name to use in atrribute output file
                                col2 file name for the fasta GFF input file for this genome
                                col3 file name for the fasta contig inputfile for this genome
                                col4 file name for the attribute output file for this genome
                                col5 file name for the gene fasta output file for this genome
 Output: An attribute file and a gene fasta file for input to panOCT

 Authors: Granger Sutton, Ph.D.
  Date: September 11, 2018; last revised  September 11, 2018
  Input: The genome list file as specified above which specifies a GFF file and genome contig fasta file for each genome

_EOB_
    exit;
}

########################################  M A I N  #####################################################
print "Getting genome names and contig sequences from $genome_file_name\n";
&get_genomes;
exit(0);
