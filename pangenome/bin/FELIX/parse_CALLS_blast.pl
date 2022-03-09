#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

##### Attribute Compare Script

use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;
use List::Util qw[min max];

my $cwd = getcwd;
my $commandline = join (" ", @ARGV);
print STDERR "$commandline\n";
my %plasmids = (); # key = plasmid name, value = plasmid sequence
my %foreign_plasmids = (); # key = plasmid name, value = 1 to indicate this plasmid was deemed a foreign plasmid
my %insertion_events = (); # key = insertion event name, value = CALLs line for the insertion event

# CONSTANTS #
use constant QSEQID => 0;
use constant SSEQID => 1;
use constant PIDENT => 2;
use constant QSTART => 3;
use constant QEND => 4;
use constant QLEN => 5;
use constant SSTART => 6;
use constant SEND => 7;
use constant SLEN => 8;
use constant EVALUE => 9;
use constant BITSCORE => 10;
use constant STITLE => 11;
use constant CSAMPLE => 0;
use constant CTYPE => 1;
use constant CQID => 2;
use constant CFIVEP => 3;
use constant CTHREEP => 4;
use constant CDELETED => 5;
use constant CDLEN => 6;
use constant CINSERTED => 7;
use constant CILEN => 8;
use constant CSID => 9;
use constant CPID => 10;
use constant CSSTART => 11;
use constant CSEND => 12;
use constant CQSTART => 13;
use constant CQEND => 14;
# END CONSTANTS #

my $help;
my $debug;
my $engineering_found;
my $files;
my $output_prefix = "OUTPUT_CALLS";

GetOptions('files=s' => \ $files,
	   'prefix=s' => \ $output_prefix,
	   'help' => \ $help,
	   'debug' => \ $debug);

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('files=s' => files,
	   'prefix=s' => output_prefix,
	   'help' => help,
	   'debug' => debug);
_EOB_
    exit(0);
}

my $out_clear_file = $out_prefix . "_clear.txt";
my $out_stop_file = $out_prefix . "_stop_codons.txt";
my $out_maybe_file = $out_prefix . "_maybe.txt";

# open output files
open (my $out_clear, ">", $out_clear_file) || die ("ERROR: cannot open output file $out_clear_file\n");
open (my $out_stop, ">", $out_stop_file) || die ("ERROR: cannot open output file $out_stop_file\n");
open (my $out_maybe, ">", $out_maybe_file) || die ("ERROR: cannot open output file $out_maybe_file\n");

# read file which specifies the Sample ID, CALLS file, CALLS INSERTIONS Blastn tabular results, Plasmids fasta file, Plasmids Blastn tabular results
open (my $infile, "<", $files) || die ("ERROR: cannot open input file $files\n");
while (my $line = <$infile>)  {
    # clear data structures for the next genome
    %plasmids = (); # key = plasmid name, value = plasmid sequence
    %foreign_plasmids = (); # key = plasmid name, value = 1 to indicate this plasmid was deemed a foreign plasmid
    %insertion_events = (); # key = insertion event name, value = CALLs line for the insertion event
    
    chomp $line;
    (my $sample_name, my $species_name, my $calls_file_name, my $calls_btab_file_name, my $plasmids_fasta_file_name, my $plasmids_btab_file_name, my $stop_codons_file_name) = split(/\t/, $line);  # split the scalar $line on tab
    $engineering_found = 0;
    
    # read in plasmids from plasmids fasta file
    my $plasmids_fasta_file;
    unless (open ($plasmids_fasta_file, "<", $plasmids_fasta_file_name) )  {
	die ("cannot open plasmids fasta file: $plasmids_fasta_file_name!\n");
    }
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    while ($line = <$plasmids_fasta_file>) {
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$id =~ s/>\s*//; # remove leading > and spaces
	$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	if ($id =~ /recover/) {
	    # ignore recovered read contigs
	} else {
	    $plasmids{$id} = $sequence;
	}
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($plasmids_fasta_file);


    # read in calls file
    my $calls_file;
    open($calls_file, "<", $calls_file_name) || die ("Couldn't open calls file $calls_file_name\n");
    my $first = 1;
    while (my $line = <$calls_file>) {
	if ($first) {
	    $first = 0;
	    next; # first line is a header
	}
	chomp($line);
	my @calls_line = split(/\t/,$line);
	my $category = $calls_line[CTYPE];
	(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	if ($contig_id =~ /recover/) {
	    next; # ignore recovered read contigs
	} 
	if (($category =~ /INSERTION/) || ($category =~ /PARTIAL INSERTION/)) {
	    # Handle an insertion event
	    $insertion_events{$calls_line[CQID]} = $line;
	} elsif ($category =~ /DELETION/) {
	    # Handle a deletion event
	    $engineering_found = 1;
	    print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tdeleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\tNA\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\t$calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	} elsif ($category =~ /MUTATION/) {
	    print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tmutation : inserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\t\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	} elsif ($category =~ /TANDEM_DUPLICATION/) {
	    print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\ttandem duplication : inserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\t\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	} else {
	    next; # ignore other types
	}
    }
    close($calls_file);

    my $plasmids_btab_file;
    open($plasmids_btab_file, "<", $plasmids_btab_file_name) || die ("Couldn't open plasmids btab file $plasmids_btab_file_name for reading\n");
    my $prev_qid = "";
    my $max_bitscore = 0;
    my $max_title = "";
    my $max_sid = "";
    my $max_length = 0;
    my $found = 0;
    my $results = "";
    while (my $line = <$plasmids_btab_file>) {
	if ($line =~ /^#/) {
	    next; # skip comment lines
	}
	chomp($line);
	my @split_line = split(/\t/,$line);
	my $qid = $split_line[0];         # query id
	my $sid = $split_line[1];         # subject id
	my $pid = $split_line[2];         # percent identity
	my $qstart = $split_line[3];      # query start
	my $qend = $split_line[4];          # query end
	my $qlen = $split_line[5];          # query length
	my $sstart = $split_line[6];      # subject start
	my $send = $split_line[7];          # subject end
	my $slen = $split_line[8];          # subject length
	my $evalue = $split_line[9];      # evalue
	my $bitscore = $split_line[10]; # bitscore
	my $stitle = $split_line[11];     # subject title
	if (scalar(@split_line) != 12) {
	    print STDERR "WARNING: malformed btab line:\n$line\n";
	    next;
	}
	if ($qid =~ /recover/) {
	    next; # ignore recovered read contigs
	}
	if (!defined $plasmids{$qid}) {
	    print STDERR "WARNING: contig $qid in the plasmids btab file but not in the plasmids fasta file\n";
	    next; # ignore contigs not in the plasmids fasta file
	}
	if (($prev_qid ne $qid) && ($prev_qid ne "")) {
	    if ($found) {
		$engineering_found = 1;
		$foreign_plasmids{$prev_qid} = 1;
		print $out_clear "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
	    } else {
		print $out_maybe "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
	    }
	    $max_bitscore = 0;
	    $max_title = "";
	    $max_sid = "";
	    $max_length = 0;
	    $found = 0;
	    $results = "";
	}
	$prev_qid = $qid;
	if ($results ne "") {
	    $results .= ":";
	}
	$results .= join(',', @split_line);
	if (($stitle =~ /vector/i) || ($stitle =~ /construct/i) || ($stitle =~ /synthetic/i)) {
	    $found = 1;
	}
	if ($bitscore > $max_bitscore) {
	    $max_bitscore = $bitscore;
	    $max_title = $stitle;
	    $max_sid = $sid;
	    $max_length = $qlen;
	}
    }
    if ($found) {
	$engineering_found = 1;
	$foreign_plasmids{$prev_qid} = 1;
	print $out_clear "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
    } elsif ($prev_qid ne "") {
	print $out_maybe "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
    }
    close($plasmids_btab_file);

    my $calls_btab_file;
    open($calls_btab_file, "<", $calls_btab_file_name) || die ("Couldn't open CALLs btab file $calls_btab_file_name for reading\n");
    $prev_qid = "";
    $max_bitscore = 0;
    $max_title = "";
    $max_sid = "";
    $max_length = 0;
    $found = 0;
    $results = "";
    while (my $line = <$calls_btab_file>) {
	if ($line =~ /^#/) {
	    next; # skip comment lines
	}
	chomp($line);
	my @split_line = split(/\t/,$line);
	my $qid = $split_line[0];         # query id
	my $sid = $split_line[1];         # subject id
	my $pid = $split_line[2];         # percent identity
	my $qstart = $split_line[3];      # query start
	my $qend = $split_line[4];          # query end
	my $qlen = $split_line[5];          # query length
	my $sstart = $split_line[6];      # subject start
	my $send = $split_line[7];          # subject end
	my $slen = $split_line[8];          # subject length
	my $evalue = $split_line[9];      # evalue
	my $bitscore = $split_line[10]; # bitscore
	my $stitle = $split_line[11];     # subject title
	if (scalar(@split_line) != 12) {
	    print STDERR "WARNING: malformed btab line:\n$line\n";
	    next;
	}
	if ($qid =~ /recover/) {
	    next; # ignore recovered read contigs
	}
	if (!defined $insertion_events{$qid}) {
	    print STDERR "WARNING: contig $qid in the CALLs btab file but not in the CALLs file\n";
	    next; # ignore contigs not in the CALLs file
	}
	if (($prev_qid ne $qid) && ($prev_qid ne "")) {
	    if ($found) {
		my @calls_line = split(/\t/,$insertion_events{$prev_qid});
		(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
		if (!defined $foreign_plasmids{$contig_id}) {
		    # only do this if this was not already called a foreign plasmid
		    $engineering_found = 1;
		    print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
		}
	    } elsif (!defined $foreign_plasmids{$contig_id}) {
		print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	    }
	    $max_bitscore = 0;
	    $max_title = "";
	    $max_sid = "";
	    $max_length = 0;
	    $found = 0;
	    $results = "";
	}
	$prev_qid = $qid;
	if ($results ne "") {
	    $results .= ":";
	}
	$results .= join(',', @split_line);
	if (($stitle =~ /vector/i) || ($stitle =~ /construct/i) || ($stitle =~ /synthetic/i)) {
	    $found = 1;
	}
	if ($bitscore > $max_bitscore) {
	    $max_bitscore = $bitscore;
	    $max_title = $stitle;
	    $max_sid = $sid;
	    $max_length = $qlen;
	}
    }
    if ($found) {
	my @calls_line = split(/\t/,$insertion_events{$prev_qid});
	(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	if (!defined $foreign_plasmids{$contig_id}) {
	    # only do this if this was not already called a foreign plasmid
	    $engineering_found = 1;
	    print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	}
    } elsif ((!defined $foreign_plasmids{$contig_id} && ($prev_qid ne "")) {
	print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
    }
    if (!$engineering_found) {
	print $out_clear "$sample_name\tyes\tPGG\tno\tNA\t\tcomparison to PGG\t$species_name\tNA\tNA\tNA\tNA\tNA\tNA\t\n";
    }
    close($calls_btab_file);

    my $stop_codons_file;
    open($stop_codons_file, "<", $stop_codons_file_name) || die ("Couldn't open stop codons file $stop_codons_file_name for reading\n");
    while (my $line = <$stop_codons_file>) {
	chomp($line);
	my @split_line = split(/\t/,$line);
	my $target_id = $split_line[0];         # target genome id
	my $contig_id = $split_line[1];         # contig id
	my $type_stop_codon = $split_line[2];         # currently always possible_stop_codon
	my $qstart = $split_line[3];      # query start
	my $qend = $split_line[4];          # query end
	my $qlen = $split_line[5];          # query length
	my $feat_name = $split_line[6];          # feature name
	my $frame = $split_line[7];      # which reading frame
	my $stop_codon_coord = $split_line[8];      # stop codon coordinate on contig
	my $feat_seq = $split_line[9]; # feature sequence
	if ($contig_id =~ /recover/) {
	    next; # ignore recovered read contigs
	}
	print $out_stop "$sample_name\tyes\tPGG feature medoids\tyes\t$type_stop_codon in $feat_name : sequence : $feat_seq\t\tcomparison to PGG feature medoids\t$species_name\tNA\tContig $contig_id: feature $feat_name coordinates $qstart,$qend: stop codon contig coordinate $stop_codon_coord\tfeature length $qlen : stop codon length 1-3bp\t\t\t\t\n";
    }
    close($stop_codons_file);
}

close($out_clear);
close($out_stop);
close($out_maybe);

exit(0);
