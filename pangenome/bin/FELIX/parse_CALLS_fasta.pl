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
my %event_counts = (); # key1 = species, key2 = sample name / total, value = array of counts
my @event_counts_total = [0,0,0,0,0,0,0,0]; # array of total counts

# CONSTANTS #
use constant TINSC => 0;
use constant TFPC => 1;
use constant TDEL => 2;
use constant TINSM => 3;
use constant TFPM => 4;
use constant TTD => 5;
use constant TMUT => 6;
use constant TSC => 7;
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

# subroutine to print sequence out in fasta format
sub print_fasta {

    my ($file_handle, $seq_name, $seq) = @_;
    #print STDERR "$seq_name\n";
    print $file_handle ">$seq_name\n";
    my $tmp_pos;
    my $tmp_seq_len = length($seq);
    for ( $tmp_pos = 0 ; $tmp_seq_len > 60 ; $tmp_pos += 60 ) {
	print $file_handle substr($seq, $tmp_pos, 60), "\n";
	$tmp_seq_len -= 60;
    }
    print $file_handle substr($seq, $tmp_pos, $tmp_seq_len), "\n";
    return;
}

my $out_clear_file = $output_prefix . "_clear.txt";
my $out_stop_file = $output_prefix . "_stop_codons.txt";
my $out_maybe_file = $output_prefix . "_maybe.txt";
my $out_fasta_file = $output_prefix . "_seqs.fasta";
my $out_short_file = $output_prefix . "_all_summary.txt";
my $out_summary_file = $output_prefix . "_counts_summary.txt";

# open output files
open (my $out_clear, ">", $out_clear_file) || die ("ERROR: cannot open output file $out_clear_file\n");
open (my $out_stop, ">", $out_stop_file) || die ("ERROR: cannot open output file $out_stop_file\n");
open (my $out_maybe, ">", $out_maybe_file) || die ("ERROR: cannot open output file $out_maybe_file\n");
open (my $out_fasta, ">", $out_fasta_file) || die ("ERROR: cannot open output file $out_fasta_file\n");
open (my $out_short, ">", $out_short_file) || die ("ERROR: cannot open output file $out_short_file\n");
open (my $out_summary, ">", $out_summary_file) || die ("ERROR: cannot open output file $out_summary_file\n");

# read file which specifies the Sample ID, CALLS file, CALLS INSERTIONS Blastn tabular results, Plasmids fasta file, Plasmids Blastn tabular results
open (my $infile, "<", $files) || die ("ERROR: cannot open input file $files\n");
while (my $line = <$infile>)  {
    # clear data structures for the next genome
    %plasmids = (); # key = plasmid name, value = plasmid sequence
    %foreign_plasmids = (); # key = plasmid name, value = 1 to indicate this plasmid was deemed a foreign plasmid
    %insertion_events = (); # key = insertion event name, value = CALLs line for the insertion event
    
    chomp $line;
    (my $sample_name, my $species_name, my $calls_file_name, my $calls_btab_file_name, my $plasmids_fasta_file_name, my $plasmids_btab_file_name, my $stop_codons_file_name) = split(/\t/, $line);  # split the scalar $line on tab
    my $species_name_ns = $species_name;
    $species_name_ns =~ s/\s+/_/g;
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
	$calls_line[CTYPE] =~ s/\s+/_/g;
	(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	if ($contig_id =~ /recover/) {
	    next; # ignore recovered read contigs
	} 
	if (($category =~ /INSERTION/) || ($category =~ /PARTIAL INSERTION/)) {
	    # Handle an insertion event
	    $insertion_events{$calls_line[CQID]} = $line;
	    #print STDERR "SAW: $line\n";
	} elsif ($category =~ /DELETION/) {
	    # Handle a deletion event
	    $engineering_found = 1;
	    print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tdeleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\tNA\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\t$calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tclear\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TDEL]++;
	    $event_counts{$species_name}->{"Total"}->[TDEL]++;
	    $event_counts_total[TDEL]++;
	    if ($calls_line[CDLEN] >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] clear deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
	    }
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] clear 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] clear 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
	} elsif ($category =~ /MUTATION/) {
	    print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tmutation : inserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\t\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tmaybe\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TMUT]++;
	    $event_counts{$species_name}->{"Total"}->[TMUT]++;
	    $event_counts_total[TMUT]++;
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] maybe deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] maybe inserted sequence $calls_line[CQID]$contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
	} elsif ($category =~ /TANDEM_DUPLICATION/) {
	    print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\ttandem duplication : inserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\t\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tmaybe\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TTD]++;
	    $event_counts{$species_name}->{"Total"}->[TTD]++;
	    $event_counts_total[TTD]++;
	    if ($calls_line[CILEN] >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] maybe tandem duplicated sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
	    }
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] maybe 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] maybe 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
	} else {
	    print STDERR "UNEXPECTED CALL TYPE: $category $line\n";
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
		print $out_short "$sample_name\t$species_name\tclear\tFOREIGN PLASMID\t$prev_qid\t1\t$max_length\n";
		if (!defined $event_counts{$species_name}) {
		    $event_counts{$species_name} = {};
		    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
		}
		if (!defined $event_counts{$species_name}->{$sample_name}) {
		    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
		}
		$event_counts{$species_name}->{$sample_name}->[TFPC]++;
		$event_counts{$species_name}->{"Total"}->[TFPC]++;
		$event_counts_total[TFPC]++;
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$prev_qid clear foreign plasmid $max_title", $plasmids{$prev_qid});
	    } else {
		$foreign_plasmids{$prev_qid} = 1;
		print $out_maybe "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
		print $out_short "$sample_name\t$species_name\tmaybe\tFOREIGN PLASMID\t$prev_qid\t1\t$max_length\n";
		if (!defined $event_counts{$species_name}) {
		    $event_counts{$species_name} = {};
		    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
		}
		if (!defined $event_counts{$species_name}->{$sample_name}) {
		    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
		}
		$event_counts{$species_name}->{$sample_name}->[TFPM]++;
		$event_counts{$species_name}->{"Total"}->[TFPM]++;
		$event_counts_total[TFPM]++;
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$prev_qid maybe foreign plasmid $max_title", $plasmids{$prev_qid});
	    }
	    delete $plasmids{$prev_qid};
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
	if (($stitle =~ /vector/i) || ($stitle =~ /construct/i) || ($stitle =~ /synthetic/i) || ($stitle =~ /marker plasmid/i)) {
	    $found = 1;
	}
	if ($bitscore > $max_bitscore) {
	    $max_bitscore = $bitscore;
	    $max_title = $stitle;
	    $max_sid = $sid;
	    $max_length = $qlen;
	}
    }
    if ($prev_qid ne "") {
	if ($found) {
	    $engineering_found = 1;
	    $foreign_plasmids{$prev_qid} = 1;
	    print $out_clear "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tclear\tFOREIGN PLASMID\t$prev_qid\t1\t$max_length\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TFPC]++;
	    $event_counts{$species_name}->{"Total"}->[TFPC]++;
	    $event_counts_total[TFPC]++;
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$prev_qid clear foreign plasmid $max_title", $plasmids{$prev_qid});
	} else {
	    $foreign_plasmids{$prev_qid} = 1;
	    print $out_maybe "$sample_name\tyes\t$max_sid $max_title\tyes\tforeign plasmid : $plasmids{$prev_qid}\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $prev_qid\t$max_length\tplasmid\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tmaybe\tFOREIGN PLASMID\t$prev_qid\t1\t$max_length\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TFPM]++;
	    $event_counts{$species_name}->{"Total"}->[TFPM]++;
	    $event_counts_total[TFPM]++;
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$prev_qid maybe foreign plasmid $max_title", $plasmids{$prev_qid});
	}
	delete $plasmids{$prev_qid};
    }
    close($plasmids_btab_file);
    foreach my $key (keys %plasmids) {
	$max_length = length($plasmids{$key});
	print $out_maybe "$sample_name\tyes\tunknown plasmid\tyes\tforeign plasmid : $plasmids{$key}\t\tcomparison to PGG\t$species_name\tunknown plasmid\tContig $key\t$max_length\tplasmid\t\t\t\n";
	print $out_short "$sample_name\t$species_name\tmaybe\tFOREIGN PLASMID\t$key\t1\t$max_length\n";
	if (!defined $event_counts{$species_name}) {
	    $event_counts{$species_name} = {};
	    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	}
	if (!defined $event_counts{$species_name}->{$sample_name}) {
	    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	}
	$event_counts{$species_name}->{$sample_name}->[TFPM]++;
	$event_counts{$species_name}->{"Total"}->[TFPM]++;
	$event_counts_total[TFPM]++;
	print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$key maybe foreign plasmid no match", $plasmids{$key});
    }

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
	    print STDERR "SKIPPING recovered reads contig: $line\n";
	    next; # ignore recovered read contigs
	}
	if (!defined $insertion_events{$qid}) {
	    print STDERR "WARNING: contig $qid in the CALLs btab file but not in the CALLs file\n";
	    next; # ignore contigs not in the CALLs file
	}
	if (($prev_qid ne $qid) && ($prev_qid ne "")) {
	    my @calls_line = split(/\t/,$insertion_events{$prev_qid});
	    my $category = $calls_line[CTYPE];
	    $calls_line[CTYPE] =~ s/\s+/_/g;
	    (my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	    if ($found) {
		if (!defined $foreign_plasmids{$contig_id}) {
		    # only do this if this was not already called a foreign plasmid
		    $engineering_found = 1;
		    print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
		    print $out_short "$sample_name\t$species_name\tclear\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
		    if (!defined $event_counts{$species_name}) {
			$event_counts{$species_name} = {};
			$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
		    }
		    if (!defined $event_counts{$species_name}->{$sample_name}) {
			$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
		    }
		    $event_counts{$species_name}->{$sample_name}->[TINSC]++;
		    $event_counts{$species_name}->{"Total"}->[TINSC]++;
		    $event_counts_total[TINSC]++;
		    if ($calls_line[CILEN] >= 20) {
			print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] clear inserted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
		    }
		    if ($calls_line[CDLEN] >= 20) {
			print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] clear deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
		    }
		    if (length($calls_line[CFIVEP]) >= 20) {
			print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] clear 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
		    }
		    if (length($calls_line[CTHREEP]) >= 20) {
			print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] clear 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
		    }
		}
	    } elsif (!defined $foreign_plasmids{$contig_id}) {
		print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
		print $out_short "$sample_name\t$species_name\tmaybe\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
		if (!defined $event_counts{$species_name}) {
		    $event_counts{$species_name} = {};
		    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
		}
		if (!defined $event_counts{$species_name}->{$sample_name}) {
		    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
		}
		$event_counts{$species_name}->{$sample_name}->[TINSM]++;
		$event_counts{$species_name}->{"Total"}->[TINSM]++;
		$event_counts_total[TINSM]++;
		if ($calls_line[CILEN] >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] maybe inserted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
		}
		if ($calls_line[CDLEN] >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] maybe deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
		}
		if (length($calls_line[CFIVEP]) >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] maybe 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
		}
		if (length($calls_line[CTHREEP]) >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] maybe 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
		}
	    } else {
		print STDERR "Not outputting previoulsy output as a foreign plasmid: $calls_line[CQID]\n"
	    }
	    delete $insertion_events{$prev_qid};
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
	if (($stitle =~ /vector/i) || ($stitle =~ /construct/i) || ($stitle =~ /synthetic/i) || ($stitle =~ /marker plasmid/i)) {
	    $found = 1;
	}
	if ($bitscore > $max_bitscore) {
	    $max_bitscore = $bitscore;
	    $max_title = $stitle;
	    $max_sid = $sid;
	    $max_length = $qlen;
	}
    }
    if ($prev_qid ne "") {
	my @calls_line = split(/\t/,$insertion_events{$prev_qid});
	my $category = $calls_line[CTYPE];
	$calls_line[CTYPE] =~ s/\s+/_/g;
	(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	if ($found) {
	    if (!defined $foreign_plasmids{$contig_id}) {
		# only do this if this was not already called a foreign plasmid
		$engineering_found = 1;
		print $out_clear "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
		print $out_short "$sample_name\t$species_name\tclear\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
		if (!defined $event_counts{$species_name}) {
		    $event_counts{$species_name} = {};
		    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
		}
		if (!defined $event_counts{$species_name}->{$sample_name}) {
		    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
		}
		$event_counts{$species_name}->{$sample_name}->[TINSC]++;
		$event_counts{$species_name}->{"Total"}->[TINSC]++;
		$event_counts_total[TINSC]++;
		if ($calls_line[CILEN] >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] clear inserted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
		}
		if ($calls_line[CDLEN] >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] clear deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
		}
		if (length($calls_line[CFIVEP]) >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] clear 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
		}
		if (length($calls_line[CTHREEP]) >= 20) {
		    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] clear 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
		}
	    }
	} elsif (!defined $foreign_plasmids{$contig_id}) {
	    print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG: BLAST: $results\t$species_name\t$max_title\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	    print $out_short "$sample_name\t$species_name\tmaybe\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
	    if (!defined $event_counts{$species_name}) {
		$event_counts{$species_name} = {};
		$event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	    }
	    if (!defined $event_counts{$species_name}->{$sample_name}) {
		$event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	    }
	    $event_counts{$species_name}->{$sample_name}->[TINSM]++;
	    $event_counts{$species_name}->{"Total"}->[TINSM]++;
	    $event_counts_total[TINSM]++;
	    if ($calls_line[CILEN] >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] maybe inserted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
	    }
	    if ($calls_line[CDLEN] >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] maybe deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
	    }
	    if (length($calls_line[CFIVEP]) >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] maybe 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
	    }
	    if (length($calls_line[CTHREEP]) >= 20) {
		print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] maybe 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
	    }
	} else {
	    print STDERR "Not outputting previoulsy output as a foreign plasmid: $calls_line[CQID]\n"
	}
	delete $insertion_events{$prev_qid};
    }
    close($calls_btab_file);
    foreach my $key (keys %insertion_events) {
	my @calls_line = split(/\t/,$insertion_events{$key});
	my $category = $calls_line[CTYPE];
	$calls_line[CTYPE] =~ s/\s+/_/g;
	(my $contig_id, my $details) = split(/_DIV/, $calls_line[CQID]);
	print $out_maybe "$sample_name\tyes\t$calls_line[CSID]\tyes\tinserted sequence : $calls_line[CINSERTED] : deleted sequence : $calls_line[CDELETED]\t\tcomparison to PGG\t$species_name\tunknown insertion\tContig $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]:5' flank $calls_line[CFIVEP]:3' flank $calls_line[CTHREEP]: Reference $calls_line[CSID]:coordinates $calls_line[CSSTART],$calls_line[CSEND]\tinsertion : $calls_line[CILEN] : deletion : $calls_line[CDLEN]\t$contig_id : $calls_line[CSID]\t\t\t\n";
	print $out_short "$sample_name\t$species_name\tmaybe\t$category\t$contig_id\t$calls_line[CQSTART]\t$calls_line[CQEND]\n";
	if (!defined $event_counts{$species_name}) {
	    $event_counts{$species_name} = {};
	    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	}
	if (!defined $event_counts{$species_name}->{$sample_name}) {
	    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	}
	$event_counts{$species_name}->{$sample_name}->[TINSM]++;
	$event_counts{$species_name}->{"Total"}->[TINSM]++;
	$event_counts_total[TINSM]++;
	if ($calls_line[CILEN] >= 20) {
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_ins_$calls_line[CTYPE] maybe inserted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CINSERTED]);
	}
	if ($calls_line[CDLEN] >= 20) {
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_del_$calls_line[CTYPE] maybe deleted sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CDELETED]);
	}
	if (length($calls_line[CFIVEP]) >= 20) {
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_5p_$calls_line[CTYPE] maybe 5 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CFIVEP]);
	}
	if (length($calls_line[CTHREEP]) >= 20) {
	    print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$calls_line[CQSTART]_$calls_line[CQEND]_3p_$calls_line[CTYPE] maybe 3 prime flank sequence $calls_line[CQID] $contig_id:coordinates $calls_line[CQSTART],$calls_line[CQEND]", $calls_line[CTHREEP]);
	}
    }
    if (!$engineering_found) {
	print $out_clear "$sample_name\tyes\tPGG\tno\tNA\t\tcomparison to PGG\t$species_name\tNA\tNA\tNA\tNA\tNA\tNA\t\n";
    }
    
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
	print $out_short "$sample_name\t$species_name\tmaybe\tSTOP CODON\t$contig_id\t$stop_codon_coord\t$stop_codon_coord\n";
	if (!defined $event_counts{$species_name}) {
	    $event_counts{$species_name} = {};
	    $event_counts{$species_name}->{"Total"} = [0,0,0,0,0,0,0,0];
	}
	if (!defined $event_counts{$species_name}->{$sample_name}) {
	    $event_counts{$species_name}->{$sample_name} = [0,0,0,0,0,0,0,0];
	}
	$event_counts{$species_name}->{$sample_name}->[TSC]++;
	$event_counts{$species_name}->{"Total"}->[TSC]++;
	$event_counts_total[TSC]++;
	print_fasta($out_fasta, "$sample_name" . "$species_name_ns" . "_$contig_id" . "_$stop_codon_coord $type_stop_codon $feat_name $frame", $feat_seq);
    }
    close($stop_codons_file);
}

print $out_summary "Species\tSample\tClear Insertions\tClear Foreign Plasmids\tDeletions\tPossible Insertions\tPossible Foreign Plasmids\tTandem Duplications\tMutations\tStop Codons\n";

my $sort_by_name_total = sub { # sort by name but put Total last
    if ($a eq "Total") {
	return (1);
    } elsif ($b eq "Total") {
	return (-1);
    } else {
	return ($a cmp $b);
    }
};

my @sorted_species = sort $sort_by_name_total (keys %event_counts);

foreach my $species (@sorted_species) {
    my @sorted_samples = sort $sort_by_name_total (keys %{ $event_counts{$species} });
    foreach my $sample (@sorted_samples) {
	my $counts_string = join('\t', @{ $event_counts{$species}->{$sample} });
	print $out_summary "$species\t$sample\t$counts_string\n";
    }
}

my $total_counts_string = join('\t', @event_counts_total);
print $out_summary "All\tTotal\t$total_counts_string\n"; 

close($out_clear);
close($out_stop);
close($out_maybe);
close($out_fasta);
close($out_short);
close($out_summary);

exit(0);
