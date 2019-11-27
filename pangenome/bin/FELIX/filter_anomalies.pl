#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

##### Attribute Compare Script

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;
use List::Util qw[min max];

my @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
my @pgg_blast_results = ();  # These are the blast results of the query sequences against the PGG genomes
my @ordered = ();      # Same as above but sorted by CONTIG then START
my %contig_len = ();   # key = contig id, value = length of contig
my %contigs = ();      # key = contig id, value = sequence of contig
my %query_seqs = ();   # key = query id, value = sequence of query
my %pggdb_contig_len = ();   # key = contig id, value = length of contig
my %pggdb_contigs = ();      # key = contig id, value = sequence of contig
my %seen_contig = ();  # key = contig id, value = 1 (placeholder to show we've seen this contig)
my %max_bitscore = (); # key1 = query id, key2 = subject id, array = [MB5P,MB3P,MBFULL], value = max bitscore : pgg_blast_results index
my @categories = ("identical","gapped","divergent","conserved"); # literals used for output in ranges file

# CONSTANTS #
use constant MBBIT => 0;
use constant MBIND => 1;
use constant MBSEC => 2;
use constant MB5P => 0;
use constant MB3P => 1;
use constant MBFULL => 2;
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
use constant IDENTICAL => 0;
use constant GAPPED => 1;
use constant DIVERGED => 2;
use constant MAYBE => 3;
use constant CONTIG => 0;
use constant LOCUS => 1;
use constant START => 2;
use constant STOP => 3;
use constant TYPE => 4;
use constant GENOME => 5;
use constant INVERT => 6;
use constant LENGTH => 7;
use constant DELETE => 8;
use constant CATEGORY => 9;
# END CONSTANTS #

GetOptions('genomes=s' => \my $genomes,
	'help' => \my $help,
	'debug' => \my $debug,
	'nrdb=s' => \my $nrdb,
	'pggdb=s' => \my $PGGdb,
	'engdb=s' => \my $engdb);

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genomes=s' => genomes,
	'help' => help,
	'debug' => debug,
	'nrdb=s' => nrdb,
	'pggdb=s' => PGGdb,
	'engdb=s' => engdb);
_EOB_
    exit(0);
}

# subroutine to print contig segments out in fasta format
sub print_fasta { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file_handle, $contig_name, $seq, $beg, $end) = @_;
    $beg--;
    $end--;
    print $file_handle ">$contig_name\n";
    my $pos;
    my $seq_len = ($end - $beg) + 1;
    $query_seqs{$contig_name} = substr($seq, $beg, $seq_len);
    for ( $pos = $beg ; $seq_len > 60 ; $pos += 60 ) {
	print $file_handle substr($seq, $pos, 60), "\n";
	$seq_len -= 60;
    }
    print $file_handle substr($seq, $pos, $seq_len), "\n";
    return;
}

# read in contigs from genome file
my $pggdbfile;
unless (open ($pggdbfile, "<", $PGGdb) )  {
    die ("cannot open PGG database file: $PGGdb!\n");
}
my $pggdb_save_input_separator = $/;
my $pggdb_line;
$/="\n>";
while ($pggdb_line = <$pggdbfile>) {
    (my $title, my $sequence) = split(/\n/, $pggdb_line, 2); # split the header line and sequence (very cool)
    my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
    my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
    $id =~ s/>\s*//; # remove leading > and spaces
    #$id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
    $pggdb_contigs{$id} = $sequence;
    $pggdb_contig_len{$id} = length($sequence);
    #print STDERR "\n";
    $title = ""; # clear the title for the next contig
    $sequence = ""; #clear out the sequence for the next contig
}
$/ = $pggdb_save_input_separator; # restore the input separator
close ($pggdbfile);

# read file which specifies the output file prefix, assembly fasta file, and anomalies file
open (my $infile, "<", $genomes) || die ("ERROR: cannot open file $genomes\n");
while (my $line = <$infile>)  {
    # clear data structures for the next genome
    @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
    @ordered = ();      # Same as above but sorted by CONTIG then START
    %contigs = ();      # key = contig id, value = sequence of contig
    %contig_len = ();   # key = contig id, value = length of contig
    %seen_contig = ();  # key = contig id, value = 1 (placeholder to show we've seen this contig)
    @pgg_blast_results = ();  # These are the blast results of the query sequences against the PGG genomes
    %query_seqs = ();   # key = query id, value = sequence of query
    %max_bitscore = (); # key1 = query id, key2 = subject id, array = [MB5P,MB3P,MBFULL], value = max bitscore : pgg_blast_results index
    
    chomp $line;
    (my $output, my $genome, my $anomalies) = split(/\t/, $line);  # split the scalar $line on tab

    # read in contigs from genome file
    my $contigfile;
    unless (open ($contigfile, "<", $genome) )  {
	die ("cannot open genome file: $genome!\n");
    }
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    while ($line = <$contigfile>) {
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$id =~ s/>\s*//; # remove leading > and spaces
	$id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	my $contig_length = length($sequence);
	#print STDERR "$id\t$contig_length";
	if (($contig_length >= 1000) || ($id !~ /_cov_[01]$/)) { #only store contigs whose depth of coverage is > 1 or whose length is >= 1000
	    $contigs{$id} = $sequence;
	    $contig_len{$id} = $contig_length;
	    #print STDERR "\t$id";
	}
	#print STDERR "\n";
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($contigfile);

    # read in anomalies file
    open(ANOM_FILE, "<", $anomalies) || die ("Couldn't open $anomalies\n");
    my $count = 0;
    my $first = 1;
    while (my $line = <ANOM_FILE>) {
	if ($first) {
	    $first = 0;
	    next; # first line is a header
	}
	chomp($line);
	my @split_line = split(/\t/,$line);
	if ($split_line[5] <= 0) {
	    next; # ignore negative or zero length edge anomalies
	}
	my $category = $split_line[2];
	if (($category =~ /^divergent/) || ($category =~ /^very/) || ($category eq "uniq_edge")) {
	    $category = DIVERGED;
	} elsif ($category =~ /^identical/) {
	    $category = IDENTICAL;
	} elsif ($category =~ /^gapped/) {
	    $category = GAPPED;
	} elsif ($category =~ /^uniq_.*_allele$/) {
	    $category = MAYBE;
	} else {
	    next; # ignore other types
	}
	$annotations[$count][CATEGORY] = $category;       # category
	$annotations[$count][TYPE] = $split_line[2];       # type
	my $contig = $split_line[1];
	$contig =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	$annotations[$count][CONTIG] = $contig;     # contig
	$annotations[$count][LOCUS] = $split_line[6];      # locus_id
	$annotations[$count][GENOME] = $split_line[0];     # genome
	$annotations[$count][LENGTH] = $split_line[5];     # length
	$annotations[$count][DELETE] = 0;     # subsumed by another segment
	if ($split_line[3] < $split_line[4]) {
	    $annotations[$count][START] = $split_line[3]; # start
	    $annotations[$count][STOP] = $split_line[4]; # stop
	    $annotations[$count][INVERT] = 0;              # invert
	} else {
	    $annotations[$count][START] = $split_line[4]; # start
	    $annotations[$count][STOP] = $split_line[3]; # stop
	    $annotations[$count][INVERT] = 1;              # invert
	}
	if ((($annotations[$count][STOP] - $annotations[$count][START]) + 1) != $annotations[$count][LENGTH]) {
	    next; # for now ignore annotations going around the end of the contig - also gets rid of negative and zero length edges
	}
	$count++;
    }
    $count--;
    $#annotations = $count; #make sure that if last annotation was not supposed to be saved that it is tuncated off
    close(ANOM_FILE);
    
    if ($debug) {
	print "DEBUG***annotations\n";
	for (my $j=0; $j < @annotations; $j++) {
	    print ("$j: $annotations[$j][CONTIG]\t$annotations[$j][LOCUS]\t$annotations[$j][START]\t$annotations[$j][STOP]\t$annotations[$j][TYPE]\t$ordered[$j][CATEGORY]\t$annotations[$j][GENOME]\t$annotations[$j][INVERT]\n");
	}
    }

    # sort attribute files by contig then by start, store in ordered data-structure. 

    @ordered = sort { $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[CATEGORY] <=> $b->[CATEGORY] } @annotations; # sort on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***ordered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][CATEGORY]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }

    # eliminate overlaps between annotated segments by truncation - delete those that are contained by others - order of precedence: IDENTICAL, DIVERGED, MAYBE, GAPPED
    for (my $i=0; $i < @ordered; $i++) {
	if ($ordered[$i][DELETE]) {
	    next;
	}
	for (my $j=$i+1; $j < @ordered; $j++) {
	    if ($ordered[$j][DELETE]) {
		last;
	    }
	    if ($ordered[$i][CONTIG] lt $ordered[$j][CONTIG]) {
		last;
	    }
	    if ($ordered[$i][CONTIG] gt $ordered[$j][CONTIG]) {
		die ("ERROR: Annotations are not properly sorted!\n");
	    }
	    if ($ordered[$i][STOP] < $ordered[$j][START]) {
		last;
	    }
	    if ($ordered[$i][STOP] < $ordered[$j][STOP]) {
		if ($ordered[$i][CATEGORY] <= $ordered[$j][CATEGORY]) {
		    $ordered[$j][START] = $ordered[$i][STOP] + 1;
		} else {
		    $ordered[$i][STOP] = $ordered[$j][START] - 1;
		    if ($ordered[$i][STOP] < $ordered[$i][START]) {
			$ordered[$i][DELETE] = 1;
			last;
		    }
		}
	    } else {
		if ($ordered[$i][CATEGORY] <= $ordered[$j][CATEGORY]) {
		    $ordered[$j][DELETE] = 1;
		} else {
		    $ordered[$i][STOP] = $ordered[$j][START] - 1;
		    if ($ordered[$i][STOP] < $ordered[$i][START]) {
			$ordered[$i][DELETE] = 1;
			last;
		    }
		}
	    }
	}
    }

    # open file for extracted interesting sequences
    my $out_fasta_seqs = $output . "_QUERY_SEQS.fasta";
    my $file_fasta_seqs;
    unless (open ($file_fasta_seqs, ">", $out_fasta_seqs) )  {
	die ("cannot open file $out_fasta_seqs!\n");
    }

    # open file for extracted interesting sequences
    my $out_ranges = $output . "_ranges.txt";
    my $file_ranges;
    unless (open ($file_ranges, ">", $out_ranges) )  {
	die ("cannot open file $out_ranges!\n");
    }

    # resort attribute files by deleted or not, then by contig, then by start, store in ordered data-structure. 

    @ordered = sort { $a->[DELETE] <=> $b->[DELETE] || $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[CATEGORY] <=> $b->[CATEGORY] } @ordered; # sort on deleted or not, then on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***filtered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][CATEGORY]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }

    # output segments which are diverged, or at the ends of contigs
    my $range_beg = 0;
    my $range_end = 0;
    my $cur_beg;
    my $cur_end;
    my $cur_category;
    my $cur_contig;
    my $prev_beg = 1;
    my $prev_end = 0;
    my $prev_category = -1;
    my $prev_contig = "";
    my $last_output = 0;
    my $beg_diverged = 0;
    my $end_diverged = 0;
    my $cur_type;
    my $cur_locus;
    my $diverged_type = "";
    for (my $i=0; $i < @ordered; $i++) {
	if ($ordered[$i][DELETE]) {
	    last;
	}
	$cur_contig = $ordered[$i][CONTIG];
	if (!defined $contigs{$cur_contig}) {
	    #print STDERR "triaged: $cur_contig\n";
	    next; # skip triaged contigs
	}
	$cur_beg = $ordered[$i][START];
	$cur_end = $ordered[$i][STOP];
	$cur_type = $ordered[$i][TYPE];
	$cur_locus = $ordered[$i][LOCUS];
	$cur_category = $ordered[$i][CATEGORY];
	$seen_contig{$cur_contig} = 1;
	if ($prev_contig ne $cur_contig) { # first segment for this contig
	    if (($prev_contig ne "") && ($beg_diverged != 0)) {
		&print_fasta($file_fasta_seqs, ($prev_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type), $contigs{$prev_contig}, $beg_diverged, $end_diverged);
		$beg_diverged = 0;
		$end_diverged = 0;
		$diverged_type = "";
	    }
	    if (($prev_contig ne "") && (($contig_len{$prev_contig} - $prev_end) > 20)) { # include unannotated contig ends > 20 bp
		my $tmp_seq = substr($contigs{$prev_contig}, $prev_end, ($contig_len{$prev_contig} - $prev_end));
		if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		    &print_fasta($file_fasta_seqs, ($prev_contig . "_END_" . ($prev_end + 1) . "_" . $contig_len{$prev_contig}), $contigs{$prev_contig}, ($prev_end + 1), $contig_len{$prev_contig});
		}
	    }
	    if ($cur_beg > 20) { # include unannotated contig ends > 20 bp
		my $tmp_seq = substr($contigs{$cur_contig}, 0, ($cur_beg - 1));
		if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		    &print_fasta($file_fasta_seqs, ($cur_contig . "_BEG_1_" . ($cur_beg - 1)), $contigs{$cur_contig}, 1, ($cur_beg - 1));
		}
	    }
	    if (($prev_contig ne "") && ($range_beg != 0)) {
		print $file_ranges "$prev_contig\t$range_beg\t$range_end\t$categories[$prev_category]\n";
	    }
	    if (($prev_contig ne "") && (($contig_len{$prev_contig} - $range_end) > 0)) { # include unannotated contig ends > 0 bp
		print $file_ranges "$prev_contig\t", ($range_end + 1), "\t$contig_len{$prev_contig}\tunannotated\n";
	    }
	    if ($cur_beg > 1) { # include unannotated contig ends > 0 bp
		print $file_ranges "$cur_contig\t1\t", ($cur_beg - 1), "\tunannotated\n";
	    }
	    $range_beg = 0;
	    $range_end = 0;
	}
	
	if ($range_beg == 0) {
	    $range_beg =  $cur_beg;
	    $range_end =  $cur_end;
	} elsif ($cur_category == $prev_category) {
	    $range_end = $cur_end;
	} else {
	    print $file_ranges "$cur_contig\t$range_beg\t$range_end\t$categories[$prev_category]\n";
	    $range_beg =  $cur_beg;
	    $range_end =  $cur_end;
	}
	
	if (($cur_category == IDENTICAL) || ($cur_category == MAYBE) || ($cur_category == GAPPED)) {
	    if ($beg_diverged != 0) {
		&print_fasta($file_fasta_seqs, ($cur_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type), $contigs{$cur_contig}, $beg_diverged, $end_diverged);
		$beg_diverged = 0;
		$end_diverged = 0;
		$diverged_type = "";
	    }
	} elsif ($cur_category == DIVERGED) {
	    if ($beg_diverged == 0) {
		$beg_diverged = ($cur_beg > 100) ? ($cur_beg - 100) : 1; # instead of $cur_beg to include context around the diverged region
		$diverged_type = $cur_type . "_" . $cur_locus . "_" . $cur_beg . "_" . $cur_end;
	    } else {
		$diverged_type .= "_" . $cur_type . "_" . $cur_locus . "_" . $cur_beg . "_" . $cur_end;
	    }
	    $end_diverged = (($contig_len{$cur_contig} - $cur_end) > 100) ? ($cur_end + 100) : $contig_len{$cur_contig}; #instead of #cur_end to include context around the diverged region
	} else {
	    die ("ERROR: Unexpected type found: $cur_category\n");
	}
	$prev_beg = $cur_beg;
	$prev_end = $cur_end;
	$prev_category = $cur_category;
	$prev_contig = $cur_contig;
    }
    if (($prev_contig ne "") && ($beg_diverged != 0)) {
	&print_fasta($file_fasta_seqs, ($prev_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type), $contigs{$prev_contig}, $beg_diverged, $end_diverged);
    }
    if (($prev_contig ne "") && (($contig_len{$prev_contig} - $prev_end) > 20)) { # include unannotated contig ends > 20 bp
	my $tmp_seq = substr($contigs{$prev_contig}, $prev_end, ($contig_len{$prev_contig} - $prev_end));
	if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
	    &print_fasta($file_fasta_seqs, ($prev_contig . "_END_" . ($prev_end + 1) . "_" . $contig_len{$prev_contig}), $contigs{$prev_contig}, ($prev_end + 1), $contig_len{$prev_contig});
	}
    }
    if (($prev_contig ne "") && ($range_beg != 0)) {
	print $file_ranges "$prev_contig\t$range_beg\t$range_end\t$categories[$prev_category]\n";
    }
    if (($prev_contig ne "") && (($contig_len{$prev_contig} - $range_end) > 0)) { # include unannotated contig ends > 0 bp
	print $file_ranges "$prev_contig\t", ($range_end + 1), "\t$contig_len{$prev_contig}\tunannotated\n";
    }

    # output entire contigs with no annotation
    foreach my $id (keys %contigs)  { # go through all contigs
	if (!defined $seen_contig{$id}) {
	    #print STDERR "not seen: $id\t$contig_len{$id}\n";
	    my $tmp_seq = substr($contigs{$id}, 0, $contig_len{$id});
	    if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		&print_fasta($file_fasta_seqs, ($id . "_WHOLE_1_" . $contig_len{$id}), $contigs{$id}, 1, $contig_len{$id});
	    }
	}
    }

    close($file_fasta_seqs);
    close($file_ranges);

    my $out_predictions = $output . "_CALLS";
    my $out_genome = $output . "_GENOME";
    my $out_genome_nsq = $out_genome . ".nsq";
    my $out_genome_nhr = $out_genome . ".nhr";
    my $out_genome_nin = $out_genome . ".nin";
    my $out_genome_blast = $out_genome . ".btab";
    my $out_PGG_blast = $output . "_PGG.btab";
    my $out_nrdb_blast = $output . "_NRDB.btab";
    my $out_engdb_blast = $output . "_ENGDB.btab";
    my $combined_btab = $output . "_COMBINED.btab";
    my $out_features = $output . "_FEATURES";
    my $cmd;
    `cp $genome $out_genome`;
    $cmd = 'export LD_LIBRARY_PATH=/usr/local/packages/glibc-2.14/lib:$LD_LIBRARY_PATH; ' . "/usr/local/projdata/99999/IFX/CommonDB/ncbi-blast+/ncbi-blast-2.9.0+/bin/makeblastdb -in $out_genome -dbtype nucl";
    `$cmd`;
    $cmd = 'export LD_LIBRARY_PATH=/usr/local/packages/glibc-2.14/lib:$LD_LIBRARY_PATH; ' . "/usr/local/projdata/99999/IFX/CommonDB/ncbi-blast+/ncbi-blast-2.9.0+/bin/blastn -query $out_fasta_seqs -db $out_genome -out $out_genome_blast -task blastn -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 3";
    `$cmd`;
    `rm $out_genome $out_genome_nsq $out_genome_nin $out_genome_nhr`;
    #$cmd = 'export LD_LIBRARY_PATH=/usr/local/packages/glibc-2.14/lib:$LD_LIBRARY_PATH; ' . "/usr/local/projdata/99999/IFX/CommonDB/ncbi-blast+/ncbi-blast-2.9.0+/bin/blastn -query $out_fasta_seqs -db $engdb -out $out_engdb_blast -task blastn -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 2";
    #`$cmd`;
    $cmd = 'export LD_LIBRARY_PATH=/usr/local/packages/glibc-2.14/lib:$LD_LIBRARY_PATH; ' . "/usr/local/projdata/99999/IFX/CommonDB/ncbi-blast+/ncbi-blast-2.9.0+/bin/blastn -query $out_fasta_seqs -db $PGGdb -out $out_PGG_blast -task blastn -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\"";
    `$cmd`;
    # read in anomalies file
    open(PGG_BLAST_FILE, "<", $out_PGG_blast) || die ("Couldn't open $out_PGG_blast for reading\n");
    $count = 0;
    while (my $line = <PGG_BLAST_FILE>) {
	chomp($line);
	my @split_line = split(/\t/,$line);
	my $qid = $pgg_blast_results[$count][QSEQID] = $split_line[0];         # query id
	my $sid = $pgg_blast_results[$count][SSEQID] = $split_line[1];         # subject id
	my $pid = $pgg_blast_results[$count][PIDENT] = $split_line[2];         # percent identity
	my $qstart = $pgg_blast_results[$count][QSTART] = $split_line[3];      # query start
	my $qend = $pgg_blast_results[$count][QEND] = $split_line[4];          # query end
	my $qlen = $pgg_blast_results[$count][QLEN] = $split_line[5];          # query length
	my $sstart = $pgg_blast_results[$count][SSTART] = $split_line[6];      # subject start
	my $send = $pgg_blast_results[$count][SEND] = $split_line[7];          # subject end
	my $slen = $pgg_blast_results[$count][SLEN] = $split_line[8];          # subject length
	my $evalue = $pgg_blast_results[$count][EVALUE] = $split_line[9];      # evalue
	my $bitscore = $pgg_blast_results[$count][BITSCORE] = $split_line[10]; # bitscore
	my $stitle = $pgg_blast_results[$count][STITLE] = $split_line[11];     # subject title
	if (($pid >= 95) && (($qstart <= 6) || ($qend >= ($qlen - 5)))) {
	    if (($qstart <= 6) && ($qend >= ($qlen - 5))) {
		if (defined $max_bitscore{$qid}{$sid}) {
		    $max_bitscore{$qid}{$sid}[MB5P][MBBIT] = 0;
		    $max_bitscore{$qid}{$sid}[MB5P][MBSEC] = 0;
		    $max_bitscore{$qid}{$sid}[MB5P][MBIND] = -1;
		    $max_bitscore{$qid}{$sid}[MB3P][MBBIT] = 0;
		    $max_bitscore{$qid}{$sid}[MB3P][MBSEC] = 0;
		    $max_bitscore{$qid}{$sid}[MB3P][MBIND] = -1;
		    if ($bitscore > $max_bitscore{$qid}{$sid}[MBFULL][MBBIT]) {
			$max_bitscore{$qid}{$sid}[MBFULL][MBSEC] = $max_bitscore{$qid}{$sid}[MBFULL][MBBIT];
			$max_bitscore{$qid}{$sid}[MBFULL][MBBIT] = $bitscore;
			$max_bitscore{$qid}{$sid}[MBFULL][MBIND] = $count;
		    }
		} else {
		    $max_bitscore{$qid}{$sid} = [ [ 0, -1, 0 ], [ 0, -1, 0 ], [ $bitscore, $count, 0 ] ];
		}
	    } elsif ($qstart <= 6) {
		if (defined $max_bitscore{$qid}{$sid}) {
		    if (($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] == 0) && ($bitscore > $max_bitscore{$qid}{$sid}[MB5P][MBBIT])) {
			$max_bitscore{$qid}{$sid}[MB5P][MBSEC] = $max_bitscore{$qid}{$sid}[MB5P][MBBIT];
			$max_bitscore{$qid}{$sid}[MB5P][MBBIT] = $bitscore;
			$max_bitscore{$qid}{$sid}[MB5P][MBIND] = $count;
		    }
		} else {
		    $max_bitscore{$qid}{$sid} = [ [ $bitscore, $count, 0 ], [ 0, -1, 0 ], [ 0, -1, 0 ] ];
		}
	    } else {
		if (defined $max_bitscore{$qid}{$sid}) {
		    if (($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] == 0) && ($bitscore > $max_bitscore{$qid}{$sid}[MB3P][MBBIT])) {
			$max_bitscore{$qid}{$sid}[MB3P][MBSEC] = $max_bitscore{$qid}{$sid}[MB3P][MBBIT];
			$max_bitscore{$qid}{$sid}[MB3P][MBBIT] = $bitscore;
			$max_bitscore{$qid}{$sid}[MB3P][MBIND] = $count;
		    }
		} else {
		    $max_bitscore{$qid}{$sid} = [ [ 0, -1, 0 ], [ $bitscore, $count, 0 ], [ 0, -1, 0 ] ];
		}
	    }
	    $count++;
	}
    }
    $count--;
    $#pgg_blast_results = $count; #make sure that if last blast result was not supposed to be saved that it is tuncated off
    close(PGG_BLAST_FILE);
    open(PGG_BLAST_FILE, ">", $out_PGG_blast) || die ("Couldn't open $out_PGG_blast for writing\n");
    open(CALLS_FILE, ">", $out_predictions) || die ("Couldn't open $out_predictions for writing\n");
    print CALLS_FILE "Sample\tType\tQuery ID\t5pflank\t3pflank\tDeleted\tDeletion Length\tInserted\tInsertion Length\tSubject ID\t% Identity\tSubject Start\tSubject End\tQuery Start\tQuery End\n";
    my $mutations = 0;
    my $deletions = 0;
    my $insertions = 0;
    my $tandem_duplications = 0;
    foreach my $qid (keys %max_bitscore)  { # go through all query sequence blast matches
	my $full_max = 0;
	my $sum_max = 0;
	my $full_index = -1;
	my $fivep_index = -1;
	my $threep_index = -1;
	my $max_sid = "";
	my $repeat = 0;
	foreach my $sid (keys %{ $max_bitscore{$qid} })  { # go through all subject sequence blast matches looking for max
	    if ($full_max > 0) {
		if ($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] > $full_max) {
		    $full_max = $max_bitscore{$qid}{$sid}[MBFULL][MBBIT];
		    $full_index = $max_bitscore{$qid}{$sid}[MBFULL][MBIND];
		    $max_sid = $sid;
		    if ($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] == $max_bitscore{$qid}{$sid}[MBFULL][MBSEC]) {
			$repeat = 1;
		    } else {
			$repeat = 0;
		    }
		}
	    } else {
		if ($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] > $full_max) {
		    $full_max = $max_bitscore{$qid}{$sid}[MBFULL][MBBIT];
		    $full_index = $max_bitscore{$qid}{$sid}[MBFULL][MBIND];
		    $max_sid = $sid;
		    if ($max_bitscore{$qid}{$sid}[MBFULL][MBBIT] == $max_bitscore{$qid}{$sid}[MBFULL][MBSEC]) {
			$repeat = 1;
		    } else {
			$repeat = 0;
		    }
		    $sum_max = 0;
		    $fivep_index = -1;
		    $threep_index = -1;
		} elsif (($max_bitscore{$qid}{$sid}[MB5P][MBBIT] + $max_bitscore{$qid}{$sid}[MB3P][MBBIT]) > $sum_max) {
		    $sum_max = $max_bitscore{$qid}{$sid}[MB5P][MBBIT] + $max_bitscore{$qid}{$sid}[MB3P][MBBIT];
		    $fivep_index = $max_bitscore{$qid}{$sid}[MB5P][MBIND];
		    $threep_index = $max_bitscore{$qid}{$sid}[MB3P][MBIND];
		    $max_sid = $sid;
		    if (($max_bitscore{$qid}{$sid}[MB5P][MBBIT] == $max_bitscore{$qid}{$sid}[MB5P][MBSEC]) || ($max_bitscore{$qid}{$sid}[MB3P][MBBIT] == $max_bitscore{$qid}{$sid}[MB3P][MBSEC])) {
			$repeat = 1;
		    } else {
			$repeat = 0;
		    }
		}
	    }
	}
	#print STDERR "$full_max:$full_index:$sum_max:$fivep_index:$threep_index:$#pgg_blast_results:$max_sid\n";
	if ($max_sid eq "") {
	    next; #no good matches for this query sequence
	}
	if (($full_index > 0) && ($max_sid ne $pgg_blast_results[$full_index][SSEQID])) {
	    die ("ERROR: for full_index $full_index subject id $max_sid not the same as $pgg_blast_results[$full_index][SSEQID]\n");
	} elsif (($fivep_index > 0) && ($max_sid ne $pgg_blast_results[$fivep_index][SSEQID])) {
	    die ("ERROR: for fivep_index $fivep_index subject id $max_sid not the same as $pgg_blast_results[$fivep_index][SSEQID]\n");
	} elsif (($threep_index > 0) && ($max_sid ne $pgg_blast_results[$threep_index][SSEQID])) {
	    die ("ERROR: for threep_index $threep_index subject id $max_sid not the same as $pgg_blast_results[$threep_index][SSEQID]\n");
	}
	if ($full_max > 0) {# possible mutation
	    if (($full_index > $#pgg_blast_results) || ($full_index < 0)) {
		die ("ERROR: full_index out of range $full_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$full_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    my $pid = $pgg_blast_results[$full_index][PIDENT];
	    if (($pid < 99.5) && ($qid !~ /_WHOLE_/) && ($qid !~ /_BEG_/) && ($qid !~ /_END_/)) {
		my $cur_sid = $pgg_blast_results[$full_index][SSEQID];
		my $qstart = $pgg_blast_results[$full_index][QSTART];
		my $qend = $pgg_blast_results[$full_index][QEND];
		my $sstart = $pgg_blast_results[$full_index][SSTART];
		my $send = $pgg_blast_results[$full_index][SEND];
		my $revcomp = ($sstart > $send);
		my $insertion = "";
		$insertion = substr($query_seqs{$qid}, ($qstart - 1), (($qend - $qstart) + 1));
		my $deletion = "";
		if ($revcomp) {
		    my $tmp_seq = substr($pggdb_contigs{$cur_sid}, ($send - 1), (($sstart - $send) + 1));
		    $deletion = reverse($tmp_seq);
		    $deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		} else {
		    $deletion = substr($pggdb_contigs{$cur_sid}, ($sstart - 1), (($send - $sstart) + 1));
		}
		$mutations++;
		my $del_len = length($deletion);
		my $ins_len = length($insertion);
		print CALLS_FILE "$output\tMUTATION\t$qid\t\t\t$deletion\t$del_len\t$insertion\t$ins_len\t$cur_sid\t$pid\t$sstart\t$send\t$qstart\t$qend\n";
	    }
	} elsif (($fivep_index >= 0) && ($threep_index >= 0)) {# possible deletion or insertion/replacement
	    if (($fivep_index > $#pgg_blast_results) || ($fivep_index < 0)) {
		die ("ERROR: fivep_index out of range $fivep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    if (($threep_index > $#pgg_blast_results) || ($threep_index < 0)) {
		die ("ERROR: threep_index out of range $threep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$fivep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    $btab = join("\t", @{ $pgg_blast_results[$threep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    my $cur_sid = $pgg_blast_results[$fivep_index][SSEQID];
	    my $pid = ($pgg_blast_results[$fivep_index][PIDENT] + $pgg_blast_results[$threep_index][PIDENT]) / 2;
	    my $revcomp5p = ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$fivep_index][SEND]);
	    my $revcomp3p = ($pgg_blast_results[$threep_index][SSTART] > $pgg_blast_results[$threep_index][SEND]);
	    if (($pgg_blast_results[$fivep_index][QLEN] <= 500) && (($qid =~ /_BEG_/) || ($qid =~ /_END_/))) { # do not believe short sequences on the ends of contigs
		#this could change but for now these seem to be assembly/sequencing errors
	    } elsif ($repeat || ($revcomp5p != $revcomp3p) || ($revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] < $pgg_blast_results[$threep_index][SSTART])) || (!$revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$threep_index][SSTART]))) {
		#doesn't work for repeats - need to add check for circular contig causing third or fourth condition
	    } elsif ((($pgg_blast_results[$threep_index][QSTART] - $pgg_blast_results[$fivep_index][QEND]) >= -10) && (($pgg_blast_results[$threep_index][QSTART] - $pgg_blast_results[$fivep_index][QEND]) <= 10)) {
		my $fivep_start = $pgg_blast_results[$fivep_index][QSTART];
		my $fivep_end = $pgg_blast_results[$fivep_index][QEND];
		my $threep_start = $pgg_blast_results[$threep_index][QSTART];
		my $threep_end = $pgg_blast_results[$threep_index][QEND];
		my $del_start = $pgg_blast_results[$fivep_index][SEND];
		my $del_end = $pgg_blast_results[$threep_index][SSTART];
		my $revcomp = ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$threep_index][SEND]);
		if ($revcomp) {
		    $del_start--;
		    $del_end++;
		} else {
		    $del_start++;
		    $del_end--;
		}
		if ($fivep_end >= $threep_start) {
		    my $overlap = ($fivep_end - $threep_start) + 1;
		    if ($overlap % 2) {
			#odd
			$fivep_end -= ($overlap + 1) / 2;
			$threep_start += ($overlap - 1) / 2;
			if ($revcomp) {
			    $del_start += ($overlap + 1) / 2;
			    $del_end -= ($overlap - 1) / 2;
			} else {
			    $del_start -= ($overlap + 1) / 2;
			    $del_end += ($overlap - 1) / 2;
			}
		    } else {
			#even
			$fivep_end -= $overlap / 2;
			$threep_start += $overlap / 2;
			if ($revcomp) {
			    $del_start += $overlap / 2;
			    $del_end -= $overlap / 2;
			} else {
			    $del_start -= $overlap / 2;
			    $del_end += $overlap / 2;
			}
		    }
		}
		my $fivep_flank = substr($query_seqs{$qid}, ($fivep_start - 1), (($fivep_end - $fivep_start) + 1));
		my $threep_flank = substr($query_seqs{$qid}, ($threep_start - 1), (($threep_end - $threep_start) + 1));
		my $insertion = "";
		if ((($threep_start - $fivep_end) - 1) > 0) {
		    $insertion = substr($query_seqs{$qid}, $fivep_end, (($threep_start - $fivep_end) - 1));
		}
		my $deletion = "";
		my $tandem_duplication = 0;
		if ($revcomp) {
		    my $tmp_seq = "";
		    if ($del_start >= $del_end) {
			$tmp_seq = substr($pggdb_contigs{$cur_sid}, ($del_end - 1), (($del_start - $del_end) + 1));
		    } else {
			$tandem_duplication = 1;
			$tmp_seq = substr($pggdb_contigs{$cur_sid}, ($del_start - 1), (($del_end - $del_start) + 1));
		    }
		    $deletion = reverse($tmp_seq);
		    $deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		} else {
		    if ($del_start <= $del_end) {
			$deletion = substr($pggdb_contigs{$cur_sid}, ($del_start - 1), (($del_end - $del_start) + 1));
		    } else {
			$tandem_duplication = 1;
			$deletion = substr($pggdb_contigs{$cur_sid}, ($del_end - 1), (($del_start - $del_end) + 1));
		    }
		}
		if (($deletion !~ /NNNNN/) && (length($insertion) <= 200000) && (length($deletion) <= 200000)) { # do not believe sequences with gaps in them or are really long
		    my $del_len = length($deletion);
		    my $ins_len = length($insertion);
		    if ($tandem_duplication) {
			# swap insertion with deletion sequence for tandem duplication - should possibly be grabbing the sequence from the query rather than the matching bit from the subject
			$tandem_duplications++;
			print CALLS_FILE "$output\tTANDEM_DUPLICATION\t$qid\t$fivep_flank\t$threep_flank\t$insertion\t$ins_len\t$deletion\t$del_len\t$cur_sid\t$pid\t$del_start\t$del_end\t$fivep_end\t$threep_start\n";
		    } else {
			$deletions++;
			print CALLS_FILE "$output\tDELETION\t$qid\t$fivep_flank\t$threep_flank\t$deletion\t$del_len\t$insertion\t$ins_len\t$cur_sid\t$pid\t$del_start\t$del_end\t$fivep_end\t$threep_start\n";
		    }
		}
	    } elsif (($pgg_blast_results[$threep_index][QSTART] - $pgg_blast_results[$fivep_index][QEND]) > 10) {
		my $fivep_start = $pgg_blast_results[$fivep_index][QSTART];
		my $fivep_end = $pgg_blast_results[$fivep_index][QEND];
		my $threep_start = $pgg_blast_results[$threep_index][QSTART];
		my $threep_end = $pgg_blast_results[$threep_index][QEND];
		my $del_start = $pgg_blast_results[$fivep_index][SEND];
		my $del_end = $pgg_blast_results[$threep_index][SSTART];
		my $revcomp = ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$threep_index][SEND]);
		if ($revcomp) {
		    $del_start--;
		    $del_end++;
		} else {
		    $del_start++;
		    $del_end--;
		}
		my $fivep_flank = substr($query_seqs{$qid}, ($fivep_start - 1), (($fivep_end - $fivep_start) + 1));
		my $threep_flank = substr($query_seqs{$qid}, ($threep_start - 1), (($threep_end - $threep_start) + 1));
		my $insertion = substr($query_seqs{$qid}, $fivep_end, (($threep_start - $fivep_end) - 1));
		my $deletion = "";
		if ($revcomp) {
		    if ((($del_start - $del_end) + 1) > 0) {
			my $tmp_seq = substr($pggdb_contigs{$cur_sid}, ($del_end - 1), (($del_start - $del_end) + 1));
			$deletion = reverse($tmp_seq);
			$deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		} else {
		    if ((($del_end - $del_start) + 1) > 0) {
			$deletion = substr($pggdb_contigs{$cur_sid}, ($del_start - 1), (($del_end - $del_start) + 1));
		    }
		}
		if (($deletion !~ /NNNNN/) && (length($insertion) <= 200000) && (length($deletion) <= 200000)) { # do not believe sequences with gaps in them or are really long
		    my $del_len = length($deletion);
		    my $ins_len = length($insertion);
		    $insertions++;
		    print CALLS_FILE "$output\tINSERTION\t$qid\t$fivep_flank\t$threep_flank\t$deletion\t$del_len\t$insertion\t$ins_len\t$cur_sid\t$pid\t$del_start\t$del_end\t$fivep_end\t$threep_start\n";
		}
	    }
	} elsif ($fivep_index >= 0) {# possible deletion
	    if (($fivep_index > $#pgg_blast_results) || ($fivep_index < 0)) {
		die ("ERROR: fivep_index out of range $fivep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$fivep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	} elsif ($threep_index >= 0) {# possible deletion
	    if (($threep_index > $#pgg_blast_results) || ($threep_index < 0)) {
		die ("ERROR: threep_index out of range $threep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$threep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	}
    }
    close(PGG_BLAST_FILE);
    close(CALLS_FILE);
    open(OUT_FEATURES, ">", $out_features) || die ("Couldn't open $out_features for writing\n");
    print OUT_FEATURES "Mutations\tDeletions\tInsertions\tTandem Duplications\n$mutations\t$deletions\t$insertions\t$tandem_duplications\n";
    close (OUT_FEATURES);
    #$cmd = 'export LD_LIBRARY_PATH=/usr/local/packages/glibc-2.14/lib:$LD_LIBRARY_PATH; ' . "/usr/local/projdata/99999/IFX/CommonDB/ncbi-blast+/ncbi-blast-2.9.0+/bin/blastn -num_threads 4 -query $out_fasta_seqs -db $nrdb -out $out_nrdb_blast -task megablast -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 2";
    #`$cmd`;
    open(COMBINED_BTAB, ">", $combined_btab) || die ("Couldn't open $combined_btab for writing\n");
    print COMBINED_BTAB "qseqid\tsseqid\tpident\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tstitle\n";
    close (COMBINED_BTAB);
    #`cat $out_genome_blast $out_nrdb_blast $out_engdb_blast $out_PGG_blast | sort -k1,1 -k11,11nr >> $combined_btab`;
    #`rm $out_genome_blast $out_nrdb_blast $out_engdb_blast $out_PGG_blast`;
    `cat $out_genome_blast $out_PGG_blast | sort -k1,1 -k11,11nr >> $combined_btab`;
    `rm $out_genome_blast $out_PGG_blast`;
}

exit(0);
