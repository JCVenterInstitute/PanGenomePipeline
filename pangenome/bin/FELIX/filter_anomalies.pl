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
my @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there are BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
my @ordered = ();      # Same as above but sorted by CONTIG then START
my @query_coords = (); # Same as above but first combining intervals and dropping most fields - sorted by QCCTG then QCBEG 
my @pgg_blast_results = ();  # These are the blast results of the query sequences against the PGG genomes
my %is_circular = (); # key = contig ID, value = 1 if contig is circular, 0 otherwise
my %contig_len = ();   # key = contig id, value = length of contig
my %contigs = ();      # key = contig id, value = sequence of contig
my %query_seqs = ();   # key = query id, value = sequence of query
my %pggdb_contig_len = ();   # key = contig id, value = length of contig
my %pggdb_contigs = ();      # key = contig id, value = sequence of contig
my %pggdb_is_circular = (); # key = contig ID, value = 1 if contig is circular, 0 otherwise
my %seen_contig = ();  # key = contig id, value = 1 (placeholder to show we've seen this contig)
my %max_bitscore = (); # key1 = query id, key2 = subject id, array = [MB5P,MB3P,MBFULL], value = max bitscore : pgg_blast_results index
my @categories = ("identical","conserved","gapped","divergent"); # literals used for output in ranges file

# CONSTANTS #
use constant MINALLOWPID => 98.0;
use constant CONTEXTLEN => 1000;
use constant MINCONTIGLEN => 1000;
use constant MAXEDGELEN => 100000;
use constant MAXINDELLEN => 200000;
use constant MINENDLEN => 500;
use constant QCCTG => 0;
use constant QCBEG => 1;
use constant QCEND => 2;
use constant QCNAME => 3;
use constant QCDEL => 4;
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
use constant CONSERVED => 1;
use constant GAPPED => 2;
use constant DIVERGED => 3;
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

my $bin_directory = "/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/";
my $genomes;
my $help;
my $debug;
my $blast_directory;
my $local_blast_directory;
my $ld_load_directory;
my $blast_task = "blastn";
my $nrdb;
my $PGGdb;
my $engdb;
my $pggdb_topology_file;
my $strip_version;
my $combine_topology_ids = 0;
my $use_existing_db = 0;
my $soft_mask_id = "";

GetOptions('genomes=s' => \ $genomes,
	   'strip_version' => \ $strip_version,
	   'combine_topology_ids' => \ $combine_topology_ids,
	   'use_existing_db' => \ $use_existing_db,
	   'help' => \ $help,
	   'debug' => \ $debug,
	   'bin_directory=s' => \ $bin_directory,
	   'blast_directory=s' => \ $blast_directory,
	   'ld_load_directory=s' => \ $ld_load_directory,
	   'PGGdb_topology=s' => \ $pggdb_topology_file,
	   'blast_task=s' => \ $blast_task,
	   'soft_mask_id=s' => \ $soft_mask_id,
	   'nrdb=s' => \ $nrdb,
	   'pggdb=s' => \ $PGGdb,
	   'engdb=s' => \ $engdb);

if (-d $bin_directory) {
    if (substr($bin_directory, 0, 1) ne "/") {
	$bin_directory = $cwd . "/$bin_directory";
    }
} else {
    print STDERR "The specified bin directory: $bin_directory does not exist or is not a directory!\n";
    $help = 1;
}
if ($blast_directory) {
    if (-d $blast_directory) {
	if (substr($blast_directory, -1, 1) ne "/") {
	    $blast_directory .= "/";
	}
	if (substr($blast_directory, 0, 1) ne "/") {
	    $blast_directory = $cwd . "/$blast_directory";
	}
    } else {
	print STDERR "Error with -blast_directory $blast_directory\n";
	$help = 1;
    }
} else {
    $blast_directory = "";
}

if ($ld_load_directory) {
    if (-d $ld_load_directory) {
	if (substr($ld_load_directory, -1, 1) ne "/") {
	    $ld_load_directory .= "/";
	}
	if (substr($ld_load_directory, 0, 1) ne "/") {
	    $ld_load_directory = $cwd . "/$ld_load_directory";
	}
	$local_blast_directory = 'export LD_LIBRARY_PATH=' . $ld_load_directory . ':$LD_LIBRARY_PATH; ' . $blast_directory;
    } else {
	print STDERR "Error with -ld_load_directory $ld_load_directory\n";
	$help = 1;
    }
} else {
    $ld_load_directory = "";
    $local_blast_directory = $blast_directory;
}

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('genomes=s' => genomes,
	   'help' => help,
	   'debug' => debug,
	   'combine_topology_ids' => \ combine_topology_ids,
	   'use_existing_db' => \ use_existing_db,
	   'strip_version' => \ strip_version,
	   'bin_directory=s' => \ bin_directory,
	   'blast_directory=s' => \ blast_directory,
	   'ld_load_directory=s' => \ ld_load_directory,
	   'PGGdb_topology=s' => \ pggdb_topology_file,
	   'blast_task=s' => \ blast_task,
	   'soft_mask_id=s' => \ soft_mask_id,
	   'nrdb=s' => nrdb,
	   'pggdb=s' => PGGdb,
	   'engdb=s' => engdb);
_EOB_
    exit(0);
}

######################################COMPONENT PROGRAM PATHS################################
my $medoid_blast_path = "$bin_directory/medoid_blast_search.pl";
#############################################################################################

# subroutine to handle substr for circular contigs
sub circ_substr { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file, $file_offset, $beg, $end , $ctg_len) = @_; # beg and end are in 1 base coordinates
    my $seq;
    my $len = ($end - $beg) + 1;
    my $local_save_input_separator = $/;
    seek($file, $file_offset, 0);
    $/="\n>";
    $seq = <$file>;
    $seq =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
    $/ = $local_save_input_separator; # restore the input separator
    if (length($seq) != $ctg_len) {
	die ("ERROR: Length of contig from seek position different than from initial reading.\n");
    }
    if ($end > $ctg_len) {
	if ($beg > $ctg_len) {
	    $beg -= $ctg_len;
	    return(substr($seq, ($beg - 1), $len));
	} else {
	    my $tmp_seq = substr($seq, ($beg - 1));
	    $tmp_seq .= substr($seq, 0, ($end - $ctg_len));
	    return($tmp_seq);
	}
    } else {
	return(substr($seq, ($beg - 1), $len));
    }
}

# subroutine to print contig segments out in fasta format
sub print_fasta { # have to adjust coordinates because they are in 1 base based coordinates and perl strings start at 0

    my ($file_handle, $seq_name, $seq, $beg, $end) = @_;
    #my $full_len = length($seq);
    #print STDERR "$seq_name $beg $end $full_len\n";
    $beg--;
    $end--;
    print $file_handle ">$seq_name\n";
    my $pos;
    my $seq_len = ($end - $beg) + 1;
    for ( $pos = $beg ; $seq_len > 60 ; $pos += 60 ) {
	print $file_handle substr($seq, $pos, 60), "\n";
	$seq_len -= 60;
    }
    print $file_handle substr($seq, $pos, $seq_len), "\n";
    return;
}

# read in contigs from pggdb genome file
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
    if ($strip_version) {
	$id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
    }
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
    $pggdb_contig_len{$id} = length($sequence);
    $title = ""; # clear the title for the next contig
    $sequence = ""; #clear out the sequence for the next contig
}
$/ = $pggdb_save_input_separator; # restore the input separator
close ($pggdbfile);

unless (open ($pggdbfile, "<", $PGGdb) )  {
    die ("cannot open PGG database file: $PGGdb!\n");
}
while ($pggdb_line = <$pggdbfile>) {
    if ($pggdb_line =~ /^>/) {
	my @fields = split(/\s+/, $pggdb_line);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$id =~ s/>\s*//; # remove leading > and spaces
	if ($strip_version) {
	    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	$pggdb_contigs{$id} = tell($pggdbfile); # save location of start of sequence
    }
}

#read in pggdb topolgy inforamtion for genome
unless (open (CIRCFILE, "<", "$pggdb_topology_file") )  {
    die ("ERROR: can not open pggdb contig topology file $pggdb_topology_file.\n");
}
while (<CIRCFILE>) {
    my $tag = "";
    my $asmbl_id = "";
    my $type = "";

    chomp;
    ($tag, $asmbl_id, $type) = split(/\t/, $_);  # split the scalar $line on tab
    if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the pggdb contig topology file $pggdb_topology_file.\nLine:\n$_\n");
    }
    if ($strip_version) {
	$asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
    }
    if ($combine_topology_ids) {
	$asmbl_id = $tag . "_" . $asmbl_id;
    }
    if (!defined $pggdb_contigs{$asmbl_id}) {
	die ("ERROR: $asmbl_id is a contig in the pggdb contig topology file $pggdb_topology_file but not in the pggdb fasta file $PGGdb!\nLine:\n$_\n");
    }
    if (defined $pggdb_is_circular{$asmbl_id}) {
	die ("ERROR: $asmbl_id occurs multiple times in the topology file $pggdb_topology_file\n");
    }
    if ($type eq "circular") {
	$pggdb_is_circular{$asmbl_id} = 1;
    } elsif ($type eq "linear") {
	$pggdb_is_circular{$asmbl_id} = 0;
    } else {
	die ("ERROR: type $type must be either circular or linear in the pggdb contig topology file $pggdb_topology_file.\nLine:\n$_\n");
    }
}
close (CIRCFILE);

# read file which specifies the output file prefix, assembly fasta file, assembly topology file, and anomalies file
open (my $infile, "<", $genomes) || die ("ERROR: cannot open file $genomes\n");
while (my $line = <$infile>)  {
    # clear data structures for the next genome
    @annotations = ();  # These are the lines of the attribute files but with 3 changes: 1) there BEST, VALUE, and TYPE fields 2) coordinates are now smallest then largest, not start then stop 3) there is an INVERT field to indicate strand rather than STOP being smaller than START
    @ordered = ();      # Same as above but sorted by CONTIG then START
    @query_coords = ();   # Combine intervals from ordered: QCCTG is contig, QCBEG is start coordinate, QCEND is stop coordinate, QCNAME is query sequence name, QCDEL for deleted - sort by QCCTG then QCBEG
    %is_circular = (); # key = contig ID, value = 1 if contig is circular, 0 otherwise
    %contigs = ();      # key = contig id, value = sequence of contig
    %contig_len = ();   # key = contig id, value = length of contig
    %seen_contig = ();  # key = contig id, value = 1 (placeholder to show we've seen this contig)
    @pgg_blast_results = ();  # These are the blast results of the query sequences against the PGG genomes
    %query_seqs = ();   # key1 = name, value = sequence
    %max_bitscore = (); # key1 = query id, key2 = subject id, array = [MB5P,MB3P,MBFULL], value = max bitscore : pgg_blast_results index
    
    chomp $line;
    (my $output, my $genome, my $target_topology_file, my $anomalies) = split(/\t/, $line);  # split the scalar $line on tab

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
	if ($strip_version) {
	    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	my $contig_length = length($sequence);
	#print STDERR "$id\t$contig_length";
	if (($contig_length >= MINCONTIGLEN) || ($id !~ /_cov_[01]$/)) { #only store contigs whose depth of coverage is > 1 or whose length is >= MINCONTIGLEN
	    $contigs{$id} = $sequence;
	    $contig_len{$id} = $contig_length;
	    #print STDERR "\t$id";
	} else {
	    $contigs{$id} = "IGNORE";
	    $contig_len{$id} = 0;
	}
	#print STDERR "\n";
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($contigfile);

    #read in topolgy inforamtion for genome
    unless (open (CIRCFILE, "<", "$target_topology_file") )  {
	die ("ERROR: can not open contig topology file $target_topology_file.\n");
    }
    while (<CIRCFILE>) {
	my $tag = "";
	my $asmbl_id = "";
	my $type = "";

	chomp;
	($tag, $asmbl_id, $type) = split(/\t/, $_);  # split the scalar $line on tab
	if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	    die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the contig topology file $target_topology_file.\nLine:\n$_\n");
	}
	if ($strip_version) {
	    $asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	if (!defined $contigs{$asmbl_id}) {
	    die ("ERROR: $asmbl_id is a contig in the contig topology file $target_topology_file but not in the genome fasta file $genome!\nLine:\n$_\n");
	}
	if (defined $is_circular{$asmbl_id}) {
	    die ("ERROR: $asmbl_id occurs multiple times in the topology file $target_topology_file\n");
	}
	if ($type eq "circular") {
	    $is_circular{$asmbl_id} = 1;
	} elsif ($type eq "linear") {
	    $is_circular{$asmbl_id} = 0;
	} else {
	    die ("ERROR: type $type must be either circular or linear in the contig topology file $target_topology_file.\nLine:\n$_\n");
	}
    }
    close (CIRCFILE);

    # read in anomalies file
    open(ANOM_FILE, "<", $anomalies) || die ("Couldn't open $anomalies\n");
    my $count = -1;
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
	if ($split_line[5] > MAXEDGELEN) {
	    print STDERR "WARNING: not expecting length of any anomaly to be greater than 100,000bp\n$line\n"; # warn on long length edge anomalies
	}
	my $category = $split_line[2];
	if (($category =~ /^divergent/) || ($category =~ /^very/) || ($category eq "uniq_edge")) {
	    $category = DIVERGED;
	} elsif ($category =~ /^identical_clus/) { #ignore edges as context anchors
	    $category = IDENTICAL;
	} elsif ($category =~ /^gapped/) {
	    $category = GAPPED;
	} elsif ($category =~ /^conserved_clus_allele$/) { #ignore edges as context anchors
	    $category = CONSERVED;
	} else {
	    next; # ignore other types
	}
	my $contig = $split_line[1];
	$contig =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	if (!defined $contigs{$contig}) {
	    die ("ERROR: $contig is a contig in the anomalies file $anomalies but not in the genome fasta file $genome!\n$line\n");
	}
	if ($contigs{$contig} eq "IGNORE") {
	    next;
	}
	if ($is_circular{$contig}) {
	    if ($split_line[3] < $split_line[4]) {
		if (((($contig_len{$contig} - $split_line[4]) + 1 + $split_line[3]) == $split_line[5]) && ((($split_line[4] - $split_line[3]) + 1) != $split_line[5])) {
		    if ($split_line[4] > $contig_len{$contig}) {
			$split_line[4] -= $contig_len{$contig}; # correct for going around the end
		    } else {
			$split_line[3] += $contig_len{$contig}; # correct for going around the end
		    }
		}
	    } else {
		if (((($contig_len{$contig} - $split_line[3]) + 1 + $split_line[4]) == $split_line[5]) && ((($split_line[3] - $split_line[4]) + 1) != $split_line[5])) {
		    if ($split_line[3] > $contig_len{$contig}) {
			$split_line[3] -= $contig_len{$contig}; # correct for going around the end
		    } else {
			$split_line[4] += $contig_len{$contig}; # correct for going around the end
		    }
		}
	    }
	}
	$count++;
	$annotations[$count][CONTIG] = $contig;     # contig
	$annotations[$count][CATEGORY] = $category;       # category
	$annotations[$count][TYPE] = $split_line[2];       # type
	$annotations[$count][LOCUS] = $split_line[6];      # locus_id
	$annotations[$count][GENOME] = $split_line[0];     # genome
	$annotations[$count][LENGTH] = $split_line[5];     # length
	$annotations[$count][DELETE] = 0;     # subsumed by another segment
	if ((!$is_circular{$contig}) && (($split_line[3] <= 0) || ($split_line[4] <= 0) || ($split_line[3] > $contig_len{$contig}) || ($split_line[4] > $contig_len{$contig}))) {
	    die ("ERROR: for a linear contig neither start ($split_line[3]) nor stop ($split_line[4]) coordinates should be <= 0 or > contig length ($contig_len{$contig})\n$line\n");
	}
	if ((($split_line[3] <= 0) || ($split_line[3] > $contig_len{$contig})) && (($split_line[4] <= 0) || ($split_line[4] > $contig_len{$contig}))) {
	    die ("ERROR: both start ($split_line[3]) and stop ($split_line[4]) coordinates should not be <= 0 or > contig length ($contig_len{$contig})\n$line\n");
	}
	if (($split_line[3] <= 0) || ($split_line[4] <= 0) || ($split_line[3] > $contig_len{$contig}) || ($split_line[4] > $contig_len{$contig})) {
	    if (($split_line[3] <= 0) || ($split_line[4] <= 0)) {
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
		    die ("ERROR: start ($split_line[3]) and stop ($split_line[4]) coordinates are not consistent with anomaly length ($split_line[5])\n$line\n"); # die on non-normalized annotations going around the end of the contig
		}
		$count++; # for circular contigs with anomalies across the begin/end coordinates duplicate match
		# change to be off of the end of a circular contig
		$split_line[3] += $contig_len{$contig};
		$split_line[4] += $contig_len{$contig};
		$annotations[$count][CONTIG] = $contig;     # contig
		$annotations[$count][CATEGORY] = $category;       # category
		$annotations[$count][TYPE] = $split_line[2];       # type
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
		    die ("ERROR: start ($split_line[3]) and stop ($split_line[4]) coordinates are not consistent with anomaly length ($split_line[5])\n$line\n"); # die on non-normalized annotations going around the end of the contig
		}
	    } else {
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
		    die ("ERROR: start ($split_line[3]) and stop ($split_line[4]) coordinates are not consistent with anomaly length ($split_line[5])\n$line\n"); # die on non-normalized annotations going around the end of the contig
		}
		$count++; # for circular contigs with anomalies across the begin/end coordinates duplicate match
		# change to be off of the beginning of a circular contig
		$split_line[3] -= $contig_len{$contig};
		$split_line[4] -= $contig_len{$contig};
		$annotations[$count][CONTIG] = $contig;     # contig
		$annotations[$count][CATEGORY] = $category;       # category
		$annotations[$count][TYPE] = $split_line[2];       # type
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
		    die ("ERROR: start ($split_line[3]) and stop ($split_line[4]) coordinates are not consistent with anomaly length ($split_line[5])\n$line\n"); # die on non-normalized annotations going around the end of the contig
		}
	    }
	} else {
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
		die ("ERROR: start ($split_line[3]) and stop ($split_line[4]) coordinates are not consistent with anomaly length ($split_line[5])\n$line\n"); # die on non-normalized annotations going around the end of the contig
	    }
	}
    }
    close(ANOM_FILE);
    
    if ($debug) {
	print STDERR "DEBUG***annotations\n";
	for (my $j=0; $j < @annotations; $j++) {
	    print STDERR ("$j: $annotations[$j][CONTIG]\t$annotations[$j][LOCUS]\t$annotations[$j][START]\t$annotations[$j][STOP]\t$annotations[$j][TYPE]\t$ordered[$j][CATEGORY]\t$annotations[$j][GENOME]\t$annotations[$j][INVERT]\n");
	}
    }

    # sort anomalies by contig then by start, store in ordered data-structure. 

    @ordered = sort { $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[CATEGORY] <=> $b->[CATEGORY] } @annotations; # sort on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***ordered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][CATEGORY]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }

    # eliminate overlaps between annotated segments by truncation - delete those that are contained by others - order of precedence: IDENTICAL, CONSERVED, GAPPED, DIVERGED
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
		    #want to save the larger side when this happens - possibly need to split the record and insert second piece back into sorted array
		    if (($ordered[$j][START] - $ordered[$i][START]) > ($ordered[$i][STOP] - $ordered[$j][STOP])) {
			$ordered[$i][STOP] = $ordered[$j][START] - 1;
		    } else {
			$ordered[$i][START] = $ordered[$j][STOP] + 1;
		    }
		    if ($ordered[$i][STOP] < $ordered[$i][START]) {
			$ordered[$i][DELETE] = 1;
			last;
		    }
		}
	    }
	}
    }

    # resort anomalies by deleted or not, then by contig, then by start, store in ordered data-structure. 

    @ordered = sort { $a->[DELETE] <=> $b->[DELETE] || $a->[CONTIG] cmp $b->[CONTIG] || $a->[START] <=> $b->[START] || $a->[CATEGORY] <=> $b->[CATEGORY] } @ordered; # sort on deleted or not, then on contig, then on start, then on type

    if ($debug) {
	print "DEBUG***filtered\n";
	for (my $j=0; $j < @ordered; $j++) {
	    print ("$j: $ordered[$j][CONTIG]\t$ordered[$j][LOCUS]\t$ordered[$j][START]\t$ordered[$j][STOP]\t$ordered[$j][TYPE]\t$ordered[$j][CATEGORY]\t$ordered[$j][GENOME]\t$ordered[$j][INVERT]\n");
	}
    }

    # open file for extracted interesting ranges
    my $out_ranges = $output . "_ranges.txt";
    my $file_ranges;
    unless (open ($file_ranges, ">", $out_ranges) )  {
	die ("cannot open file $out_ranges!\n");
    }


    # output segments which are diverged, or at the ends of contigs
    my $range_beg = 0;
    my $range_end = 0;
    my $cur_beg;
    my $cur_end;
    my $cur_len;
    my $cur_category;
    my $cur_contig;
    my $prev_beg = 1;
    my $prev_end = 0;
    my $prev_category = -1;
    my $prev_contig = "";
    my $last_output = 0;
    my $beg_conserved = 0;
    my $end_conserved = 0;
    my $beg_diverged = 0;
    my $end_diverged = 0;
    my $cur_type;
    my $cur_locus;
    my $diverged_type = "";
    my $qc_index = 0;
    my $first_qc_index = -1;
    for (my $i=0; $i < @ordered; $i++) {
	if ($ordered[$i][DELETE]) {
	    last;
	}
	$cur_contig = $ordered[$i][CONTIG];
	if ((!defined $contigs{$cur_contig}) || ($contigs{$cur_contig} eq "IGNORE")) {
	    #print STDERR "triaged: $cur_contig\n";
	    next; # skip triaged contigs
	}
	$cur_beg = $ordered[$i][START];
	$cur_end = $ordered[$i][STOP];
	$cur_len = ($cur_end - $cur_beg) + 1;
	$cur_type = $ordered[$i][TYPE];
	$cur_locus = $ordered[$i][LOCUS];
	$cur_category = $ordered[$i][CATEGORY];
	$seen_contig{$cur_contig} = 1;
	if ($prev_contig ne $cur_contig) { # first segment for this contig
	    if (($prev_contig ne "") && ($beg_diverged != 0)) {
		#print STDERR "New contig $cur_contig ($prev_contig) $beg_diverged:$end_diverged $diverged_type\n";
		$query_coords[$qc_index][QCNAME] = $prev_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type;
		$query_coords[$qc_index][QCCTG] = $prev_contig;
		$query_coords[$qc_index][QCBEG] = $beg_diverged;
		$query_coords[$qc_index][QCEND] = $end_diverged;
		$query_coords[$qc_index][QCDEL] = 0;
		if (($first_qc_index >= 0) && ($query_coords[$first_qc_index][QCBEG] <= 0) && ($end_diverged > $contig_len{$prev_contig})) {
		    $query_coords[$first_qc_index][QCDEL] = 1; # discard the off the beginning segment in favor of the off the end segment for ciruclar contigs
		    $first_qc_index = -1;
		    $query_coords[$qc_index][QCNAME] .= "_" . $query_coords[$first_qc_index][QCNAME];
		    $query_coords[$qc_index][QCEND] = $contig_len{$prev_contig} + $query_coords[$qc_index][QCEND];
		}
		$qc_index++;
	    }
	    $beg_diverged = 0;
	    $end_diverged = 0;
	    $beg_conserved = 0;
	    $end_conserved = 0;
	    $diverged_type = "";
	    #if (($prev_contig ne "") && (($contig_len{$prev_contig} - $prev_end) > 20)) { # include unannotated contig ends > 20 bp
		#my $tmp_seq = substr($contigs{$prev_contig}, $prev_end, ($contig_len{$prev_contig} - $prev_end));
		#if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		    #$query_coords[$qc_index][QCNAME] = $prev_contig . "_END_" . ($prev_end + 1) . "_" . $contig_len{$prev_contig};
		    #$query_coords[$qc_index][QCCTG] = $prev_contig;
		    #$query_coords[$qc_index][QCBEG] = $prev_end;
		    #$query_coords[$qc_index][QCEND] = $contig_len{$prev_contig};
		    #$query_coords[$qc_index][QCDEL] = 0;
		    #$qc_index++;
		#}
	    #}
	    #if ($cur_beg > 20) { # include unannotated contig ends > 20 bp
		#my $tmp_seq = substr($contigs{$cur_contig}, 0, ($cur_beg - 1));
		#if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		    #$query_coords[$qc_index][QCNAME] = $cur_contig . "_BEG_1_" . ($cur_beg - 1);
		    #$query_coords[$qc_index][QCCTG] = $cur_contig;
		    #$query_coords[$qc_index][QCBEG] = 1;
		    #$query_coords[$qc_index][QCEND] = $cur_beg - 1;
		    #$query_coords[$qc_index][QCDEL] = 0;
		    #$qc_index++;
		#}
	    #}
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
	    $prev_category = -1;
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
	
	if (($cur_category == IDENTICAL) || ($cur_category == CONSERVED)) {
	    if ($cur_len >= 100) {
		$beg_conserved = $cur_beg;
		$end_conserved = $cur_end;
		if ($beg_diverged != 0) {
		    if ($first_qc_index < 0) {
			$first_qc_index = $qc_index;
		    }
		    if ($end_diverged > $end_conserved) {
			$end_diverged = $end_conserved;
		    } elsif ($end_diverged < (($beg_conserved + CONTEXTLEN) - 1)) {
			$end_diverged = (($contig_len{$cur_contig} - $beg_conserved) > CONTEXTLEN) ? ($beg_conserved + CONTEXTLEN) : $contig_len{$cur_contig};
		    }
		    #print STDERR "End diverged region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
		    $query_coords[$qc_index][QCNAME] = $cur_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type;
		    $query_coords[$qc_index][QCCTG] = $cur_contig;
		    $query_coords[$qc_index][QCBEG] = $beg_diverged;
		    $query_coords[$qc_index][QCEND] = $end_diverged;
		    $query_coords[$qc_index][QCDEL] = 0;
		    $qc_index++;
		    $beg_diverged = 0;
		    $end_diverged = 0;
		    $diverged_type = "";
		}
	    } else {
		if ($beg_diverged == 0) {
		    # do not include starting short conserved as part of diverged but do use it as context
		    if ($beg_conserved == 0) { # only do this if there is not a longer previous conserved region
			$beg_conserved = $cur_beg;
			$end_conserved = $cur_end;
		    }  
		    #print STDERR "Short conserved region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
		} else {
		    $diverged_type .= "_" . $cur_type . "_" . $cur_locus . "_" . $cur_beg . "_" . $cur_end;
		}
		# do not stop the diverged region for a short conserved region
	    }
	} elsif ($cur_category == DIVERGED) {
	    #print STDERR "Diverged region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
	    if ($beg_diverged == 0) {
		if (($end_conserved - $beg_conserved) > CONTEXTLEN)  {
		    $beg_diverged = $end_conserved - CONTEXTLEN;
		} elsif ($beg_conserved > 0) {
		    $beg_diverged = $beg_conserved;
		} else {
		    $beg_diverged = 1; # instead of $cur_beg to include context around the diverged region
		}
		$diverged_type = $cur_type . "_" . $cur_locus . "_" . $cur_beg . "_" . $cur_end;
		#print STDERR "Begin new diverged region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
	    } else {
		$diverged_type .= "_" . $cur_type . "_" . $cur_locus . "_" . $cur_beg . "_" . $cur_end;
	    }
	    $end_diverged = (($contig_len{$cur_contig} - $cur_end) > CONTEXTLEN) ? ($cur_end + CONTEXTLEN) : $contig_len{$cur_contig}; #instead of #cur_end to include context around the diverged region
	    #print STDERR "Diverged region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
	} elsif ($cur_category == GAPPED) {
	    # do nothing for gapped - continue diverged or conserved as the case may be
	} else {
	    die ("ERROR: Unexpected type found: $cur_category\n");
	}
	$prev_beg = $cur_beg;
	$prev_end = $cur_end;
	$prev_category = $cur_category;
	$prev_contig = $cur_contig;
    }
    if (($prev_contig ne "") && ($range_beg != 0)) {
	print $file_ranges "$prev_contig\t$range_beg\t$range_end\t$categories[$prev_category]\n";
    }
    if (($prev_contig ne "") && (($contig_len{$prev_contig} - $range_end) > 0)) { # include unannotated contig ends > 0 bp
	print $file_ranges "$prev_contig\t", ($range_end + 1), "\t$contig_len{$prev_contig}\tunannotated\n";
    }
    close($file_ranges);
    if (($prev_contig ne "") && ($beg_diverged != 0)) {
	#print STDERR "Last diverged region $cur_contig $beg_diverged:$end_diverged $beg_conserved:$end_conserved $diverged_type\n";
	$query_coords[$qc_index][QCNAME] = $prev_contig . "_DIV_" . $beg_diverged . "_" . $end_diverged . "_" . $diverged_type;
	$query_coords[$qc_index][QCCTG] = $prev_contig;
	$query_coords[$qc_index][QCBEG] = $beg_diverged;
	$query_coords[$qc_index][QCEND] = $end_diverged;
	$query_coords[$qc_index][QCDEL] = 0;
	if (($first_qc_index >= 0) && ($query_coords[$first_qc_index][QCBEG] <= 0) && ($end_diverged > $contig_len{$prev_contig})) {
	    $query_coords[$first_qc_index][QCDEL] = 1; # discard the off the beginning segment in favor of the off the end segment for ciruclar contigs
	    $first_qc_index = -1;
	    $query_coords[$qc_index][QCNAME] .= "_" . $query_coords[$first_qc_index][QCNAME];
	    $query_coords[$qc_index][QCEND] = $contig_len{$prev_contig} + $query_coords[$qc_index][QCEND];
	}
	$qc_index++;
    }
    #if (($prev_contig ne "") && (($contig_len{$prev_contig} - $prev_end) > 20)) { # include unannotated contig ends > 20 bp
	#my $tmp_seq = substr($contigs{$prev_contig}, $prev_end, ($contig_len{$prev_contig} - $prev_end));
	#if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
	    #$query_coords[$qc_index][QCNAME] = $prev_contig . "_END_" . ($prev_end + 1) . "_" . $contig_len{$prev_contig};
	    #$query_coords[$qc_index][QCCTG] = $prev_contig;
	    #$query_coords[$qc_index][QCBEG] = $prev_end;
	    #$query_coords[$qc_index][QCEND] = $contig_len{$prev_contig};
	    #$query_coords[$qc_index][QCDEL] = 0;
	    #$qc_index++;
	#}
    #}

    # output entire contigs with no annotation
    foreach my $id (keys %contigs)  { # go through all contigs
	if ((!defined $seen_contig{$id}) && ($contigs{$id} ne "IGNORE")) {
	    #print STDERR "not seen: $id\t$contig_len{$id}\n";
	    my $tmp_seq = substr($contigs{$id}, 0, $contig_len{$id});
	    if ($tmp_seq !~ /NNNNN/) { # do not use sequences with gaps in them - perhaps should split on gaps instead
		$query_coords[$qc_index][QCNAME] = $id . "_WHOLE_1_" . $contig_len{$id};
		$query_coords[$qc_index][QCCTG] = $id;
		$query_coords[$qc_index][QCBEG] = 1;
		$query_coords[$qc_index][QCEND] = $contig_len{$id};
		$query_coords[$qc_index][QCDEL] = 0;
		$qc_index++;
	    }
	}
    }

    # open file for extracted interesting sequences
    my $out_fasta_seqs = $output . "_QUERY_SEQS.fasta";
    my $file_fasta_seqs;
    unless (open ($file_fasta_seqs, ">", $out_fasta_seqs) )  {
	die ("cannot open file $out_fasta_seqs!\n");
    }
    # resort combined anomalies by deleted or not, then by contig, then by start, store in ordered data-structure. 

    @query_coords = sort { $a->[QCDEL] <=> $b->[QCDEL] || $a->[QCCTG] cmp $b->[QCCTG] || $a->[QCBEG] <=> $b->[QCBEG] } @query_coords; # sort on deleted or not, then on contig, then on start

    if ($debug) {
	print "DEBUG***query_coords\n";
	for (my $j=0; $j < @query_coords; $j++) {
	    print ("$j: $query_coords[$j][QCCTG]\t$query_coords[$j][QCBEG]\t$query_coords[$j][QCEND]\t$query_coords[$j][QCNAME]\t$query_coords[$j][QCDEL]\n");
	}
    }

    # output segments which are diverged, or at the ends of contigs
    my $cur_name = "";
    for (my $i=0; $i < @query_coords; $i++) {
	if ($query_coords[$i][QCDEL]) {
	    last;
	}
	$cur_contig = $query_coords[$i][QCCTG];
	if (!defined $contigs{$cur_contig}) {
	    die ("ERROR: contig $cur_contig is not defined here but should be\n");
	}
	$cur_beg = $query_coords[$i][QCBEG];
	$cur_end = $query_coords[$i][QCEND];
	$cur_name = $query_coords[$i][QCNAME];
	if ($cur_beg > $contig_len{$cur_contig}) {
	    die ("ERROR: contig $cur_contig($contig_len{$cur_contig}) has segment $cur_name with bad coordinates $cur_beg:$cur_end\n");
	}
	if (($cur_beg <= 0) || ($cur_end <= 0)) {
	    die ("ERROR: contig $cur_contig has segment $cur_name with bad coordinates $cur_beg:$cur_end\n");
	}
	if ($cur_end > $contig_len{$cur_contig}) {
	    $query_seqs{$cur_name} = substr($contigs{$cur_contig}, $cur_beg);
	    $query_seqs{$cur_name} .= substr($contigs{$cur_contig}, 0, ($cur_end - $contig_len{$cur_contig}));
	} else {
	    $query_seqs{$cur_name} = substr($contigs{$cur_contig}, ($cur_beg - 1), (($cur_end - $cur_beg) + 1));
	}
	&print_fasta($file_fasta_seqs, $cur_name, $query_seqs{$cur_name}, 1, length($query_seqs{$cur_name}));
    }
    close($file_fasta_seqs);

    my $out_predictions = $output . "_CALLS";
    my $out_genome = $output . "_GENOME";
    my $out_genome_blast = $output . "_SELFBLAST.btab";
    my $out_PGG_blast = $output . "_PGG.btab";
    my $filtered_PGG_blast = $output . "_filtered_PGG.btab";
    my $out_nrdb_blast = $output . "_NRDB.btab";
    my $out_engdb_blast = $output . "_ENGDB.btab";
    my $combined_btab = $output . "_COMBINED.btab";
    my $out_features = $output . "_FEATURES";
    my $makeblastdb = $blast_directory . "makeblastdb";
    my $blastn = $local_blast_directory . "blastn";
    my $cmd;
    #`cp $genome $out_genome`;
    #$cmd = "$makeblastdb -in $out_genome -dbtype nucl";
    #`$cmd`;
    #$cmd = "$blastn -query $out_fasta_seqs -db $out_genome -out $out_genome_blast -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 3";
    #`$cmd`;
    #`rm $out_genome*`;
    #$cmd = "$blastn -query $out_fasta_seqs -db $engdb -out $out_engdb_blast -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 2";
    #`$cmd`;
    #$cmd = "$blastn -query $out_fasta_seqs -db $PGGdb -out $out_PGG_blast -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\"";
    if ($use_existing_db) {
	$cmd = "$blastn -query $out_fasta_seqs -db $PGGdb -out $out_PGG_blast -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 3";
	if ($soft_mask_id) {
	    $cmd .= " -db_soft_mask $soft_mask_id ";
	}	
	`$cmd`;
    } else {
	$cmd = "$medoid_blast_path -blastout $out_PGG_blast -genome $PGGdb -topology $pggdb_topology_file -blast_directory $blast_directory -ld_load_directory $ld_load_directory -medoids $out_fasta_seqs -blast_task $blast_task -filter_anomalies";
	if ($combine_topology_ids) {
	    $cmd .= " -combine_topology_ids";
	}
	`$cmd`;
    }
    # read in anomalies file
    open(PGG_BLAST_FILE, "<", $out_PGG_blast) || die ("Couldn't open $out_PGG_blast for reading\n");
    $count = 0;
    while (my $line = <PGG_BLAST_FILE>) {
	if ($line =~ /^#/) {
	    next; # skip comment lines
	}
	chomp($line);
	my @split_line = split(/\t/,$line);
	my $qid = $pgg_blast_results[$count][QSEQID] = $split_line[0];         # query id
	my $sid = $pgg_blast_results[$count][SSEQID] = $split_line[1];         # subject id
	if ($strip_version) {
	    $sid =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    $pgg_blast_results[$count][SSEQID] = $sid;
	}
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
	#print STDERR "$pid $qstart $qend $qlen\n";
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
    open(PGG_BLAST_FILE, ">", $filtered_PGG_blast) || die ("Couldn't open $filtered_PGG_blast for writing\n");
    open(CALLS_FILE, ">", $out_predictions) || die ("Couldn't open $out_predictions for writing\n");
    print CALLS_FILE "Sample\tType\tQuery ID\t5pflank\t3pflank\tDeleted\tDeletion Length\tInserted\tInsertion Length\tSubject ID\tPercent Identity\tSubject Start\tSubject End\tQuery Start\tQuery End\n";
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
	print STDERR "Findmax:$qid:$max_sid:$full_max:$full_index:$sum_max:$fivep_index:$threep_index:$#pgg_blast_results\n";
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
	my $start_offset;
	my $end_offset;
	if ($qid =~ /.*_DIV_([0-9]+)_([0-9]+)_/) {
	    $start_offset = $1 - 1;
	    $end_offset = $2;
	} else {
	    die ("ERROR qid $qid is not in expected format!\n");
	}
	if ($full_max > 0) {# possible mutation
	    if (($full_index > $#pgg_blast_results) || ($full_index < 0)) {
		die ("ERROR: full_index out of range $full_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$full_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    #print STDERR "$btab\n";
	    my $pid = $pgg_blast_results[$full_index][PIDENT];
	    if (($pid < MINALLOWPID) && ($qid !~ /_WHOLE_/) && ($qid !~ /_BEG_/) && ($qid !~ /_END_/)) {
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
		    my $tmp_seq = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $send, $sstart, $pggdb_contig_len{$cur_sid});
		    $deletion = reverse($tmp_seq);
		    $deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		} else {
		    $deletion = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $sstart, $send, $pggdb_contig_len{$cur_sid});
		}
		$mutations++;
		my $del_len = length($deletion);
		my $ins_len = length($insertion);
		$qstart += $start_offset;
		$qend += $start_offset;
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
	    #print STDERR "$btab\n";
	    $btab = join("\t", @{ $pgg_blast_results[$threep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    #print STDERR "$btab\n";
	    my $cur_sid = $pgg_blast_results[$fivep_index][SSEQID];
	    my $pid = ($pgg_blast_results[$fivep_index][PIDENT] + $pgg_blast_results[$threep_index][PIDENT]) / 2;
	    my $revcomp5p = ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$fivep_index][SEND]);
	    my $revcomp3p = ($pgg_blast_results[$threep_index][SSTART] > $pgg_blast_results[$threep_index][SEND]);
	    if (($pggdb_is_circular{$cur_sid}) && ($revcomp5p == $revcomp3p)) {
		if (!$revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$threep_index][SSTART]) && ($pgg_blast_results[$threep_index][SSTART] < ($pggdb_contig_len{$cur_sid} - $pgg_blast_results[$threep_index][SSTART]))) {
		    $pgg_blast_results[$threep_index][SSTART] += $pggdb_contig_len{$cur_sid};
		    $pgg_blast_results[$threep_index][SEND] += $pggdb_contig_len{$cur_sid};
		} elsif ($revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] < $pgg_blast_results[$threep_index][SSTART]) && ($pgg_blast_results[$fivep_index][SSTART] < ($pggdb_contig_len{$cur_sid} - $pgg_blast_results[$fivep_index][SSTART]))) {
		    $pgg_blast_results[$fivep_index][SSTART] += $pggdb_contig_len{$cur_sid};
		    $pgg_blast_results[$fivep_index][SEND] += $pggdb_contig_len{$cur_sid};
		}
	    }
	    if (($pgg_blast_results[$fivep_index][QLEN] <= MINENDLEN) && (($qid =~ /_BEG_/) || ($qid =~ /_END_/))) { # do not believe short sequences on the ends of contigs
		#this could change but for now these seem to be assembly/sequencing errors
	    } elsif ($repeat || ($revcomp5p != $revcomp3p) || ($revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] < $pgg_blast_results[$threep_index][SSTART])) || (!$revcomp5p && ($pgg_blast_results[$fivep_index][SSTART] > $pgg_blast_results[$threep_index][SSTART]))) {
		#doesn't work for repeats - need to add check for circular contig causing third or fourth condition
	    } elsif ((($pgg_blast_results[$threep_index][QSTART] - $pgg_blast_results[$fivep_index][QEND]) >= -50) && (($pgg_blast_results[$threep_index][QSTART] - $pgg_blast_results[$fivep_index][QEND]) <= 10)) {
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
			$tmp_seq = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_end, $del_start, $pggdb_contig_len{$cur_sid});
		    } else {
			$tandem_duplication = 1;
			$tmp_seq = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_start, $del_end, $pggdb_contig_len{$cur_sid});
		    }
		    $deletion = reverse($tmp_seq);
		    $deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		} else {
		    if ($del_start <= $del_end) {
			$deletion = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_start, $del_end, $pggdb_contig_len{$cur_sid});
		    } else {
			$tandem_duplication = 1;
			$deletion = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_end, $del_start, $pggdb_contig_len{$cur_sid});
		    }
		}
		if (($deletion !~ /NNNNN/) && (length($insertion) <= MAXINDELLEN) && (length($deletion) <= MAXINDELLEN)) { # do not believe sequences with gaps in them or are really long
		    my $del_len = length($deletion);
		    my $ins_len = length($insertion);
		    $fivep_start += $start_offset;
		    $fivep_end += $start_offset;
		    $threep_start += $start_offset;
		    $threep_end += $start_offset;
		    if ($tandem_duplication) {
			# swap insertion with deletion sequence for tandem duplication - should possibly be grabbing the sequence from the query rather than the matching bit from the subject
			$tandem_duplications++;
			$fivep_end -= $del_len;
			$threep_start += $del_len;
			print CALLS_FILE "$output\tTANDEM_DUPLICATION\t$qid\t$fivep_flank\t$threep_flank\t$insertion\t$ins_len\t$deletion\t$del_len\t$cur_sid\t$pid\t$del_end\t$del_start\t$fivep_end\t$threep_start\n";
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
		if ($revcomp) {
		    if ((($del_start - $del_end) + 1) < 0) {
			my $overlap = ($del_end - $del_start) - 1;
			if ($overlap % 2) {
			    #odd
			    $fivep_end -= ($overlap + 1) / 2;
			    $threep_start += ($overlap - 1) / 2;
			    $del_start += ($overlap + 1) / 2;
			    $del_end -= ($overlap - 1) / 2;
			} else {
			    #even
			    $fivep_end -= $overlap / 2;
			    $threep_start += $overlap / 2;
			    $del_start += $overlap / 2;
			    $del_end -= $overlap / 2;
			}
		    }
		} else {
		    if ((($del_end - $del_start) + 1) < 0) {
			my $overlap = ($del_start - $del_end) - 1;
			if ($overlap % 2) {
			    #odd
			    $fivep_end -= ($overlap + 1) / 2;
			    $threep_start += ($overlap - 1) / 2;
			    $del_start -= ($overlap + 1) / 2;
			    $del_end += ($overlap - 1) / 2;
			} else {
			    #even
			    $fivep_end -= $overlap / 2;
			    $threep_start += $overlap / 2;
			    $del_start -= $overlap / 2;
			    $del_end += $overlap / 2;
			}
		    }
		}
		my $fivep_flank = substr($query_seqs{$qid}, ($fivep_start - 1), (($fivep_end - $fivep_start) + 1));
		my $threep_flank = substr($query_seqs{$qid}, ($threep_start - 1), (($threep_end - $threep_start) + 1));
		my $insertion = substr($query_seqs{$qid}, $fivep_end, (($threep_start - $fivep_end) - 1));
		my $deletion = "";
		if ($revcomp) {
		    if ((($del_start - $del_end) + 1) > 0) {
			my $tmp_seq = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_end, $del_start, $pggdb_contig_len{$cur_sid});
			$deletion = reverse($tmp_seq);
			$deletion =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		} else {
		    if ((($del_end - $del_start) + 1) > 0) {
			$deletion = &circ_substr($pggdbfile, $pggdb_contigs{$cur_sid}, $del_start, $del_end, $pggdb_contig_len{$cur_sid});
		    }
		}
		if (($deletion !~ /NNNNN/) && (length($insertion) <= MAXINDELLEN) && (length($deletion) <= MAXINDELLEN)) { # do not believe sequences with gaps in them or are really long
		    my $del_len = length($deletion);
		    my $ins_len = length($insertion);
		    $insertions++;
		    $fivep_start += $start_offset;
		    $fivep_end += $start_offset;
		    $threep_start += $start_offset;
		    $threep_end += $start_offset;
		    print CALLS_FILE "$output\tINSERTION\t$qid\t$fivep_flank\t$threep_flank\t$deletion\t$del_len\t$insertion\t$ins_len\t$cur_sid\t$pid\t$del_start\t$del_end\t$fivep_end\t$threep_start\n";
		}
	    }
	} elsif ($fivep_index >= 0) {# possible deletion
	    if (($fivep_index > $#pgg_blast_results) || ($fivep_index < 0)) {
		die ("ERROR: fivep_index out of range $fivep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$fivep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    #print STDERR "$btab\n";
	} elsif ($threep_index >= 0) {# possible deletion
	    if (($threep_index > $#pgg_blast_results) || ($threep_index < 0)) {
		die ("ERROR: threep_index out of range $threep_index not in [0 - $#pgg_blast_results]\n");
	    }
	    my $btab = join("\t", @{ $pgg_blast_results[$threep_index] });
	    print PGG_BLAST_FILE "$btab\n";
	    #print STDERR "$btab\n";
	} else {
	    die ("ERROR: query sequence $qid no blast matches found against $PGGdb\n");
	}
    }
    close ($pggdbfile);
    close(PGG_BLAST_FILE);
    close(CALLS_FILE);
    open(OUT_FEATURES, ">", $out_features) || die ("Couldn't open $out_features for writing\n");
    print OUT_FEATURES "Mutations\tDeletions\tInsertions\tTandem Duplications\n$mutations\t$deletions\t$insertions\t$tandem_duplications\n";
    close (OUT_FEATURES);
    #$cmd = "$blastn -num_threads 4 -query $out_fasta_seqs -db $nrdb -out $out_nrdb_blast -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\" -culling_limit 2";
    #`$cmd`;
    open(COMBINED_BTAB, ">", $combined_btab) || die ("Couldn't open $combined_btab for writing\n");
    print COMBINED_BTAB "qseqid\tsseqid\tpident\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tstitle\n";
    close (COMBINED_BTAB);
    #`cat $out_genome_blast $out_nrdb_blast $out_engdb_blast $filtered_PGG_blast | sort -k1,1 -k11,11nr >> $combined_btab`;
    #`rm $out_genome_blast $out_nrdb_blast $out_engdb_blast $filtered_PGG_blast`;
    `cat $out_genome_blast $filtered_PGG_blast | sort -k1,1 -k11,11nr >> $combined_btab`;
    `rm $out_genome_blast $filtered_PGG_blast`;
}

exit(0);
