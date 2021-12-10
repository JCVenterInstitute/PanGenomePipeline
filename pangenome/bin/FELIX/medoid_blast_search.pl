#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use Cwd;
use FileHandle;
use Getopt::Long;
use Carp;
use strict;

my $commandline = join (" ", @ARGV);
print STDERR "$commandline\n";
my $cwd = getcwd;
my $blast_task = "blastn";
my $help_text = 'This program BLASTs a FASTA file of medoids against a genome.
This work is done in a directory called PGGPDBS_XXXXXX, which is deleted when the program finishes.
The directory is created by mktemp and the XXXXXX is replaced by random characters.
The directory is created in the current working directory unless the use_local_disk flag is specified.
If using local disk the directory is created in the $HOME(~) directory on the grid node. This is risky if
there is not enough available space in this area but was implemented because of NFS issues for the
makebalstdb command on some systems. Should only use as a last resort.

Input Flags:
-medoids - The nucleotide multiFASTA file containing the medoids centroids.fasta for PanOCT (required)
-genome - A nucleotide multiFASTA file with of a target genome (required)
-blastout - The output file for the blast results (required)
-topology - topology file for the genome(s) being searched
-blast_directory - directory name for where blast executables are located - default is not to use a directory
-ld_load_directory - directory name for where blast libraries are located - default is not to use a directory
-blast_task - blast task to use, typically blastn or megablast - default is blastn
-strip_version - strip versions from contig names
-filter_anomalies - called for filter anomalies rather than medoid search
-combine_topology_ids - flag to use when the ids in the topology file need to be combined to create the id in contig file
-use_local_disk - flag to use local disk on grid nodes in the $HOME(~) directory instead of the current working directory
-help - Outputs this help text';

GetOptions('medoids=s' => \my $medoids,
	   'genome=s' => \my $genome,
	   'blast_directory=s' => \my $blast_directory,
	   'ld_load_directory=s' => \my $ld_load_directory,
	   'blast_task=s' => \ $blast_task,
	   'topology=s' => \my $topology_file,
	   'blastout=s' => \my $blastout,
	   'filter_anomalies' => \my $filter_anomalies,
	   'combine_topology_ids' => \ my $combine_topology_ids,
	   'use_local_disk' => \ my $use_local_disk,
	   'strip_version' => \my $strip_version,
	   'help' => \my $help);
	
if($help){
    print("$help_text\n");
    exit;
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
	print("$help_text\n");
	exit;
    } # if no value for option s (seq_file), quit with help menu
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
	$blast_directory = 'export LD_LIBRARY_PATH=' . $ld_load_directory . ':$LD_LIBRARY_PATH; ' . $blast_directory;
    } else {
	print STDERR "Error with -ld_load_directory $ld_load_directory\n";
	print("$help_text\n");
	exit;
    }
} else {
    $ld_load_directory = "";
}
	
if(!$genome or !$medoids or !$blastout or !$topology_file){
    die("Error: One or more of the required file arguments are missing: genome($genome) medoids($medoids) blastout($blastout) topology($topology_file)\n$help_text\n");
}

#Globals
my @blast_matches_raw = ();# contains all of the blast_matches input to the program via the -blast input file
my @matches_by_cluster = ();# contains all of the blast_matches input to the program via the -blast input file; sorted by cluster followed by bit score
my %contigs = ();  # key = contig ID, value = sequence
my %is_circular = (); # key = contig ID, value = 1 if contig is circular, 0 otherwise
my %headers = ();  # key = contig ID, value = fasta header line
my %ctglen = ();   # key = contig ID, value = contig length
my %added = ();    # key = contig ID, value = length added at beginning and end of circular contig

######################################################################################################################################################################
sub read_topology {

    unless (open (CIRCFILE, "<", "$topology_file") )  {
	die ("ERROR: can not open contig topology file $topology_file.\n");
    }
    while (<CIRCFILE>) {
	my $tag = "";
	my $asmbl_id = "";
	my $type = "";

	chomp;
	($tag, $asmbl_id, $type) = split(/\t/, $_);  # split the scalar $line on tab
	if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	    die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the contig topology file $topology_file.\nLine:\n$_\n");
	}
	if ($strip_version) {
	    $asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	if ($combine_topology_ids) {
	    $asmbl_id = $tag . "_" . $asmbl_id;
	}
	if (!defined $contigs{$asmbl_id}) {
	    die ("ERROR: $asmbl_id is a contig in the contig topology file $topology_file but not in the genome fasta file $genome!\nLine:\n$_\n");
	}
	if (defined $is_circular{$asmbl_id}) {
	    die ("ERROR: $asmbl_id occurs multiple times in the topology file $topology_file\n");
	}
	if ($type eq "circular") {
	    $is_circular{$asmbl_id} = 1;
	    $added{$asmbl_id} = $ctglen{$asmbl_id} > 20000 ? 20000 : $ctglen{$asmbl_id}; # add 20,000 bp or as much as possible to both ends of circular contigs
	} elsif ($type eq "linear") {
	    $is_circular{$asmbl_id} = 0;
	} else {
	    die ("ERROR: type $type must be either circular or linear in the  contig topology file $topology_file.\nLine:\n$_\n");
	}
    }
    close (CIRCFILE);
    return;
}

sub read_genome {  # read in the contigs for a genome and add 100,000 bp or as much as possible to both ends of circular contigs

    my $contigfile;
    unless (open ($contigfile, "<", $genome) )  {
	die ("cannot open file $genome!\n");
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
	$contigs{$id} = $sequence;
	if ($title =~ /^>/) {
	    $headers{$id} = $title . "\n";
	} else {
	    $headers{$id} = ">" . $title . "\n";
	}
	$ctglen{$id} = length ($sequence);
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($contigfile);
    return;
}

sub print_fasta { # print 60 character fasta lines

    my ($sequence, $length, $file) = @_;
    my $pos;
    for ( $pos = 0 ; $length > 60 ; $pos += 60 ) {
	print $file substr($sequence, $pos, 60), "\n";
	$length -= 60;
    }
    print $file substr($sequence, $pos, $length), "\n";
    return;
}

sub write_genome {  # read in the contigs for a genome and add a fixed amount of bp or as much as possible to both ends of circular contigs

    my $write_file = shift (@_);
    my $contigfile;
    unless (open ($contigfile, ">", $write_file) )  {
	die ("cannot open file $write_file!\n");
    }
    foreach my $contig (keys %contigs) {
	print $contigfile $headers{$contig};
	if ($is_circular{$contig}) {
	    &print_fasta((substr($contigs{$contig}, (-$added{$contig})) . $contigs{$contig} . substr($contigs{$contig}, 0, $added{$contig})), ($ctglen{$contig} + (2 * $added{$contig})), $contigfile);
	} else {
	    &print_fasta($contigs{$contig}, $ctglen{$contig}, $contigfile);
	}
    }
	close ($contigfile);
	return;
}

sub mod_blast { # eliminate blast matches to the added regions but keep one copy of overlapping matches
    my $blastin = shift (@_);
    my @btab_line = (); # array variable to store split btab lines
    my $qid = ""; # query id (cluster id)
    my $sid = ""; # subject id (contig from genome)
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $evalue; # blast evalue
    my $pid = ""; # percent identity
    my $score = ""; # BLAST bit score
    my $qlength = ""; # length of query sequence
    my $slength = ""; # length of subject sequence
    my $stitle = ""; # title of subject sequence
    my $line = ""; # raw input line
    my $blast_match_num = 0; #current blast match used for array index
    my $blast_in;
    my $blast_out_fp;
    
    unless (open($blast_in,"<", $blastin)) {
	die ("cannot open file $blastin!\n");
    }
    my $maxq_bit_thresh = 0;
    my $cur_clus = "";
    while ($line = <$blast_in>) {
	if ($line =~ /^#/)                                           # Skip header line
	{
	    next;
	}
	chomp $line;
	@btab_line = split(/\t/, $line);
	# ========================================================
	# btab output from NCBI blast+ blastn customized: -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\"
	# column number Description
	# 0      Query_id
	# 1	 subject_id (Hit from db)
	# 2	 % Identity
	# 3	 start of alignment on query (5' nucleotide match in query)
	# 4	 end of alignment on query (3' nucleotide match in query)
	# 5	 query length
	# 6	 start of alignment on subject (5' for query)
	# 7	 end of alignment on subject (3' for query)
	# 8	 subject length
	# 9      e-value
	# 10     score (bits)
	# 11     subject title
	# ========================================================
	$qid = $btab_line[0];
	$sid = $btab_line[1];
	if (!$filter_anomalies) {
	    $qid =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	    if ($strip_version) {
		$sid =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    }
	}
	$pid = $btab_line[2];
	$qbegin = $btab_line[3];
	$qend = $btab_line[4];
	$qlength = $btab_line[5];
	$sbegin = $btab_line[6];
	$send = $btab_line[7];
	$slength = $btab_line[8];
	$evalue = $btab_line[9];
	$score = $btab_line[10];
	$stitle = $btab_line[11];
	if ($cur_clus ne $qid) {
	    $maxq_bit_thresh = 0.5 * $score; #assumes blast output is sorted by qid then by bitscore per sid with highest bitscore sid first
	    $cur_clus = $qid;
	}
	if ($is_circular{$sid}) {
	    my $ctgbeg = $added{$sid} + 1;
	    my $ctgend = $added{$sid} + $ctglen{$sid};
	    if ((($sbegin < $ctgbeg) && ($send < $ctgbeg)) || (($sbegin > $ctgend) && ($send > $ctgend))) {
		next; # skip matches entirely in the added regions
	    }
	}
	if (($pid >= 90.0) && ($score > $maxq_bit_thresh)) {
	    $blast_matches_raw[$blast_match_num++] = {'clus' => $qid,      # cluster number
						      'ctg' => $sid,       # contig identifier
						      'pid' => $pid,       # percent identity
						      'qbeg' => $qbegin,   # start coordinate of cluster medoid
						      'qend' => $qend,     # end  coordinate of cluster medoid
						      'qlen' => $qlength,  # length of cluster medoid
						      'sbeg' => $sbegin,   # start coordinate on contig
						      'send' => $send,     # end  coordinate on contig
						      'ctglen' => $slength,# length of contig
						      'evalue' => $evalue, # blast evalue
						      'bits' => $score,    # bit score
						      'stitle' => $stitle, # subject title
						      'keep' => 1,         # set this to zero if not to be output
						      'line' => $line      # 0 if not best score for contig region, 1 otherwise
	    }
	}
    }
    close ($blast_in);
    my $sort_by_cluster_score = sub { # sort by cluster number then bit score
	my $cluster_test = $a->{'clus'} <=> $b->{'clus'};
		
	if ($cluster_test) {
	    return ($cluster_test);
	} else {
	    my $contig_test = $a->{'ctg'} cmp $b->{'ctg'};
	    
	    if ($contig_test) {
		return ($contig_test);
	    } else {
		my $bits_test = $b->{'bits'} <=> $a->{'bits'};
		if ($bits_test) {
		    return ($bits_test);
		} else {
		    return($a->{'sbeg'} <=> $b->{'sbeg'});
		}
	    }
	}
    };
    
    my $sort_by_name_score = sub { # sort by query name then bit score
	my $cluster_test = $a->{'clus'} cmp $b->{'clus'};
		
	if ($cluster_test) {
	    return ($cluster_test);
	} else {
	    my $contig_test = $a->{'ctg'} cmp $b->{'ctg'};
	    
	    if ($contig_test) {
		return ($contig_test);
	    } else {
		my $bits_test = $b->{'bits'} <=> $a->{'bits'};
		if ($bits_test) {
		    return ($bits_test);
		} else {
		    return($a->{'sbeg'} <=> $b->{'sbeg'});
		}
	    }
	}
    };
    
    if ($filter_anomalies) {
	@matches_by_cluster = sort $sort_by_name_score (@blast_matches_raw);
    } else {
	@matches_by_cluster = sort $sort_by_cluster_score (@blast_matches_raw);
    }

    unless (open($blast_out_fp,">", $blastout)) {
	die ("cannot open file $blastout!\n");
    }
    print $blast_out_fp "#qid\tsid\t%id\tqbeg\tqend\tqlen\tsbeg\tsend\tslen\tevalue\tbitscore\tstitle\n";
    foreach my $i (0 .. $#matches_by_cluster) {
	my $match = $matches_by_cluster[$i];
	my $cur_cluster = $match->{'clus'};
	my $cur_ctg = $match->{'ctg'};
	my $cur_bitscore = $match->{'bits'};
	if ($is_circular{$cur_ctg}) {
	    my $ctgadd = $added{$cur_ctg};
	    my $ctglen = $ctglen{$cur_ctg};
	    my $found = 0;
	    my $sbeg1 = $match->{'sbeg'};
	    my $send1 = $match->{'send'};
	    if (($sbeg1 <= $ctgadd) || ($send1 <= $ctgadd)) {
		foreach my $j (($i + 1) .. $#matches_by_cluster) {
		    my $new_match = $matches_by_cluster[$j];
		    if (($cur_cluster == $new_match->{'clus'}) && ($cur_ctg eq $new_match->{'ctg'}) && ($cur_bitscore == $new_match->{'bits'})) {
			my $sbeg2 = $new_match->{'sbeg'};
			my $send2 = $new_match->{'send'};
			if (($sbeg2 == ($sbeg1 + $ctglen)) && ($send2 == ($send1 + $ctglen))) {
			    $found = 1;
			    # keep only the match at the end of the circular contig which means a coordinate > contig length
			    $match->{'keep'} = 0;
			    # code to keep the larger part of the match between 1-contig length but this allows for both nonpositive coordinates and coordinates > contig length
#			    if ($sbeg1 <= $ctgadd) {
#				if (($send1 - $ctgadd) > (($send1 - $sbeg1) / 2)) {
#				    $new_match->{'keep'} = 0;
#				} else {
#				    $match->{'keep'} = 0;
#				}
#			    } else {
#				if (($sbeg1 - $ctgadd) > (($sbeg1 - $send1) / 2)) {
#				    $new_match->{'keep'} = 0;
#				} else {
#				    $match->{'keep'} = 0;
#				}
#			    }
			}
		    } else {
			last;
		    }
		}
		if (!$found) {
		    print STDERR "WARNING: equivalent match was not found at both ends of circular contig: $match->{'line'}\n";
		}
	    }
	    $match->{'sbeg'} = $sbeg1 - $added{$cur_ctg};
	    $match->{'send'} = $send1 - $added{$cur_ctg};
	    $match->{'ctglen'} = $ctglen{$cur_ctg};
	}
	if ($match->{'keep'}) {
	    print $blast_out_fp "$match->{'clus'}\t$match->{'ctg'}\t$match->{'pid'}\t$match->{'qbeg'}\t$match->{'qend'}\t$match->{'qlen'}\t$match->{'sbeg'}\t$match->{'send'}\t$match->{'ctglen'}\t$match->{'evalue'}\t$match->{'bits'}\t$match->{'stitle'}\n";
	}
    }
    close ($blast_out_fp);
    return;
}

{#main
    my $tmp_blast_dir;
    my $tmp_blast_out = $blastout . "_TMP";
    if ($use_local_disk) {
	$tmp_blast_dir = `mktemp -d -p ~ PGGPDBS_XXXXXX`;
    } else {
	$tmp_blast_dir = `mktemp -d PGGPDBS_XXXXXX`;
    }
    chomp $tmp_blast_dir;
    if (!(-d $tmp_blast_dir)) {
	die ("ERROR: could not create temporary blast directory PGGPDBS_XXXXXX($tmp_blast_dir) using mktemp\n")
    }
    &read_genome;
    &read_topology;
    &write_genome("$tmp_blast_dir/temp_fasta.ftmp");
    my $makeblastdb = $blast_directory . "makeblastdb";
    `$makeblastdb -in $tmp_blast_dir/temp_fasta.ftmp -dbtype nucl -out $tmp_blast_dir/temp_fasta.ftmp`;
    my $blastn = $blast_directory . "blastn";
    `$blastn -query $medoids -db $tmp_blast_dir/temp_fasta.ftmp -out $tmp_blast_out -task $blast_task -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore stitle\"`;
    &mod_blast("$tmp_blast_out");
    `rm -f -r $tmp_blast_dir $tmp_blast_out`;
}
