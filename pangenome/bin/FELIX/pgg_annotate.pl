#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

#PanGenome Annotation (Vector Approach) [ATTEMPT 1]

use FileHandle;
use Getopt::Long;
use strict;
my %pgg_core_edges = ();   # key = pgg core edge, value = 1 just a place holder
my @edges = ();            # the identifiers from the first column of the PGG file for each edge
my %hash_edges = ();       # store the edge identifiers for the target genome in a hash for easier lookup
my @cluster_scores = ();   # the best score for the cluster
my @cluster_bits = ();     # the best bitscore for the cluster to break ties when the score is the same
my @cluster_colindex = (); # the contig and column index of the best score for the cluster to break ties when the score and bitscore are the same
my @neighbors = ();        # array indexed by cluster number of array for NEAR_5, NEAR_3, CORE_5, CORE_3 which each have a hash of cluster_end where end is 5 or 3
my %is_circular = ();      # key = contig id, value = 1 if circular 0 otherwise
my @blast_matches_raw = ();# contains all of the blast_matches input to the program via the -blast input file
my @matches_by_cluster = ();# contains all of the blast_matches input to the program via the -blast input file; sorted by cluster followed by bit score
my @matches_by_region = ();# contains all of the blast_matches input to the program via the -blast input file; sorted by contig followed by bit scorecontig start coordinate
my @reduced_by_region = ();# contains only the blast_matches which are a best match for a cluster medoid or best match for a region; sorted by contig followed by bit scorecontig start coordinate
my @cluster_matches = ();  # array of hashes array index = cluster id, key = field name (best, last_best, last), values are array indices in @matches_by_cluster for the sorted best matches for the cluster's medoid
my %overlaps = ();         # key1 = match array index, key2 = match array index, value = max amount of overlap as a fraction 0-1
my %columns = ();          # hash of arrays of possible cluster assignments, key = contig id, value = array of columns, format hierarchical but only one level clus1;clus2;(clus3,clus4);clus5;(clus6,clus7,clus8);(clus5,clus9);clus10
my %columns_core = ();     # hash of arrays of hashes, key = contig id, value = array of hashes parallel to the columns hash of arrays, format 'is_core' = 1 if column contains a core 0 otherwise, each possible cluster in the column is a key in the hash = 1 if core 0 otherwise
my %columns_status = ();   # hash of array, key = contig id, value = array of hashes parallel to the columns hash of arrays, value = 1 if cluster ortholog, 2 if cluster paralog, 0 otherwise
my %nearest = ();          # key = contig id, parallesl %columns, value = array(number of columns) of array(four NEAR_5, NEAR_3, CORE_5, CORE_3)of hash key = cluster_inv, value = 1
my @cluster_size = ();     # array of the size of a cluster  which is the number of genomes the cluster members are present in
my @cluster_hits = ();     # array of the number of blast matches in the reduced_by_region array for a cluster
my %contigs = ();          # key = contig id, value = sequence of contig
my @medoids = ();          # index = cluster id, value = sequence of medoid
my @medoids_anno = ();     # index = cluster id, value = annotation from fasta header of medoid
my @core_clusters = ();    # index = cluster number, value = 1 if core 0 otherwise
my $num_clusters = 0;      # this is the number of clusters in the pan-genome graph determined from the cluster_weights file
my $num_genomes = 0;       # this is the number of genomes in the pan-genome graph determined from the cluster_weights file
my @match_table = ();      # this is a column for the matchtable for this genome
my @new_match_table = ();  # this is a column for the matchtable for this genome for new clusters not in the PGG
my %column_scores = ();    # this is a parallele data structure to %columns which records the scores of matches in each column
my @best_scores = ();      # array of best scores to calcuate median and quartile
my @second_scores = ();    # array of second best scores to calculte quartile
my $score_median = 6;      # this is the median score for best scores in a column - initialized variable is a place holder until computed
my $score_threshold = 4;   # this is the equal distance between the best score third quartile and second best first quartile - place holder until then
my $minimum_score = 3;     # this is the minmum score needed to keep a column around for further consideration
my $align_anchor_len = 25; # this is the amount of sequence from the Blast match to include as an anchor for the Needleman-Wunsch alignment
my $overlap_threshold = 0.75; # this is the fractional amount two matches need to overlap to be considerd overlapping in the reduced by region Blast processing
# indices for neighbors in %neighbors{$cluster}[$index]{$cluster_end} data structure and %nearest{cluster id}[column][index]{cluster_inv}
my $NEAR_5 = 0;
my $NEAR_3 = 1;
my $CORE_5 = 2;
my $CORE_3 = 3;

my $help_text = "THIS IS PLACEHOLDER TEXT\n";                             #<------------------------------------------------------------------------ TODO: write help text


GetOptions('neighbors=s' => \my $adjacency_file_path,                              # Location of CGN data output by "core_neighbor_finder.pl"
			'blast=s' => \my $blast_file_path,
			'help' => \my $help,
			'reannotate' => \my $reannotate,
			'core=s' => \my $core_list,
			'genome=s' => \my $genome,
	                'medoids=s' => \my $medoids_path,
	                'target=s' => \my $target,
	                'rootname=s' => \my $rootname,
	                'pgg=s' => \my $pgg,
			'clusters=s' => \my $cluster_weights_path);

			#NOTE: This program currently knows nothing about the genomes nor about the full set of edges. We think this is fine.
			#NOTE: We also are assuming 1) BLAST results are in a specific format 2) This format involves a custom midpoint field, and results sorted by coordinates           <------ IMPORTANT
			#NOTE: We also are assuming 3) CGN data is filtered in descending order of frequency
			
			
if (!$blast_file_path || !$cluster_weights_path || !$genome) 
{
    print ("This program requires gene neighborhood information, BLAST results, and list of core clusters as input\n$help_text\n");   #<----- update to include new mandatory input
    exit;
}
if ($help)
{
    print ("$help_text");
    exit;
}


######################################################################################################################################################################
sub read_neighbors                                  # reads in neighbors from core_neighbor_finder.pl output, this is how program gets CGN info
{
    unless (open(NEIGHBORS, "<", "$adjacency_file_path")) {
	die ("cannot open neighbors file: $adjacency_file_path!\n");
    }
    my $line = '';
    while ($line = <NEIGHBORS>)
	{
	    #print $line;
	    if ($line =~ /^#/)
	    {
		next;                                    # skip comment/header lines
	    }
	    if (!($line =~ /(\d+)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)/))
	    {
		print ("BAD CGN DATA:\n $line \n");                                # NOTE: SHOULD WE DIE HERE???????
		next;                                                              # Is it better to split on tab and sanitize later?
	    }
	    elsif ($line =~ /(\d+)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)\t([0123456789_, ]+|NONE)/)
	    {
		#print "$1:$2:$3:$4:$5\n";
		&store_neighbors($1,$2,$NEAR_3);                                   # Cluster, 3' nearest, index
		&store_neighbors($1,$3,$NEAR_5);                                   # Cluster, 5' nearest, index
		&store_neighbors($1,$4,$CORE_3);                                   # Cluster, 3' core, index
		&store_neighbors($1,$5,$CORE_5);                                   # Cluster, 5' core, index
	    }
	}
    close(NEIGHBORS);
    foreach my $cluster (1 .. $num_clusters) { # make the edges symmetric
	foreach my $choice (keys %{ $neighbors[$cluster][$NEAR_5] }) {
	    (my $near_5_clus, my $near_5_end) = split(/_/, $choice);
	    if ($cluster < $near_5_clus) {
		my $index = ($near_5_end == 5) ? $NEAR_5 : $NEAR_3;
		my $score = ($neighbors[$cluster][$NEAR_5]->{$choice} + $neighbors[$near_5_clus][$index]->{$cluster . "_5"}) / 2;
		#print "SYM $cluster:$NEAR_5:$choice:$neighbors[$cluster][$NEAR_5]->{$choice}:$neighbors[$near_5_clus][$index]->{$cluster . '_5'}:$score\n";
		$neighbors[$near_5_clus][$index]->{$cluster . "_5"} = $score;
		$neighbors[$cluster][$NEAR_5]->{$choice} = $score;
	    }
	}
	foreach my $choice (keys %{ $neighbors[$cluster][$NEAR_3] }) {
	    (my $near_3_clus, my $near_3_end) = split(/_/, $choice);
	    if ($cluster < $near_3_clus) {
		my $index = ($near_3_end == 5) ? $NEAR_5 : $NEAR_3;
		my $score = ($neighbors[$cluster][$NEAR_3]->{$choice} + $neighbors[$near_3_clus][$index]->{$cluster . "_3"}) / 2;
		#print "SYM $cluster:$NEAR_3:$choice:$neighbors[$cluster][$NEAR_3]->{$choice}:$neighbors[$near_3_clus][$index]->{$cluster . '_3'}:$score\n";
		$neighbors[$near_3_clus][$index]->{$cluster . "_3"} = $score;
		$neighbors[$cluster][$NEAR_3]->{$choice} = $score;
	    }
	}
	# Core neighbors cannot be made symmetric
	# Construct set of single copy core edges
	if ($core_clusters[$cluster]) {
	    foreach my $choice (keys %{ $neighbors[$cluster][$CORE_5] }) {
		$pgg_core_edges{"(" . $choice . "," . $cluster . "_5)"} = 1;
		$pgg_core_edges{"(" . $cluster . "_5" . "," . $choice . ")"} = 1;
	    }
	    foreach my $choice (keys %{ $neighbors[$cluster][$CORE_3] }) {
		$pgg_core_edges{"(" . $choice . "," . $cluster . "_3)"} = 1;
		$pgg_core_edges{"(" . $cluster . "_3" . "," . $choice . ")"} = 1;
	    }
	}
    }
    return;
}


################################################################################################################################################################
sub store_neighbors
{
	# PUTS NEIGHBORS IN DATA STRUCTURE SO THAT NEIGHBORS CAN BE IDENTIFIED LATER, TURN COUNT INTO RANK ORDER    
	# WHEN COMPUTING VECTOR, SOMETIMES DIVIDE BY 2 IF ORIENTATION WRONG  (NOT HANDLED HERE)
	# WHEN COMPUTING VECTOR, ALSO CHECK NEIGHBORS ON THE WRONG SIDE, ALSO DIVIDE BY 2 (NOT HANDLED HERE EITHER)
	# WHAT IS STORED IS SCORE EACH NEIGHBOR WOULD GIVE IF IN PROPER ORIENTATION
	
	
	my $cluster = shift;                                    # get arguments
	my $string_neighbor = shift;                             # get arguments
	my $index = shift;                              # get arguments
	my @split = ();
	if ($string_neighbor eq "NONE")
	{
		$neighbors[$cluster][$index]{"NONE"} = "NONE";                                                # When computing vectors, check for "NONE" as a key, before checking any counts
	}
	elsif (($index == $NEAR_5) || ($index == $NEAR_3))
	{
	    @split = split(/\,/,$string_neighbor);
	    my $total = 0;
		foreach my $element (@split)
		{
		    my @fields = split(/\s+/, $element);
		    my $count = $fields[1];                                                                  # Get count of neighbor
		    $total += $count;
		}
		foreach my $element (@split)
		{
		    my @fields = split(/\s+/, $element);
		    my $count = $fields[1];                                                                  # Get count of neighbor
		    my $adjacent = $fields[0];
		    $neighbors[$cluster][$index]{$adjacent} = $count / $total;            # weight by percentage of total neighbors
		    if (@fields == 3) {
			$neighbors[$cluster][$index]{$adjacent . "_d"} = $fields[2];
			#print "SCORE $cluster:$index:$adjacent:$neighbors[$cluster][$index]{$adjacent}:$neighbors[$cluster][$index]{$adjacent . '_d'}\n";
		    } else {
			#print "SCORE $cluster:$index:$adjacent:$neighbors[$cluster][$index]{$adjacent}\n";
		    }
		}
	}
	else
	{
		@split = split(/\,/,$string_neighbor);
		my $number_ranks = 0;
		my $last_count = 0;
		foreach my $element (@split)
		{
		    my @fields = split(/\s+/, $element);
		    my $count = $fields[1];                                                                   #Get count of neighbor
		    if ($count != $last_count)                                                       # If we haven't seen that count before...
		    {
			$number_ranks++;                                                             # Then we must increment the number of ranks
			$last_count = $count;                                                     #This assumes the entries are sorted by count
		    }
		}
		$last_count = 0;
		my $current_rank = $number_ranks;                                                        # Most frequent gets highest rank
		foreach my $element (@split)
		{
		    my @fields = split(/\s+/, $element);
		    my $count = $fields[1];                                                                  # Get count of neighbor
		    my $adjacent = $fields[0];
		    if (($count != $last_count) && ($last_count != 0))                               # If we hit a new count (and it's not the first rank)...
		    {
			$current_rank--;                                                             # Decrease current rank
		    }
		    $last_count = $count;
		    $neighbors[$cluster][$index]{$adjacent} = (($current_rank)/$number_ranks);            # Finally, assign score based on currnet rank/# ranks   (top rank gets n/n=100%, worst rank gets 1/n)
		    if (@fields == 3) {
			$neighbors[$cluster][$index]{$adjacent . "_d"} = $fields[2];
			#print "SCORE $cluster:$index:$adjacent:$neighbors[$cluster][$index]{$adjacent}:$neighbors[$cluster][$index]{$adjacent . '_d'}\n";
		    } else {
			#print "SCORE $cluster:$index:$adjacent:$neighbors[$cluster][$index]{$adjacent}\n";
		    }
		}
	}
	
	return;
	
}

#####################################################################################################
sub read_pgg                                               # Read in and store the set of edges in the PGG
{
    my $line;
    my @fields = ();
    unless (open(PGG, "<", $pgg)) {
	die ("cannot open PGG file: $pgg!\n");
    }
    while ($line = <PGG>) {
	chomp $line;
	my @edge_values = split(/\t/, $line);  # split the scalar $line on tab
	my $edge_id = shift @edge_values;
	if ($edge_id =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
	    push @edges, $edge_id;
	    #$edge_id = "edge".$1."_".$2."to".$3."_".$4;
	    #$cluster1 = $1;
	    #$cluster2 = $3;
	    #$whichend1 = $2;
	    #$whichend2 = $4;
	} else {
	    die ("ERROR: Bad edge formatting $edge_id in file $pgg.\n");
	}
    }
    close(PGG);
}

#####################################################################################################
sub read_core_list                                               # Read in the set of clusters which are single copy core to be used as anchors
{
    my $line;
    my @fields = ();
    foreach my $cluster (1 .. $num_clusters) {
	$core_clusters[$cluster] = 0; #initialize all clusters as not single copy core anchors
    }
    unless (open(CLUSTER_CORES, "<", $core_list)) {
	die ("cannot open cores file: $core_list!\n");
    }
    while ($line = <CLUSTER_CORES>) {
	chomp $line;
	$line =~ s/\s*//g; # remove all whitespace characters
	$line =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	$core_clusters[$line] = 1;
    }
    close(CLUSTER_CORES);
}

#####################################################################################################
sub read_cluster_sizes                                               # For each cluster store the size which is the number of genomes the cluster members are present in
{
    my $line;
    my @fields = ();
    unless (open(CLUSTER_SIZES, "<", $cluster_weights_path)) {
	die ("cannot open cluster sizes file: $cluster_weights_path!\n");
    }
    while ($line = <CLUSTER_SIZES>) {
	chomp $line;
	@fields = split(/\t/, $line);
	$cluster_size[$fields[0]] = $fields[1];
	if ($fields[1] > $num_genomes){
	    $num_genomes = $fields[1];
	}
	if ($fields[0] > $num_clusters) {
	    $num_clusters = $fields[0];
	}
    }
    close(CLUSTER_SIZES);
}

#####################################################################################################
sub read_blast                                                       # For each BLAST hit, store relevant data
{ # get tab-delimited BLAST results in custom format

    my @btab_line = (); # array variable to store split btab lines
    my $qid = ""; # query id (cluster id)
    my $sid = ""; # subject id (contig from genome)
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $inverted = 0; # 0 if match is on positive strand 1 otherwise
    my $pid = ""; # percent identity
    my $score = ""; # BLAST bit score
    my $qlength = ""; # length of query sequence
    my $slength = ""; # length of subject sequence
    my $line = ""; # raw input line
    my $blast_match_num = 0; #current blast match used for array index
    
    unless (open(BLAST_HITS, "<", "$blast_file_path")) {
	die ("cannot open blast file: $blast_file_path!\n");
    }
    while ($line = <BLAST_HITS>) {
	if ($line =~ /^#/)                                           # Skip header line
	{
	    next;
	}
	chomp $line;
	@btab_line = split(/\t/, $line);
	# ========================================================
	# btab output from NCBI blast+ blastn customized: -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore\"
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
	# ========================================================
	$qid = $btab_line[0];
	$qid =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	$sid = $btab_line[1];
	$sid =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	$pid = $btab_line[2];
	$qbegin = $btab_line[3];
	$qend = $btab_line[4];
	$qlength = $btab_line[5];
	$sbegin = $btab_line[6];
	$send = $btab_line[7];
	$slength = $btab_line[8];
	$score = $btab_line[10];
	if ($sbegin > $send) { #swap
	    my $temp_val = $send;
	    $send = $sbegin;
	    $sbegin = $temp_val;
	    $inverted = 1;
	} else {
	    $inverted = 0;
	}
	$blast_matches_raw[$blast_match_num++] = { 'clus' => $qid,      # cluster number
					       'ctg' => $sid,       # contig identifier
					       'pid' => $pid,       # percent identity
					       'qbeg' => $qbegin,   # start coordinate of cluster medoid
					       'qend' => $qend,     # end  coordinate of cluster medoid
					       'qlen' => $qlength,   # length of cluster medoid
					       'sbeg' => $sbegin,   # start coordinate on contig
					       'send' => $send,     # end  coordinate on contig
					       'sinv' => $inverted, # 0 if forward strand, 1 if reverse strand
					       'ctglen' => $slength,# length of contig
					       'bits' => $score,    # bit score
					       'weak' => 0,         # set to 1 if the match is deemed too weak
					       'keepclus' => 1,     # 0 if not best score for cluster, 1 otherwise
					       'keepctg' => 1       # 0 if not best score for contig region, 1 otherwise
	}
    }
    close (BLAST_HITS);
    return;
}


###################################################################################################
sub process_blast_by_cluster
# process the blast matches to ignore inferior blast matches for the same cluster
{
    my $sort_by_cluster_score = sub { # sort by cluster number then bit score
	my $cluster_test = $a->{'clus'} <=> $b->{'clus'};
		
	if ($cluster_test) {
	    return ($cluster_test);
	} else {
	    return ($b->{'bits'} <=> $a->{'bits'});
	}
    };
    
    @matches_by_cluster = sort $sort_by_cluster_score (@blast_matches_raw);

    my $cur_cluster = -1; # cluster numbers are assumed to start at 1 so set to impossible value
    my $best_bitscore = -1; # likewise for bitscore

    foreach my $i (0 .. $#matches_by_cluster) {
	my $match = $matches_by_cluster[$i];
	if ($cur_cluster != $match->{'clus'}) {
	    $cur_cluster = $match->{'clus'};
	    $best_bitscore = $match->{'bits'};
	    $match->{'keepclus'} = 1;
	    $cluster_matches[$cur_cluster]->{'best'} = $i;
	    $cluster_matches[$cur_cluster]->{'last_best'} = $i;
	} elsif ($match->{'bits'} > (0.9 * $best_bitscore)) {
	    $match->{'keepclus'} = 1;
	    $cluster_matches[$cur_cluster]->{'last_best'} = $i;
	} else { # ignore blast matches which are not within some threshold of the best match
	    $match->{'keepclus'} = 0;
	}
	$cluster_matches[$cur_cluster]->{'last'} = $i;
    }
}

###################################################################################################
sub calc_column
# recursively determine ordering of blast matches for a column
{
    my ($partcol, $first, $last, $core) = @_;
    my $concols = "";
    #print "CFL:$first:$last:$partcol\n";
    if ($first == $last) {
	# didn't add anything so add )
	#remove the trailing ","
	my $last_char = chop $partcol;
	if ($last_char ne ",") {
	    die "ERROR last character of $partcol was not a , but a $last_char\n";
	}
	$concols = $partcol . ")";
	#print "AA$concols\n";
	return ($concols);
    }
    my $first_hier = -1; #index of first match which doesn't overlap everything so need hierarchical matches
    foreach my $index (($first + 1) .. $last) {
	if (!defined $overlaps{$first}{$index}) {
	    if ($first_hier < 0) {
		$first_hier = $index;
	    }
	    if (($index == $first_hier) || (defined $overlaps{$first_hier}{$index})) { # too simplistic
		# only proceed if match overlaps the first hierarchical match otherwise accounted for downstream
		if ($core_clusters[$reduced_by_region[$index]->{'clus'}]) {
		    $core->{'is_core'} = 1;
		    $core->{$reduced_by_region[$index]->{'clus'} . "_" . $reduced_by_region[$index]->{'sinv'}} = 1;
		} else {
		    $core->{$reduced_by_region[$index]->{'clus'} . "_" . $reduced_by_region[$index]->{'sinv'}} = 0;
		}
		$concols .= (&calc_column (($partcol . $reduced_by_region[$index]->{'clus'} . "_" . $index . ","), $index, $last, $core)) . ";";
		#print "AB$concols\n";
	    }
	}
    }
    if ($concols eq "") {
	# didn't add anything so add )
	#remove the trailing ","
	my $last_char = chop $partcol;
	if ($last_char ne ",") {
	    die "ERROR last character of $partcol was not a , but a $last_char\n";
	}
	$concols = $partcol . ")";
    } else {
	#remove the trailing ";"
	my $last_char = chop $concols;
	if ($last_char ne ";") {
	    die "ERROR last character of $concols was not a ; but a $last_char\n";
	}
    }
    #print "AC$concols\n";
    return ($concols);
}

###################################################################################################
sub process_blast_by_region
# process the blast matches to ignore inferior matches for a region
{
    my $sort_by_contig_start = sub { # sort by cluster number then bit score
	my $contig_test = $a->{'ctg'} cmp $b->{'ctg'};
	
	if ($contig_test) {
	    return ($contig_test);
	} elsif ($a->{'sbeg'} <=> $b->{'sbeg'}) {
	    return ($a->{'sbeg'} <=> $b->{'sbeg'});
	} else {
	    return (($b->{'send'} - $b->{'sbeg'}) <=> ($a->{'send'} - $a->{'sbeg'}));
	}
    };
    
    @matches_by_region = sort $sort_by_contig_start (@blast_matches_raw);

    foreach my $i (0 .. $#matches_by_region) {
	my $ctg1 = $matches_by_region[$i]->{'ctg'};
	my $beg1 = $matches_by_region[$i]->{'sbeg'};
	my $end1 = $matches_by_region[$i]->{'send'};
	my $len1 = ($end1 - $beg1) + 1;
	foreach my $j (($i + 1) .. $#matches_by_region) {
	    my $ctg2 = $matches_by_region[$j]->{'ctg'};
	    my $beg2 = $matches_by_region[$j]->{'sbeg'};
	    my $end2 = $matches_by_region[$j]->{'send'};
	    my $len2 = ($end2 - $beg2) + 1;
	    my $overlap;
	    if (($ctg1 ne $ctg2) || ($beg2 >= $end1)){
		last;
	    }
	    if ($end2 < $end1) {
		$overlap = $len2;
	    } else {
		$overlap = ($end1 - $beg2) + 1;
	    }
	    my $cov1 = $overlap / $len1;
	    my $cov2 = $overlap / $len2;
	    my $maxcov = $cov1 > $cov2 ? $cov1 : $cov2;
	    if ($maxcov > ($overlap_threshold / 2)) {
		my $clus1 = $matches_by_region[$i]->{'clus'};
		my $clus2 = $matches_by_region[$j]->{'clus'};
		my $pid1 = $matches_by_region[$i]->{'pid'};
		my $pid2 = $matches_by_region[$j]->{'pid'};
		my $score1 = $matches_by_region[$i]->{'bits'};
		my $score2 = $matches_by_region[$j]->{'bits'};
		my $frac1 = ($matches_by_region[$i]->{'qend'} - $matches_by_region[$i]->{'qbeg'}) / $matches_by_region[$i]->{'qlen'};
		my $frac2 = ($matches_by_region[$j]->{'qend'} - $matches_by_region[$j]->{'qbeg'}) / $matches_by_region[$j]->{'qlen'};
		my $size1 = $cluster_size[$clus1];
		my $size2 = $cluster_size[$clus2];
		if ($clus1 == $clus2) {
		    if ($score1 < $score2) {
			$matches_by_region[$i]->{'keepctg'} = 0;
			$matches_by_region[$i]->{'weak'} = 1;
		    } else {
			$matches_by_region[$j]->{'keepctg'} = 0;
			$matches_by_region[$j]->{'weak'} = 1;
		    }
		}
		if ($maxcov > $overlap_threshold) {
		    if ($score1 <= 0.9 * $score2) {
			$matches_by_region[$i]->{'keepctg'} = 0;
			if (($size1 == 1) || ($size1 <= ($size2 / 10)) || (($frac1 < 0.5) && ($frac2 >= 0.9)) || ($pid1 < ($pid2 - 10)) || (($frac1 < 0.9) && ($pid1 < 90.0) && ($pid2 >= 95.0) && ($frac2 >= 0.9))) {
			    # ignore some matches from singleton or small size clusters and partial low percent identity matches
			    $matches_by_region[$i]->{'weak'} = 1;
			}
		    } elsif ($score2 <= 0.9 * $score1) {
			$matches_by_region[$j]->{'keepctg'} = 0;
			if (($size2 == 1) || ($size2 <= ($size1 / 10)) || (($frac2 < 0.5) && ($frac1 >= 0.9)) || ($pid2 < ($pid1 - 10)) || (($frac2 < 0.9) && ($pid2 < 90.0) && ($pid1 >= 95.0) && ($frac1 >= 0.9))) {
			    # ignore some matches from singleton or small size clusters and partial low percent identity matches
			    $matches_by_region[$j]->{'weak'} = 1;
			}
		    }
		}
	    }
	}
    }
    foreach my $cluster (1 .. $num_clusters) {
	$cluster_hits[$cluster] = 0; #initialize all clusters as having 0 blast matches in reduced set
	$cluster_scores[$cluster] = 0; #initialize all clusters as having 0 for best score
	$cluster_bits[$cluster] = 0; #initialize all clusters as having 0 for bitsscore
	$cluster_colindex[$cluster] = -1; #initialize all clusters as having colindex -1
    }
    my $index = 0;
    foreach my $match (@matches_by_region) {
	if ($match->{'keepctg'} || ($match->{'keepclus'} && !$match->{'weak'})) {
	    #print "PASS($index):$match->{'clus'}($cluster_size[$match->{'clus'}]):$match->{'ctg'}:$match->{'pid'}:$match->{'qbeg'}:$match->{'qend'}:$match->{'qlen'}:$match->{'sbeg'}:$match->{'send'}:$match->{'sinv'}:$match->{'ctglen'}:$match->{'bits'}:$match->{'keepclus'}:$match->{'keepctg'}:$match->{'weak'}\n";
	    push @reduced_by_region, $match;
	    $cluster_hits[$match->{'clus'}]++;
	    $index++;
	} else {
	    #print "FAIL:$match->{'clus'}($cluster_size[$match->{'clus'}]):$match->{'ctg'}:$match->{'pid'}:$match->{'qbeg'}:$match->{'qend'}:$match->{'qlen'}:$match->{'sbeg'}:$match->{'send'}:$match->{'sinv'}:$match->{'ctglen'}:$match->{'bits'}:$match->{'keepclus'}:$match->{'keepctg'}:$match->{'weak'}\n";
	}
    }
    my $first_overlap = 0;
    my $last_overlap = 0;
    my @num_overlaps = ();     # number of overlaps for a given match
    my @tmp_columns = ();      # current columns array to be added to the hash
    my $cur_ctg = "";
    my @tmp_cores = ();        # parallel array to tmp_columns to deisgnate if the current column contains a single copy core cluster
    my @tmp_scores = ();        # parallel array to tmp_columns to store scores for the matches in the column
    my @tmp_status = ();        # parallel array to tmp_columns to deisgnate if the current column contains a single copy core cluster
    foreach my $i (0 .. $#reduced_by_region) {
	my $ctg1 = $reduced_by_region[$i]->{'ctg'};
	my $beg1 = $reduced_by_region[$i]->{'sbeg'};
	my $end1 = $reduced_by_region[$i]->{'send'};
	my $score1 = $reduced_by_region[$i]->{'bits'};
	my $len1 = ($end1 - $beg1) + 1;
	if ($i > $last_overlap) {
	    my $column = "";
	    my $core = {'is_core' => 0}; # initialize to not containing a core
	    my $score = {'best_score' => 0}; #initialize to being empty
	    my $status = 0;
	    #print "FL:$first_overlap:$last_overlap\n";
	    if ($first_overlap == $last_overlap) {
		$column = $reduced_by_region[$first_overlap]->{'clus'} . "_" . $first_overlap;
		if ($core_clusters[$reduced_by_region[$first_overlap]->{'clus'}]) {
		    $core->{'is_core'} = 1;
		    $core->{$reduced_by_region[$first_overlap]->{'clus'} . "_" . $reduced_by_region[$first_overlap]->{'sinv'}} = 1;
		} else {
		    $core->{$reduced_by_region[$first_overlap]->{'clus'} . "_" . $reduced_by_region[$first_overlap]->{'sinv'}} = 0;
		}
	    } elsif ($first_overlap == ($last_overlap - 1)) {
		$column = $reduced_by_region[$first_overlap]->{'clus'} . "_" . $first_overlap . ";" . $reduced_by_region[$last_overlap]->{'clus'} . "_" . $last_overlap;
		if ($core_clusters[$reduced_by_region[$first_overlap]->{'clus'}]) {
		    $core->{'is_core'} = 1;
		    $core->{$reduced_by_region[$first_overlap]->{'clus'} . "_" . $reduced_by_region[$first_overlap]->{'sinv'}} = 1;
		} else {
		    $core->{$reduced_by_region[$first_overlap]->{'clus'} . "_" . $reduced_by_region[$first_overlap]->{'sinv'}} = 0;
		}
		if ($core_clusters[$reduced_by_region[$last_overlap]->{'clus'}]) {
		    $core->{'is_core'} = 1;
		    $core->{$reduced_by_region[$last_overlap]->{'clus'} . "_" . $reduced_by_region[$last_overlap]->{'sinv'}} = 1;
		} else {
		    $core->{$reduced_by_region[$last_overlap]->{'clus'} . "_" . $reduced_by_region[$last_overlap]->{'sinv'}} = 0;
		}
	    } else {
		my $first_hier = -1; #index of first match which doesn't overlap everything so need hierarchical matches
		foreach my $over_index ($first_overlap .. $last_overlap) {
		    $num_overlaps[$over_index] = scalar keys %{ $overlaps{$over_index} };
		    if ($num_overlaps[$over_index] == ($last_overlap - $first_overlap)) {
			# this match overlaps all overlaps in the set so output singly
			$column .= $reduced_by_region[$over_index]->{'clus'} . "_" . $over_index . ";";
			if ($core_clusters[$reduced_by_region[$over_index]->{'clus'}]) {
			    $core->{'is_core'} = 1;
			    $core->{$reduced_by_region[$over_index]->{'clus'} . "_" . $reduced_by_region[$over_index]->{'sinv'}} = 1;
			} else {
			    $core->{$reduced_by_region[$over_index]->{'clus'} . "_" . $reduced_by_region[$over_index]->{'sinv'}} = 0;
			}
			#print "BA$column\n";
		    } else {
			if ($first_hier < 0) {
			    $first_hier = $over_index;
			}
			if (($over_index == $first_hier) || (defined $overlaps{$first_hier}{$over_index})) {
			    # only proceed if match overlaps the first hierarchical match otherwise accounted for downstream
			    if ($core_clusters[$reduced_by_region[$over_index]->{'clus'}]) {
				$core->{'is_core'} = 1;
				$core->{$reduced_by_region[$over_index]->{'clus'} . "_" . $reduced_by_region[$over_index]->{'sinv'}} = 1;
			    } else {
				$core->{$reduced_by_region[$over_index]->{'clus'} . "_" . $reduced_by_region[$over_index]->{'sinv'}} = 0;
			    }
			    $column .= &calc_column (("(" . $reduced_by_region[$over_index]->{'clus'} . "_" . $over_index . ","), $over_index, $last_overlap, $core) . ";";
			    #print "BB$column\n";
			}
		    }
		}
		#remove the trailing ";"
		my $last_char = chop $column;
		if ($last_char ne ";") {
		    die "ERROR last character of $column was not a ; but a $last_char\n";
		}
	    }
	    #print "BC$column\n";
	    $first_overlap = $last_overlap = $i;
	    push @tmp_columns, $column;
	    push @tmp_cores, $core;
	    push @tmp_scores, $score;
	    push @tmp_status, $status;
	    if (($cur_ctg ne $ctg1) && ($cur_ctg ne "")) {
		$columns{$cur_ctg} = [ @tmp_columns ];
		@tmp_columns = ();
		$columns_core{$cur_ctg} = [ @tmp_cores ];
		@tmp_cores = ();
		$column_scores{$cur_ctg} = [ @tmp_scores ];
		@tmp_scores = ();
		$columns_status{$cur_ctg} = [ @tmp_status ];
		@tmp_status = ();
	    }
	}
	$cur_ctg = $ctg1;
	foreach my $j (($i + 1) .. $#reduced_by_region) {
	    my $ctg2 = $reduced_by_region[$j]->{'ctg'};
	    my $beg2 = $reduced_by_region[$j]->{'sbeg'};
	    my $end2 = $reduced_by_region[$j]->{'send'};
	    my $score2 = $reduced_by_region[$j]->{'bits'};
	    my $len2 = ($end2 - $beg2) + 1;
	    my $overlap;
	    if (($ctg1 ne $ctg2) || ($beg2 >= $end1)){
		last;
	    }
	    if ($end2 < $end1) {
		$overlap = $len2;
	    } else {
		$overlap = ($end1 - $beg2) + 1;
	    }
	    my $cov1 = $overlap / $len1;
	    my $cov2 = $overlap / $len2;
	    my $maxcov = $cov1 > $cov2 ? $cov1 : $cov2;
	    if ($maxcov > $overlap_threshold) {
		#print "OVER:$i:$j\n";
		$overlaps{$i}{$j} = $overlaps{$j}{$i} = $maxcov;
		if ($j > $last_overlap) {
		    $last_overlap = $j;
		}
	    }
	}
    }
    $columns{$cur_ctg} = [ @tmp_columns ];
    $columns_core{$cur_ctg} = [ @tmp_cores ];
    $columns_status{$cur_ctg} = [ @tmp_status ];
    $column_scores{$cur_ctg} = [ @tmp_scores ];
    return;
}

###################################################################################################
sub nw_align
# use dynamic programming to find best alignment of query to subject DNA sequences where full alignment of query is desired so all terminal query gaps are penalized
{
    my ($query, $subject) = @_;
    #print length($query), ":", length($subject), "\n";
    #print substr($query, 0, 80), "...", substr($query, -80), "\n";
    #print substr($subject, 0, 80), "...", substr($subject, -80), "\n";
    my @prev_score = (); #previous row in dynamic programming array
    my @cur_score = (); #current row in dynamic programming array
    my @prev_s_gap_len = (); #previous row in dynamic programming array
    my @cur_s_gap_len = (); #current row in dynamic programming array
    my @prev_q_gap_len = (); #previous row in dynamic programming array
    my @cur_q_gap_len = (); #current row in dynamic programming array
    my @prev_match_len = (); #previous row in dynamic programming array
    my @cur_match_len = (); #current row in dynamic programming array
    my @prev_matches = (); #previous row in dynamic programming array
    my @cur_matches = (); #current row in dynamic programming array
    my @prev_beg = (); #previous row in dynamic programming array
    my @cur_beg = (); #current row in dynamic programming array
    my $gap_open = -4; #gap open penalty
    my $gap_inc = -1; #gap increment penalty
    my $match = 3; #match reward
    my $mis = -3; #mismatch penalty
    my $beg = 0; #start coordibate of alignnment on query sequence
    my $end = 0; #stop coordinate of alignment on query sequence
    my $matches = 0; #number of non gap steps in the alignment
    my $best;

    foreach my $j (0 .. length($subject)) { #initialize @prev no terminal gap penalties for subject
	$cur_score[$j] = 0;      # score so far
	$cur_s_gap_len[$j] = 0;  # subject gap length (if 0 no gap)
	$cur_q_gap_len[$j] = 0;  # query gap length (if 0 no gap)
	$cur_match_len[$j] = 0;  # match (no gaps) length (if 0 previous position is a gap)
	$cur_matches[$j] = 0;    # number of matches (non gaps) in the alignment
	$cur_beg[$j] = $j + 1;   # start coordinate on query (need + 1 because real start is one further in current)
	$prev_score[$j] = 0;     # score so far
	$prev_s_gap_len[$j] = 0; # subject gap length (if 0 no gap)
	$prev_q_gap_len[$j] = 0; # query gap length (if 0 no gap)
	$prev_match_len[$j] = 0; # match (no gaps) length (if 0 previous position is a gap)
	$prev_matches[$j] = 0;   # number of matches (non gaps) in the alignment
	$prev_beg[$j] = $j + 1;  # start coordinate on query (need + 1 because real start is one further in current)
    }
    my $refcur_score = \@cur_score;
    my $refprev_score = \@prev_score;
    my $refcur_s_gap_len = \@cur_s_gap_len;
    my $refprev_s_gap_len = \@prev_s_gap_len;
    my $refcur_q_gap_len = \@cur_q_gap_len;
    my $refprev_q_gap_len = \@prev_q_gap_len;
    my $refcur_match_len = \@cur_match_len;
    my $refprev_match_len = \@prev_match_len;
    my $refcur_matches = \@cur_matches;
    my $refprev_matches = \@prev_matches;
    my $refcur_beg = \@cur_beg;
    my $refprev_beg = \@prev_beg;
    my $tmpref;
    foreach my $i (1 .. length($query)) {
	$tmpref = $refprev_score;
	$refprev_score = $refcur_score;
	$refcur_score = $tmpref;
	$tmpref = $refprev_s_gap_len;
	$refprev_s_gap_len = $refcur_s_gap_len;
	$refcur_s_gap_len = $tmpref;
	$tmpref = $refprev_q_gap_len;
	$refprev_q_gap_len = $refcur_q_gap_len;
	$refcur_q_gap_len = $tmpref;
	$tmpref = $refprev_match_len;
	$refprev_match_len = $refcur_match_len;
	$refcur_match_len = $tmpref;
	$tmpref = $refprev_matches;
	$refprev_matches = $refcur_matches;
	$refcur_matches = $tmpref;
	$tmpref = $refprev_beg;
	$refprev_beg = $refcur_beg;
	$refcur_beg = $tmpref;
	if ($i == 1) {
	    $refcur_score->[0] = $gap_open + $gap_inc;
	} else {
	    $refcur_score->[0] = $refprev_score->[0] + $gap_inc;
	}
	foreach my $j (1 .. length($subject)) {
	    my $gap_query = $refprev_score->[$j] + ($refprev_q_gap_len->[$j] ? $gap_inc : ($gap_open + $gap_inc));
	    my $gap_subject = $refcur_score->[$j - 1] + ($refcur_s_gap_len->[$j - 1] ? $gap_inc : ($gap_open + $gap_inc));
	    my $diag = $refprev_score->[$j - 1] + ((substr($query, ($i - 1), 1) eq substr($subject, ($j - 1), 1)) ? $match : $mis);
	    if (($diag >= $gap_query) && ($diag >= $gap_subject)) {
		$refcur_score->[$j] = $diag;
		$refcur_s_gap_len->[$j] = 0;
		$refcur_q_gap_len->[$j] = 0;
		$refcur_match_len->[$j] = $refprev_match_len->[$j - 1] + 1;
		$refcur_matches->[$j] = $refprev_matches->[$j - 1] + 1;
		$refcur_beg->[$j] = $refprev_beg->[$j - 1];
	    } elsif ($gap_subject >= $gap_query) {
		$refcur_score->[$j] = $gap_subject;
		$refcur_s_gap_len->[$j] = $refcur_s_gap_len->[$j - 1] + 1;
		$refcur_q_gap_len->[$j] = 0;
		$refcur_match_len->[$j] = 0;
		$refcur_matches->[$j] = $refcur_matches->[$j - 1];
		$refcur_beg->[$j] = $refcur_beg->[$j - 1];
	    } else {
		$refcur_score->[$j] = $gap_query;
		$refcur_s_gap_len->[$j] = 0;
		$refcur_q_gap_len->[$j] = $refprev_q_gap_len->[$j] + 1;
		$refcur_match_len->[$j] = 0;
		$refcur_matches->[$j] = $refprev_matches->[$j];
		$refcur_beg->[$j] = $refprev_beg->[$j];
	    }
		
	}
    }
    $best = $refcur_score->[0];
    #print "0:$refcur_score->[0];";
    $beg = 0;
    $end = 0;
    $matches = 0;
    foreach my $j (1 .. length($subject)) {
	#print "$j:$refcur_score->[$j];";
	if ($refcur_score->[$j] > $best) {
	    $best = $refcur_score->[$j];
	    $beg = $refcur_beg->[$j];
	    $end = $j;
	    $matches = $refcur_matches->[$j];
	}
    }
    #print "\n";
    return ($beg, $end, $matches);
}

###################################################################################################
sub nw_align_hash
# use dynamic programming to find best alignment of query to subject DNA sequences where full alignment of query is desired so all terminal query gaps are penalized
{
    my ($query, $subject) = @_;
    #print length($query), ":", length($subject), "\n";
    #print substr($query, 0, 80), "...", substr($query, -80), "\n";
    #print substr($subject, 0, 80), "...", substr($subject, -80), "\n";
    my @prev = (); #previous row in dynamic programming array
    my @cur = (); #current row in dynamic programming array
    my $gap_open = -4; #gap open penalty
    my $gap_inc = -1; #gap increment penalty
    my $match = 3; #match reward
    my $mis = -3; #mismatch penalty
    my $beg = 0; #start coordibate of alignnment on query sequence
    my $end = 0; #stop coordinate of alignment on query sequence
    my $matches = 0; #number of non gap steps in the alignment
    my $best;

    foreach my $j (0 .. length($subject)) { #initialize @prev no terminal gap penalties for subject
	$cur[$j] = { 'score' => 0,       # score so far
		      's_gap_len' => 0,  # subject gap length (if 0 no gap)
		      'q_gap_len' => 0,  # query gap length (if 0 no gap)
		      'match_len' => 0,  # match (no gaps) length (if 0 previous position is a gap)
		      'matches' => 0,    # number of matches (non gaps) in the alignment
		      'beg' => $j + 1    # start coordinate on query (need + 1 because real start is one further in current)
	};
	$prev[$j] = { 'score' => 0,      # score so far
		      's_gap_len' => 0,  # subject gap length (if 0 no gap)
		      'q_gap_len' => 0,  # query gap length (if 0 no gap)
		      'match_len' => 0,  # match (no gaps) length (if 0 previous position is a gap)
		      'matches' => 0,    # number of matches (non gaps) in the alignment
		      'beg' => 0         # start coordinate on query (need + 1 because real start is one further in current)
	};
    }
    my $refprev = \@prev;
    my $refcur = \@cur;
    my $tmpref;
    foreach my $i (1 .. length($query)) {
	$tmpref = $refprev;
	$refprev = $refcur;
	$refcur = $tmpref;
	if ($i == 1) {
	    $refcur->[0]{'score'} = $gap_open + $gap_inc;
	} else {
	    $refcur->[0]{'score'} = $refprev->[0]{'score'} + $gap_inc;
	}
	foreach my $j (1 .. length($subject)) {
	    my $gap_query = $refprev->[$j]{'score'} + ($refprev->[$j]{'q_gap_len'} ? $gap_inc : ($gap_open + $gap_inc));
	    my $gap_subject = $refcur->[$j - 1]{'score'} + ($refcur->[$j - 1]{'s_gap_len'} ? $gap_inc : ($gap_open + $gap_inc));
	    my $diag = $refprev->[$j - 1]{'score'} + ((substr($query, ($i - 1), 1) eq substr($subject, ($j - 1), 1)) ? $match : $mis);
	    if (($diag >= $gap_query) && ($diag >= $gap_subject)) {
		$refcur->[$j]{'score'} = $diag;
		$refcur->[$j]{'s_gap_len'} = 0;
		$refcur->[$j]{'q_gap_len'} = 0;
		$refcur->[$j]{'match_len'} = $refprev->[$j - 1]{'match_len'} + 1;
		$refcur->[$j]{'matches'} = $refprev->[$j - 1]{'matches'} + 1;
		$refcur->[$j]{'beg'} = $refprev->[$j - 1]{'beg'};
	    } elsif ($gap_subject >= $gap_query) {
		$refcur->[$j]{'score'} = $gap_subject;
		$refcur->[$j]{'s_gap_len'} = $refcur->[$j - 1]{'s_gap_len'} + 1;
		$refcur->[$j]{'q_gap_len'} = 0;
		$refcur->[$j]{'match_len'} = 0;
		$refcur->[$j]{'matches'} = $refcur->[$j - 1]{'matches'};
		$refcur->[$j]{'beg'} = $refcur->[$j - 1]{'beg'};
	    } else {
		$refcur->[$j]{'score'} = $gap_query;
		$refcur->[$j]{'s_gap_len'} = 0;
		$refcur->[$j]{'q_gap_len'} = $refprev->[$j]{'q_gap_len'} + 1;
		$refcur->[$j]{'match_len'} = 0;
		$refcur->[$j]{'matches'} = $refprev->[$j]{'matches'};
		$refcur->[$j]{'beg'} = $refprev->[$j]{'beg'};
	    }
		
	}
    }
    $best = $refcur->[0]{'score'};
    #print "0:$refcur->[0]{'score'};";
    $beg = 0;
    $end = 0;
    $matches = 0;
    foreach my $j (1 .. length($subject)) {
	#print "$j:$refcur->[$j]{'score'};";
	if ($refcur->[$j]{'score'} > $best) {
	    $best = $refcur->[$j]{'score'};
	    $beg = $refcur->[$j]{'beg'};
	    $end = $j;
	    $matches = $refcur->[$j]{'matches'};
	}
    }
    #print "\n";
    return ($beg, $end, $matches);
}

###################################################################################################
sub read_genome {  # read in the contigs for a genome and and store for later extraction

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
	$contigs{$id} = $sequence;
	if ($fields[1] eq "circular") {
	    $is_circular{$id} = 1;
	} else {
	    $is_circular{$id} = 0;
	}
	$title = ""; # clear the title for the next contig
	$sequence = ""; #clear out the sequence for the next contig
    }
    $/ = $save_input_separator; # restore the input separator
    close ($contigfile);
    return;
}

###################################################################################################
sub read_medoids {  # read in the medoids for the PGG

    my $medoidfile;
    unless (open ($medoidfile, "<", $medoids_path) )  {
	die ("cannot open medoid fasta file: $medoids_path!\n");
    }
    my $save_input_separator = $/;
    my $line;
    $/="\n>";
    while ($line = <$medoidfile>) {
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence (very cool)
	my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	$id =~ s/.*_//; # remove everything up to and including the last underscore
	$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	$medoids[$id] = $sequence;
	shift @fields;
	shift @fields;
	$medoids_anno[$id] = join(" ", @fields);
	$title = ""; # clear the title for the next medoid
	$sequence = ""; #clear out the sequence for the next medoid
    }
    $/ = $save_input_separator; # restore the input separator
    close ($medoidfile);
    return;
}

#####################################################################################################
sub determine_contig_nearest_neighbors                                                     # For each column, find nearest neighbors and core neighbors
{
#my %nearest = ();          # hash of arrays of hashes, key = contid id, value = array where index = column #, value = hash where key = {left|right}, value = nearest core column on that side
#my %columns_core = ();     # hash of arrays of hashes, key = contig id, value = array of hashes parallel to the columns hash of arrays, format 'is_core' = 1 if column contains a core 0 otherwise, each possible cluster in the column is a key in the hash = 1 if core 0 otherwise
    $, = ',';
    foreach my $contig (keys	%columns)                                     # Go through each contig
    {
	# traversal to the "right"              (left and right not 5' vs 3' because we aren't checking orientation until later)
	my $beg_window = 0;
	my $end_window = -1;
	my $beg_window_cores = 0;
	my $end_window_cores = -1;
	my @window_cells = (0,0,0);
	my @window_cells_cores = (0,0,0);
	my @tmp_neighbors = ();
	my @tmp_cores = ();
	my @tmp_cores_dists = ();
	my @index_cols = ();
	my @index_cores = ();
	my $last_index = $#{ $columns_core{$contig} };
	foreach my $count (0 .. $last_index)                           # Go through each column of each contig
	{
	    #print "Column $count: $columns{$contig}[$count]\n";
	    #print %{ $columns_core{$contig}[$count] };
	    #print "\n";
	    my $col_count = 0;
	    my $col_count_cores = 0;
	    my $core_found = 0;
	    $index_cols[$count] = $beg_window . ':' . $end_window;
	    $index_cores[$count] = $beg_window_cores . ':' . $end_window_cores;
	    foreach my $clus (keys %{ $columns_core{$contig}[$count] }) {
		if ($clus eq "is_core")
		{
		    next;
		}
		push @tmp_neighbors, $clus;
		$col_count++;
		if ($columns_core{$contig}[$count]{$clus}) {
		    push @tmp_cores, $clus;
		    push @tmp_cores_dists, $count;
		    $col_count_cores++;
		    $core_found = 1;
		}
	    }
	    $beg_window += shift @window_cells;
	    $end_window += $col_count;
	    push @window_cells, $col_count;
	    if ($core_found) {
		$beg_window_cores += shift @window_cells_cores;
		$end_window_cores += $col_count_cores;
		push @window_cells_cores, $col_count_cores;
	    }
	    #print "$beg_window:$end_window:$beg_window_cores:$end_window_cores\n";
	}
	#print "TMP_NEIGHBORS:";
	#print @tmp_neighbors;
	#print "\n";
	#print "TMP_CORES:";
	#print @tmp_cores;
	#print "\n";
	#print "is_circular: $is_circular{$contig}\n";
	if ($is_circular{$contig}) {
	    my $context = 0; # count when we have filled in the context
	    my $context_cores = 0; # count when we have filled in the context
	    foreach my $count (0 .. $last_index)                           # Go through each column of each contig
	    {
		#print "Column $count: $columns{$contig}[$count]\n";
		#print %{ $columns_core{$contig}[$count] };
		#print "\n";
		my $col_count = 0;
		my $col_count_cores = 0;
		my $core_found = 0;
		if ($context < 3) {
		    $index_cols[$count] = $beg_window . ':' . $end_window;
		}
		if ($context_cores < 3) {
		    $index_cores[$count] = $beg_window_cores . ':' . $end_window_cores;
		} else {
		    last;
		}
		foreach my $clus (keys %{ $columns_core{$contig}[$count] }) {
		    if ($clus eq "is_core")
		    {
			next;
		    }
		    push @tmp_neighbors, $clus;
		    $col_count++;
		    if ($columns_core{$contig}[$count]{$clus}) {
			push @tmp_cores, $clus;
			push @tmp_cores_dists, $count;
			$col_count_cores++;
			$core_found = 1;
		    }
		}
		$beg_window += shift @window_cells;
		$end_window += $col_count;
		push @window_cells, $col_count;
			$context++;
		if ($core_found) {
		    $beg_window_cores += shift @window_cells_cores;
		    $end_window_cores += $col_count_cores;
		    push @window_cells_cores, $col_count_cores;
		    $context_cores++;
		}
		#print "$beg_window:$end_window:$beg_window_cores:$end_window_cores\n";
	    }
	    #print "TMP_NEIGHBORS:";
	    #print @tmp_neighbors;
	    #print "\n";
	    #print "TMP_CORES:";
	    #print @tmp_cores;
	    #print "\n";
	}
	
	foreach my $count (0 .. $last_index)                           # Go through each column of each contig
	{
	    (my $beg, my $end) = split(':', $index_cols[$count]);
	    #print "IC$index_cols[$count]:$beg:$end\n";
	    $nearest{$contig}[$count][$NEAR_5] = {};
	    if ($end >= 0) {
		foreach my $index ($beg .. $end) {
		    (my $clus, my $inv) = split('_', $tmp_neighbors[$index]);
		    $nearest{$contig}[$count][$NEAR_5]->{$clus} = $inv;
		}
	    }
	    ($beg, $end) = split(':', $index_cores[$count]);
	    #print "IC$index_cores[$count]:$beg:$end\n";
	    $nearest{$contig}[$count][$CORE_5] = {};
	    if ($end >= 0) {
		foreach my $index ($beg .. $end) {
		    (my $clus, my $inv) = split('_', $tmp_cores[$index]);
		    #print "$index:$clus:$inv:$tmp_cores[$index]\n";
		    $nearest{$contig}[$count][$CORE_5]->{$clus} = $inv;
		    $nearest{$contig}[$count][$CORE_5]->{$clus . "_d"} = ($tmp_cores_dists[$index] < $count) ? ($count - $tmp_cores_dists[$index]) : (($count + $last_index + 1) - $tmp_cores_dists[$index]);
		    #print "LEFT:$contig:$count:$clus:$tmp_cores_dists[$index]:$last_index:$nearest{$contig}[$count][$CORE_5]->{$clus . '_d'}\n";
		}
	    }
	}
	
	# traversal to the "left"
	$beg_window = 0;
	$end_window = -1;
	$beg_window_cores = 0;
	$end_window_cores = -1;
	@window_cells = (0,0,0);
	@window_cells_cores = (0,0,0);
	@tmp_cores_dists = ();
	@tmp_neighbors = ();
	@tmp_cores = ();
	@index_cols = ();
	@index_cores = ();
	foreach my $count (reverse(0 .. $last_index))                           # Go through each column of each contig in reverse order
	{
	    #print "Column $count: $columns{$contig}[$count]\n";
	    #print %{ $columns_core{$contig}[$count] };
	    #print "\n";
	    my $col_count = 0;
	    my $col_count_cores = 0;
	    my $core_found = 0;
	    $index_cols[$count] = $beg_window . ':' . $end_window;
	    $index_cores[$count] = $beg_window_cores . ':' . $end_window_cores;
	    foreach my $clus (keys %{ $columns_core{$contig}[$count] }) {
		if ($clus eq "is_core")
		{
		    next;
		}
		push @tmp_neighbors, $clus;
		$col_count++;
		if ($columns_core{$contig}[$count]{$clus}) {
		    push @tmp_cores, $clus;
		    push @tmp_cores_dists, $count;
		    $col_count_cores++;
		    $core_found = 1;
		}
	    }
	    $beg_window += shift @window_cells;
	    $end_window += $col_count;
	    push @window_cells, $col_count;
	    if ($core_found) {
		$beg_window_cores += shift @window_cells_cores;
		$end_window_cores += $col_count_cores;
		push @window_cells_cores, $col_count_cores;
	    }
	    #print "$beg_window:$end_window:$beg_window_cores:$end_window_cores\n";
	}
	#print "TMP_NEIGHBORS:";
	#print @tmp_neighbors;
	#print "\n";
	#print "TMP_CORES:";
	#print @tmp_cores;
	#print "\n";
	#print "is_circular: $is_circular{$contig}\n";
	if ($is_circular{$contig}) {
	    my $context = 0; # count when we have filled in the context
	    my $context_cores = 0; # count when we have filled in the context
	    foreach my $count (reverse(0 .. $last_index))                           # Go through each column of each contig
	    {
		#print "Column $count: $columns{$contig}[$count]\n";
		#print %{ $columns_core{$contig}[$count] };
		#print "\n";
		my $col_count = 0;
		my $col_count_cores = 0;
		my $core_found = 0;
		if ($context < 3) {
		    $index_cols[$count] = $beg_window . ':' . $end_window;
		}
		if ($context_cores < 3) {
		    $index_cores[$count] = $beg_window_cores . ':' . $end_window_cores;
		} else {
		    last;
		}
		foreach my $clus (keys %{ $columns_core{$contig}[$count] }) {
		    if ($clus eq "is_core")
		    {
			next;
		    }
		    push @tmp_neighbors, $clus;
		    $col_count++;
		    if ($columns_core{$contig}[$count]{$clus}) {
			push @tmp_cores, $clus;
			push @tmp_cores_dists, $count;
			$col_count_cores++;
			$core_found = 1;
		    }
		}
		$beg_window += shift @window_cells;
		$end_window += $col_count;
		push @window_cells, $col_count;
		$context++;
		if ($core_found) {
		    $beg_window_cores += shift @window_cells_cores;
		    $end_window_cores += $col_count_cores;
		    push @window_cells_cores, $col_count_cores;
		    $context_cores++;
		}
		#print "$beg_window:$end_window:$beg_window_cores:$end_window_cores\n";
	    }
	    #print "TMP_NEIGHBORS:";
	    #print @tmp_neighbors;
	    #print "\n";
	    #print "TMP_CORES:";
	    #print @tmp_cores;
	    #print "\n";
	}
	foreach my $count (0 .. $last_index)                           # Go through each column of each contig
	{
	    #print "Column: $columns{$contig}[$count]\n";
	    #print %{ $columns_core{$contig}[$count] };
	    #print "\n";
	    (my $beg, my $end) = split(':', $index_cols[$count]);
	    #print "IC$index_cols[$count]:$beg:$end\n";
	    $nearest{$contig}[$count][$NEAR_3] = {};
	    if ($end >= 0) {
		foreach my $index ($beg .. $end) {
		    (my $clus, my $inv) = split('_', $tmp_neighbors[$index]);
		    $nearest{$contig}[$count][$NEAR_3]->{$clus} = $inv;
		}
	    }
	    ($beg, $end) = split(':', $index_cores[$count]);
	    #print "IC$index_cores[$count]:$beg:$end\n";
	    $nearest{$contig}[$count][$CORE_3] = {};
	    if ($end >= 0) {
		foreach my $index ($beg .. $end) {
		    (my $clus, my $inv) = split('_', $tmp_cores[$index]);
		    $nearest{$contig}[$count][$CORE_3]->{$clus} = $inv;
		    $nearest{$contig}[$count][$CORE_3]->{$clus . "_d"} = ($tmp_cores_dists[$index] > $count) ? ($tmp_cores_dists[$index] - $count) : (($tmp_cores_dists[$index] + $last_index + 1) - $count);
		    #print "RIGHT:$contig:$count:$clus:$tmp_cores_dists[$index]:$last_index:$nearest{$contig}[$count][$CORE_3]->{$clus . '_d'}\n";
		}
	    }
	    #print "NEAR_5: ";
	    #print %{ $nearest{$contig}[$count][$NEAR_5] };
	    #print "\nNEAR_3: ";
	    #print %{ $nearest{$contig}[$count][$NEAR_3] };
	    #print "\nCORE_5: ";
	    #print %{ $nearest{$contig}[$count][$CORE_5] };
	    #print "\nCORE_3: ";
	    #print %{ $nearest{$contig}[$count][$CORE_3] };
	    #print "\n";
	}
    }
    $, = '';
    return;
}

###################################################################################################
sub calc_median # calculate median and quartile scores
{
    @best_scores = sort {$a <=> $b} @best_scores;
    #print "BEST SCORES\n";
    #foreach my $score (@best_scores) {
	#print "$score\n";
    #}
    @second_scores = sort {$a <=> $b} @second_scores;
    #print "SECOND SCORES\n";
    #foreach my $score (@second_scores) {
	#print "$score\n";
    #}
    $score_median = $best_scores[(scalar @best_scores) / 2];
    $score_threshold = ($best_scores[(scalar @best_scores) / 4] + $second_scores[(scalar @second_scores) / 2]) / 2;
    print "Median $score_median Threshold $score_threshold\n";
    return;
}

###################################################################################################
sub reset_columns # reset the columns data structures flattening out hierarchical columns if possible
{
    foreach my $cluster (0 .. $#cluster_hits)
    {
	$cluster_hits[$cluster] = 0;
	$cluster_scores[$cluster] = 0; #initialize all clusters as having 0 for best score
	$cluster_bits[$cluster] = 0; #initialize all clusters as having 0 for bitsscore
	$cluster_colindex[$cluster] = -1; #initialize all clusters as having colindex -1
    }
    foreach my $contig (keys	%columns)                                     # Go through each contig
    {
	my @tmp_columns = ();
	my $status_count = 0;
	foreach my $column (@{ $columns{$contig} })                           # Go through each column of each contig
	{
	    my @choices = split(';', $column);
	    if ($columns_status{$contig}[$status_count] == -1) {
	    } elsif (((scalar @choices) == 1) && ($column =~ /^\(/)) {
		    my $last_char = chop $column;
		    if ($last_char ne ")") 
		    {
			die "ERROR $column ends in $last_char instead of )\n";
		    }
		    $column =~ s/^.//;                                         #remove first character
		    my @sec_choices = split(',', $column);
		    foreach my $sec_choice (@sec_choices) {
			push @tmp_columns, $sec_choice;
		    }
	    } else {
		push @tmp_columns, $column;
	    }
	    $status_count++;
	}
	$columns{$contig}  = [ @tmp_columns ];
	foreach my $count (0 .. $#{ $columns{$contig} })
	{
	    $columns_core{$contig}[$count] = {'is_core' => 0};     # hash of arrays of hashes, key = contig id, value = array of hashes parallel to the columns hash of arrays, format 'is_core' = 1 if column contains a core 0 otherwise, each possible cluster in the column is a key in the hash = 1 if core 0 otherwise
	    $columns_status{$contig}[$count] = 0;   # hash of array, key = contig id, value = array of hashes parallel to the columns hash of arrays, value = 1 if cluster ortholog, 2 if cluster paralog, 0 otherwise
	    $column_scores{$contig}[$count] = {'best_score' => 0};    # this is a parallele data structure to %columns which records the scores of matches in each column
	    my @choices = split(';', $columns{$contig}[$count]);
	    foreach my $choice (@choices)                                     # Go through each option of each column
	    {
		if ($choice =~ /^\(/)                                         # When there are (), each option needs to be further split
		{
		    my $last_char = chop $choice;
		    if ($last_char ne ")") 
		    {
			die "ERROR $choice ends in $last_char instead of )\n";
		    }
		    $choice =~ s/^.//;                                         #remove first character
		    my @sec_choices = split(',', $choice);
		    foreach my $sec_choice (@sec_choices) 
		    {
			(my $clus, my $index) = split('_', $sec_choice);
			$cluster_hits[$clus]++;
			if ($core_clusters[$clus]) {
			    $columns_core{$contig}[$count]->{'is_core'} = 1;
			    $columns_core{$contig}[$count]->{$clus . "_" . $reduced_by_region[$index]->{'sinv'}} = 1;
			} else {
			    $columns_core{$contig}[$count]->{$clus . "_" . $reduced_by_region[$index]->{'sinv'}} = 0;
			}
		    }
		} else {
		    (my $clus, my $index) = split('_', $choice);
		    $cluster_hits[$clus]++;
		    if ($core_clusters[$clus]) {
			$columns_core{$contig}[$count]->{'is_core'} = 1;
			$columns_core{$contig}[$count]->{$clus . "_" . $reduced_by_region[$index]->{'sinv'}} = 1;
		    } else {
			$columns_core{$contig}[$count]->{$clus . "_" . $reduced_by_region[$index]->{'sinv'}} = 0;
		    }
		}
	    }
	}
    }
    return;
}

###################################################################################################
sub score_from_neighbors
# For each contig, traverse both directions finding nearest neighbor/core neighbor [check all possibilities, retain best score]
# If there are () then nearest might not even be in a different column
# Check orientation so we can decide if 3' or 5' is where score goes - both traversals will fill both, but some will fill 3' on pass one, some on pass two
# This orientation information comes from the 'sinv' field - if it is 0, 3' is the larger index, if it is 1, 5' is the higher index
# Accessible via $reduced_by_region[$index]->{'sinv'} where $index is the portion after the underscore
# Save best choice within column, and, if good enough, also save second best
# Once we get down to the "foreach choice" level, find choice with best score and save it, and save second best if it meets some threshold
# We will look for neighbors on both sides - eventually, depending on which side finds the neighbor (and orientation) we can apply a penalty. For now, full score
# Store scores attached to BLAST index, since other scores are attached to index - this will let all scores be in one place - neighbor score can easily be attached to index, but index doesn't point back to column



{
    my $final = shift; #if this is the final run then choose only one match for a column
    my $final_hierarchical = shift; #if this is the final hierarchcial run then only allow a column to have one match if the best match is hierarchical and don't allow hierarchcial matches as secondary matches
    my $fraction = shift; # fractional threshold of which scores to keep in addition to the best score
    @best_scores = ();
    @second_scores = ();
    foreach my $contig (keys	%columns)                                     # Go through each contig
    {
	my $count = 0;
	foreach my $column (@{ $columns{$contig} })                           # Go through each column of each contig
	{
	    #print "Column $column: $columns{$contig}[$count]\n"; 
	    my @choices = split(';', $column);
	    my $choice_count = 0;
	    my $best = "";
	    my $best_score = 0;
	    my $best_bits = 0;
	    my $second_score = 0;
	    my $best_sum = 0;
	    my $best_clus = 0;
	    my $best_neighbor_score = 0;
	    my $best_index = -1;
#	    if (((scalar @choices) == 1) && (defined ($column_scores{$contig}[$count]->{$column}))) {
#		return;
#	    }
	    foreach my $choice (@choices)                                     # Go through each option of each column
	    {
		#print "$choice\n";
		if ($choice =~ /^\(/)                                         # When there are (), each option needs to be further split
		{
		    my $last_char = chop $choice;
		    if ($last_char ne ")") 
		    {
			die "ERROR $choice ends in $last_char instead of )\n";
		    }
		    $choice =~ s/^.//;                                         #remove first character
		    my @sec_choices = split(',', $choice);
		    my $sec_best = "";
		    my $sec_best_score = 0;                                  # store best score from within a set of ()
		    my $sec_best_bits = 0;                                  # store best score from within a set of ()
		    my $sec_sum = 0;
		    my $sec_clus = 0;
		    my $sec_neighbor_score = 0;
		    my $sec_index = -1;
		    my $nearhash = {};
		    foreach my $sec_choice (@sec_choices) 
		    {
			(my $clus, my $index) = split('_', $sec_choice);
			$nearhash->{$clus} = $reduced_by_region[$index]->{'sinv'};
		    }
		    foreach my $sec_choice (@sec_choices) 
		    {
			#print "$sec_choice\n";
			(my $clus, my $index) = split('_', $sec_choice);
			(my $left, my $left_neighbor, my $left_core, my $left_core_neighbor, my $left_core_dist) = &score_neighbors($reduced_by_region[$index]->{'sinv'},1,$contig,$count,$clus,$index,1,$nearhash);                       # get scores and identities of best neighbors on left
			(my $right, my $right_neighbor, my $right_core, my $right_core_neighbor, my $right_core_dist) = &score_neighbors($reduced_by_region[$index]->{'sinv'},0,$contig,$count,$clus,$index,1,$nearhash);                  # get scores and identities of best neighbors on right
			my $total = ($left + $left_core + $right + $right_core + $left_core_dist + $right_core_dist) + &homology_score($clus, $index);                                                                                      # get total neighbor score
			#print "BLAST $sec_choice:$reduced_by_region[$index]->{'clus'}($cluster_size[$reduced_by_region[$index]->{'clus'}]):$reduced_by_region[$index]->{'ctg'}:$reduced_by_region[$index]->{'pid'}:$reduced_by_region[$index]->{'qbeg'}:$reduced_by_region[$index]->{'qend'}:$reduced_by_region[$index]->{'qlen'}:$reduced_by_region[$index]->{'sbeg'}:$reduced_by_region[$index]->{'send'}:$reduced_by_region[$index]->{'sinv'}:$reduced_by_region[$index]->{'ctglen'}:$reduced_by_region[$index]->{'bits'}:$reduced_by_region[$index]->{'keepclus'}:$reduced_by_region[$index]->{'keepctg'}:$reduced_by_region[$index]->{'weak'}\n";
			#print "$total:$left:$left_neighbor:$right:$right_neighbor:$left_core:$left_core_neighbor:$right_core:$right_core_neighbor:$left_core_dist:$right_core_dist\n";
			if ($total > $sec_best_score) { # perhaps keep more than just the best if the scores are close
			    $sec_best_score = $total;
			    $sec_best_bits = $reduced_by_region[$index]->{'bits'};
			    $sec_best = $sec_choice;
			    $sec_clus = $clus;
			    $sec_neighbor_score = $left + $left_core + $right + $right_core + $left_core_dist + $right_core_dist;
			    $sec_index = $index;
			}
			$sec_sum += $total - $score_median;
		    }
		    $column_scores{$contig}[$count]{$choice} = $sec_best_score;
		    if (($sec_best_score > $best_score) || (($sec_best_score == $best_score) && ($sec_sum > $best_sum))) {
			$second_score = $best_score;
			$best_score = $sec_best_score;
			$best_bits = $sec_best_bits;
			$best = $choice;
			$best_clus = $sec_clus;
			$best_sum = $sec_sum;
			$best_neighbor_score = $sec_neighbor_score;
			$best_index = $sec_index;
		    }
		} else {
		    (my $clus, my $index) = split('_', $choice);
		    (my $left, my $left_neighbor, my $left_core, my $left_core_neighbor, my $left_core_dist) = &score_neighbors($reduced_by_region[$index]->{'sinv'},1,$contig,$count,$clus,$index,0,0);                       # get scores and identities of best neighbors on left
		    (my $right, my $right_neighbor, my $right_core, my $right_core_neighbor, my $right_core_dist) = &score_neighbors($reduced_by_region[$index]->{'sinv'},0,$contig,$count,$clus,$index,0,0);                  # get scores and identities of best neighbors on right
		    my $total = ($left + $left_core + $right + $right_core + $left_core_dist + $right_core_dist) + &homology_score($clus, $index);                                                                                      # get total neighbor score
		    #print "BLAST $choice:$reduced_by_region[$index]->{'clus'}($cluster_size[$reduced_by_region[$index]->{'clus'}]):$reduced_by_region[$index]->{'ctg'}:$reduced_by_region[$index]->{'pid'}:$reduced_by_region[$index]->{'qbeg'}:$reduced_by_region[$index]->{'qend'}:$reduced_by_region[$index]->{'qlen'}:$reduced_by_region[$index]->{'sbeg'}:$reduced_by_region[$index]->{'send'}:$reduced_by_region[$index]->{'sinv'}:$reduced_by_region[$index]->{'ctglen'}:$reduced_by_region[$index]->{'bits'}:$reduced_by_region[$index]->{'keepclus'}:$reduced_by_region[$index]->{'keepctg'}:$reduced_by_region[$index]->{'weak'}\n";
		    #print "$total:$left:$left_neighbor:$right:$right_neighbor:$left_core:$left_core_neighbor:$right_core:$right_core_neighbor:$left_core_dist:$right_core_dist\n";
		    $column_scores{$contig}[$count]{$choice} = $total;
		    if ($total > $best_score) { # perhaps keep more than just the best if the scores are close
			$second_score = $best_score;
			$best_score = $total;
			$best_bits = $reduced_by_region[$index]->{'bits'};
			$best = $choice;
			$best_clus = $clus;
			$best_sum = $total - $score_median;
			$best_neighbor_score = $left + $left_core + $right + $right_core + $left_core_dist + $right_core_dist;
			$best_index = $index;
		    }
		}
		$choice_count++;                                                  # BE SURE TO DO THIS AFTER SCORE COMPUTATION AND STORING
	    }
	    if (!$core_clusters[$best_clus]) {
		push @best_scores, $best_score;
		if ($second_score) {
		    push @second_scores, $second_score;
		}
	    }
	    $column_scores{$contig}[$count]{'best_score'} = $best_score;
	    if ($best_score > $cluster_scores[$best_clus]) {
		$cluster_scores[$best_clus] = $best_score;
		$cluster_bits[$best_clus] = $best_bits;
		$cluster_colindex[$best_clus] = $contig . "_" . $count;
	    } elsif (($best_score == $cluster_scores[$best_clus]) && ($best_bits > $cluster_bits[$best_clus])) {
		$cluster_scores[$best_clus] = $best_score;
		$cluster_bits[$best_clus] = $best_bits;
		$cluster_colindex[$best_clus] = $contig . "_" . $count;
	    }
	    $columns{$contig}[$count] = $best;
	    #print "COLUMN $contig:$count:$best_score:$best\n";
	    if ((!$final) && (!($final_hierarchical && ($best =~ /^\(/)))) {
		foreach my $choice (@choices)                                     # Go through each option of each column
		{
		    if (($best ne $choice) && (!($final_hierarchical && ($choice =~ /^\(/))) && ($column_scores{$contig}[$count]{$choice} >= ($fraction * $best_score))) {
			$columns{$contig}[$count] .= ';' . $choice;
			#print "$column_scores{$contig}[$count]{$choice}:$choice\n";
		    }
		}
	    }
	    #print "$columns{$contig}[$count]\n";
	    if ($best_score > $score_threshold) {
		$columns_status{$contig}->[$count] = 1;
	    } elsif ($best_score > $minimum_score) {
		if ($final && ($core_clusters[$best_clus] || ($best_score > ($minimum_score + 1)) || (($best_neighbor_score >= 1) && ($reduced_by_region[$best_index]->{'pid'} >= 90)))) {
		    $columns_status{$contig}->[$count] = 1;
		} else {
		    $columns_status{$contig}->[$count] = 0;
		}
	    } else {
		$columns_status{$contig}->[$count] = -1;
	    }
	    $count++;
	}
    }
    return;
}

###################################################################################################
sub score_neighbors
{
    # TODO - THIS WILL NOT WORK AS WRITTEN, SINCE NEIGHBORS DIFFERENTIATES BETWEEN 3' and 5'. Need to concatenate _3 or _5 on end to make it work
    # FOR NOW, DO BOTH. EVENTUALLY, DO BOTH, BUT APPLY A MODIFIER IF THE WRONG ONE
    # In addition to checking both orientations of the neighbor, we need to check for neighbors on both sides - so, need to ultimately check [1] and [3] in addition to [0] and [2]
    my $is_rev = shift;                                  #is the cluster on the forward or reverse strand
    my $is_left = shift;                                 #scoring the left side or the right
    my $contig = shift;                                  #contig id
    my $column = shift;                                  #column # within contig
    my $cluster = shift;                                 #cluster id
    my $blast = shift;                                   #blast match id
    my $is_hierarchical = shift;                         #is this assignment part of a hierarchical assignment
    my $nearhash = shift;                                #nearest style hash of clusters in the hierarchical set
    
    my $NEAR = $is_left ? $NEAR_5 : $NEAR_3;             #set index depending on left or right
    my $CORE = $is_left ? $CORE_5 : $CORE_3;             #set index depending on left or right
    
    # find and retain best neighbor on specified side
    my $best = 0;
    my $best_id = "";
    my $best_core = 0;
    my $best_core_dist = 0;
    my $best_core_tot = 0;
    my $best_id_core = "";
    my $PGG_NEAR = ($is_rev xor $is_left) ? $NEAR_5 : $NEAR_3; #determine which side of the cluster applies in this case
    my $PGG_CORE = ($is_rev xor $is_left) ? $CORE_5 : $CORE_3; #determine which side of the cluster applies in this case
    my $PGG_NEAR_REV = ($is_rev xor $is_left) ? $NEAR_3 : $NEAR_5; #determine which side of the cluster applies in this case
    my $PGG_CORE_REV = ($is_rev xor $is_left) ? $CORE_3 : $CORE_5; #determine which side of the cluster applies in this case
    if ($is_hierarchical)
    {
	foreach my $option (keys %{$neighbors[$cluster][$PGG_NEAR]})                             # go through all observed neighbors of cluster in PGG
	{
	    (my $clus, my $which_end) = split('_', $option);
	    if (defined($nearhash->{$clus}))                   # check if observed neighbor is a possible identity
	    {
		my $inv = $nearhash->{$clus};
		my $score = $neighbors[$cluster][$PGG_NEAR]{$option};
		if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		    $score /= 2; # decrease score if inverted
		}
		if ($score > $best)                       # if neighbor is observed, and is the best score observed
		{
		    $best = $score;                       # set it as new best score
		    $best_id = $option;                                             # retain identity
		}
	    }
	}
    }
    #print "NEAR";
    foreach my $option (keys %{$neighbors[$cluster][$PGG_NEAR]})                             # go through all observed neighbors of cluster in PGG
    {
	#print ": $option ";
	(my $clus, my $which_end) = split('_', $option);
	if (defined($nearest{$contig}[$column][$NEAR]{$clus}))                   # check if observed neighbor is a possible identity
	{
	    my $inv = $nearest{$contig}[$column][$NEAR]{$clus};
	    #print "$inv";
	    my $score = $neighbors[$cluster][$PGG_NEAR]{$option};
	    if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		$score /= 2; # decrease score if inverted
	    }
	    if ($score > $best)                       # if neighbor is observed, and is the best score observed
	    {
		$best = $score;                       # set it as new best score
		$best_id = $option;                                             # retain identity
	    }
	}
    }
    #print "\n";
    my $best_rev = 0;
    my $best_id_rev = "";
    foreach my $option (keys %{$neighbors[$cluster][$PGG_NEAR_REV]})                             # go through all observed neighbors of cluster in PGG
    {
	(my $clus, my $which_end) = split('_', $option);
	if (defined($nearest{$contig}[$column][$NEAR]{$clus}))                   # check if observed neighbor is a possible identity
	{
	    my $inv = $nearest{$contig}[$column][$NEAR]{$clus};
	    my $score = $neighbors[$cluster][$PGG_NEAR_REV]{$option};
	    if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		$score /= 2; # decrease score if inverted
	    }
	    if ($score > $best_rev)                       # if neighbor is observed, and is the best score observed
	    {
		$best_rev = $score;                       # set it as new best score
		$best_id_rev = $option;                                             # retain identity
	    }
	}
    }
    if ((2 * $best) < $best_rev) {
	$best = $best_rev / 2;
	$best_id = $best_id_rev
    }
    
    if ($is_hierarchical)
    {
	foreach my $option (keys %{$neighbors[$cluster][$PGG_CORE]})                             # go through all observed neighbors of cluster in PGG
	{
	    if ($option =~ /.*_d$/) {
		next; #skip the distance constraints
	    }
	    (my $clus, my $which_end) = split('_', $option);
	    if (defined($nearhash->{$clus}))                   # check if observed neighbor is a possible identity
	    {
		my $inv = $nearhash->{$clus};
		my $score = $neighbors[$cluster][$PGG_CORE]{$option};
		my $exp_dist = $neighbors[$cluster][$PGG_CORE]{$option . "_d"};
		my $dist = 1;
		my $dist_score = abs($exp_dist - $dist) / $exp_dist;
		if ($dist_score > .99) {
		    $dist_score = 0.01;
		} else {
		    $dist_score = 1 - $dist_score;
		}
		if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		    $score /= 2; # decrease score if inverted
		}
		my $tot_score = $score + $dist_score;
		if ($tot_score > $best_core_tot)                       # if neighbor is observed, and is the best score observed
		{
		    $best_core_tot = $tot_score;                       # set it as new best score
		    $best_core_dist = $dist_score;                       # set it as new best score
		    $best_core = $score;                       # set it as new best score
		    $best_id_core = $option;                                             # retain identity
		}
	    }
	}
    }
    #print "CORE";
    foreach my $option (keys %{$neighbors[$cluster][$PGG_CORE]})                             # go through all observed core neighbors of cluster in PGG
    {
	if ($option =~ /.*_d$/) {
	    next; #skip the distance constraints
	}
	#print ": $option ";
	(my $clus, my $which_end) = split('_', $option);
	if (defined($nearest{$contig}[$column][$CORE]{$clus}))                   # check if observed neighbor is a possible identity
	{
	    my $inv = $nearest{$contig}[$column][$CORE]{$clus};
	    my $score = $neighbors[$cluster][$PGG_CORE]{$option};
	    my $exp_dist = $neighbors[$cluster][$PGG_CORE]{$option . "_d"};
	    my $dist = $nearest{$contig}[$column][$CORE]{$clus . "_d"};
	    #print "$inv $score $exp_dist $dist";
	    my $dist_score = abs($exp_dist - $dist) / $exp_dist;
	    if ($dist_score > .99) {
		$dist_score = 0.01;
	    } else {
		$dist_score = 1 - $dist_score;
	    }
	    if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		    $score /= 2; # decrease score if inverted
	    }
	    my $tot_score = $score + $dist_score;
	    if ($tot_score > $best_core_tot)                       # if neighbor is observed, and is the best score observed
	    {
		$best_core_tot = $tot_score;                       # set it as new best score
		$best_core_dist = $dist_score;                       # set it as new best score
		$best_core = $score;                       # set it as new best score
		$best_id_core = $option;                                             # retain identity
	    }
	}
    }
    #print "\n";
    my $best_core_rev = 0;
    my $best_core_tot_rev = 0;
    my $best_core_dist_rev = 0;
    my $best_id_core_rev = "";
    foreach my $option (keys %{$neighbors[$cluster][$PGG_CORE_REV]})                             # go through all observed neighbors of cluster in PGG
    {
	if ($option =~ /.*_d$/) {
	    next; #skip the distance constraints
	}
	(my $clus, my $which_end) = split('_', $option);
	if (defined($nearest{$contig}[$column][$CORE]{$clus}))                   # check if observed neighbor is a possible identity
	{
	    my $inv = $nearest{$contig}[$column][$CORE]{$clus};
	    my $score = $neighbors[$cluster][$PGG_CORE_REV]{$option};
	    my $exp_dist = $neighbors[$cluster][$PGG_CORE_REV]{$option . "_d"};
	    my $dist = $nearest{$contig}[$column][$CORE]{$clus . "_d"};
	    my $dist_score = abs($exp_dist - $dist) / $exp_dist;
	    if ($dist_score > .99) {
		$dist_score = 0.01;
	    } else {
		$dist_score = 1 - $dist_score;
	    }
	    if (($is_left && (($inv && ($which_end == 3)) || (!$inv && ($which_end == 5)))) || (!$is_left && (($inv && ($which_end == 5)) || (!$inv && ($which_end == 3))))) {
		$score /= 2; # decrease score if inverted
	    }
	    my $tot_score = $score + $dist_score;
	    if ($tot_score > $best_core_tot_rev)                       # if neighbor is observed, and is the best score observed
	    {
		$best_core_tot_rev = $tot_score;                       # set it as new best score
		$best_core_dist_rev = $dist_score;                       # set it as new best score
		$best_core_rev = $score;                       # set it as new best score
		$best_id_core_rev = $option;                                             # retain identity
	    }
	}
    }
    if ($best_core_tot < (($best_core_rev / 2) + $best_core_dist_rev)) {
	$best_core_tot = (($best_core_rev / 2) + $best_core_dist_rev);
	$best_core_dist = $best_core_dist_rev;
	$best_core = $best_core_rev / 2;
	$best_id_core = $best_id_core_rev
    }
    
    return ($best, $best_id, $best_core, $best_id_core, $best_core_dist);
}

###################################################################################################
sub homology_score
# calculate homology based protion of score
{
    (my $clus, my $index)  = @_;
    my $per_len = (($reduced_by_region[$index]->{'qend'} - $reduced_by_region[$index]->{'qbeg'}) + 1) / $reduced_by_region[$index]->{'qlen'};
    my $per_id = $reduced_by_region[$index]->{'pid'};
    my $coreness = $cluster_size[$clus] / $num_genomes;
    my $uniqueness = 1 / $cluster_hits[$clus];
    my $single_copy = $core_clusters[$clus];
    if ($per_len >= .95) {
	$per_len = 1;
    } elsif ($per_len >= .90) {
	$per_len = .75;
    } elsif ($per_len >= .70) {
	$per_len = .5;
    } elsif ($per_len >= .50) {
	$per_len = .25;
    } else {
	$per_len = 0;
    }
    if ($per_id >= .95) {
	$per_id = 1;
    } elsif ($per_id >= .90) {
	$per_id = .75;
    } elsif ($per_id >= .70) {
	$per_id = .5;
    } elsif ($per_id >= .50) {
	$per_id = .25;
    } else {
	$per_id = 0;
    }
    #print "$single_copy:$coreness:$uniqueness:$per_id:$per_len\n";
    return($single_copy + $coreness + $uniqueness + $per_id + $per_len);
}

###################################################################################################
sub buffer_align_length
# determine extra amount of sequence for Needleman-Wunsch alignment
{
    my $length = shift; #the unbufferd alignment length
    if ($length < 50) {
	$length = sprintf("%d", ((2 * $length) + 1));
    } elsif ($length < 100) {
	$length =  sprintf("%d", ((1.5 * $length) + 26));
    } elsif ($length < 1000) {
	$length =  sprintf("%d", ((1.2 * $length) + 56));
    } else {
	$length =  sprintf("%d", ((1.1 * $length) + 156));
    }
    return ($length);
}

###################################################################################################
sub refine_alignments
# output the newly calculated attribute file and matchtable column file
{
    foreach my $contig (keys %columns) {
	my $col_index = 0;
	foreach my $column (@{ $columns{$contig} }) {
	    if ($columns_status{$contig}->[$col_index] >= 0) {
		my $beg_align = -1;
		my $end_align = -1;
		my $ignore = -1;
		my $num_matches = -1;
		my $s_offset = -1;
		my $s_extra = -1;
		my $sequence = "";
		my $contig_sequence = "";
		(my $clus, my $index) = split('_', $column);
		if (!$reduced_by_region[$index]->{'sinv'} && ($reduced_by_region[$index]->{'qbeg'} != 1)) {
		    $s_extra = ($reduced_by_region[$index]->{'sbeg'} - 1) + ($align_anchor_len - 1); # give $align_anchor_len bp to anchor the alignment but only really looking for the beginning need -1 for string beginning at 0 and 'beg" starting at 1
		    if ($s_extra >= $reduced_by_region[$index]->{'ctglen'}) {$s_extra = $reduced_by_region[$index]->{'ctglen'} - 1;}
		    $s_offset = ($reduced_by_region[$index]->{'sbeg'} - buffer_align_length($reduced_by_region[$index]->{'qbeg'} - 1)) - 1; # last -1 is for strings being indexed starting at 0 and contigs coordiantes starting at 1
		    if ($s_offset < 0) {
			if ($is_circular{$contig}) {
			    $contig_sequence = substr($contigs{$contig}, $s_offset);
			    $contig_sequence .= substr($contigs{$contig}, 0, ($s_extra  + 1));
			} else {
			    $s_offset = 0;
			    $contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
			}
		    } else {
			$contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
		    }
		    $sequence = substr($medoids[$clus], 0, (($reduced_by_region[$index]->{'qbeg'} - 1) + $align_anchor_len));
		    ($beg_align, $ignore, $num_matches) = &nw_align($sequence, $contig_sequence);
		    $beg_align += $s_offset;
		    #print "NWB:$beg_align:$ignore:$num_matches:$s_offset:$s_extra\n";
		} elsif ($reduced_by_region[$index]->{'sinv'} && ($reduced_by_region[$index]->{'qlen'} != $reduced_by_region[$index]->{'qend'})) {
		    $s_extra = ($reduced_by_region[$index]->{'sbeg'} - 1) + ($align_anchor_len - 1); # give $align_anchor_len bp to anchor the alignment but only really looking for the beginning need -1 for string beginning at 0 and 'beg" starting at 1
		    if ($s_extra >= $reduced_by_region[$index]->{'ctglen'}) {$s_extra = $reduced_by_region[$index]->{'ctglen'} - 1;}
		    $s_offset = ($reduced_by_region[$index]->{'sbeg'} - buffer_align_length($reduced_by_region[$index]->{'qlen'} - $reduced_by_region[$index]->{'qend'})) - 1; # last -1 is for strings being indexed starting at 0 and contigs coordiantes starting at 1
		    if ($s_offset < 0) {
			if ($is_circular{$contig}) {
			    $contig_sequence = substr($contigs{$contig}, $s_offset);
			    $contig_sequence .= substr($contigs{$contig}, 0, ($s_extra  + 1));
			} else {
			    $s_offset = 0;
			    $contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
			}
		    } else {
			$contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
		    }
		    $sequence = substr(reverse($medoids[$clus]), 0, ($reduced_by_region[$index]->{'qlen'} - $reduced_by_region[$index]->{'qend'}) + $align_anchor_len);
		    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    ($beg_align, $ignore, $num_matches) = &nw_align($sequence, $contig_sequence);
		    $beg_align += $s_offset;
		    #print "NWBC:$beg_align:$ignore:$num_matches:$s_offset:$s_extra\n";
		} else {
		    $beg_align = $reduced_by_region[$index]->{'sbeg'};
		}
		if (!$reduced_by_region[$index]->{'sinv'} && ($reduced_by_region[$index]->{'qlen'} != $reduced_by_region[$index]->{'qend'})) {
		    $s_offset = $reduced_by_region[$index]->{'send'} - $align_anchor_len; # give $align_anchor_len bp to anchor the alignment but only really looking for the ending
		    if ($s_offset < 0) {$s_offset = 0;}
		    $s_extra = ($reduced_by_region[$index]->{'send'} + buffer_align_length($reduced_by_region[$index]->{'qlen'} - $reduced_by_region[$index]->{'qend'})) - 1; # last -1 is for strings being indexed starting at 0 and contigs coordiantes starting at 1
		    if ($s_extra >= $reduced_by_region[$index]->{'ctglen'}) {
			if ($is_circular{$contig}) {
			    $contig_sequence = substr($contigs{$contig}, $s_offset);
			    $contig_sequence .= substr($contigs{$contig}, 0, (($s_extra  + 1) - $reduced_by_region[$index]->{'ctglen'}));
			} else {
			    $s_extra = $reduced_by_region[$index]->{'ctglen'} - 1;
			    $contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
			}
		    } else {
			$contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
		    }
		    $sequence = substr($medoids[$clus], (-1 *(($reduced_by_region[$index]->{'qlen'} - $reduced_by_region[$index]->{'qend'}) + $align_anchor_len)));
		    ($ignore, $end_align, $num_matches) = &nw_align($sequence, $contig_sequence);
		    $end_align += $s_offset;
		    #print "NWE:$ignore:$end_align:$num_matches:$s_offset:$s_extra\n";
		} elsif ($reduced_by_region[$index]->{'sinv'} && ($reduced_by_region[$index]->{'qbeg'} != 1)) {
		    $s_offset = $reduced_by_region[$index]->{'send'} - $align_anchor_len; # give $align_anchor_len bp to anchor the alignment but only really looking for the ending
		    if ($s_offset < 0) {$s_offset = 0;}
		    $s_extra = ($reduced_by_region[$index]->{'send'} + buffer_align_length($reduced_by_region[$index]->{'qbeg'} - 1)) - 1; # last -1 is for strings being indexed starting at 0 and contigs coordiantes starting at 1
		    if ($s_extra >= $reduced_by_region[$index]->{'ctglen'}) {
			if ($is_circular{$contig}) {
			    $contig_sequence = substr($contigs{$contig}, $s_offset);
			    $contig_sequence .= substr($contigs{$contig}, 0, (($s_extra  + 1) - $reduced_by_region[$index]->{'ctglen'}));
			} else {
			    $s_extra = $reduced_by_region[$index]->{'ctglen'} - 1;
			    $contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
			}
		    } else {
			$contig_sequence = substr($contigs{$contig}, $s_offset, (($s_extra - $s_offset) + 1));
		    }
		    $sequence = substr(reverse($medoids[$clus]), (-1 *(($reduced_by_region[$index]->{'qbeg'} - 1) + $align_anchor_len)));
		    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    ($ignore, $end_align, $num_matches) = &nw_align($sequence, $contig_sequence);
		    $end_align += $s_offset;
		    #print "NWEC:$ignore:$end_align:$num_matches:$s_offset:$s_extra\n";
		} else {
		    $end_align = $reduced_by_region[$index]->{'send'};
		}
		#print "COL$col_index:$column:$columns_status{$contig}->[$col_index]:$reduced_by_region[$index]->{'sinv'}\n";
		#print "NWFINAL:$beg_align:$end_align:$num_matches:$s_offset:$s_extra\n";
		#print "BLAST($column_scores{$contig}[$col_index]->{'best_score'}) $column:$reduced_by_region[$index]->{'clus'}($cluster_size[$reduced_by_region[$index]->{'clus'}]):$reduced_by_region[$index]->{'ctg'}:$reduced_by_region[$index]->{'pid'}:$reduced_by_region[$index]->{'qbeg'}:$reduced_by_region[$index]->{'qend'}:$reduced_by_region[$index]->{'qlen'}:$reduced_by_region[$index]->{'sbeg'}:$reduced_by_region[$index]->{'send'}:$reduced_by_region[$index]->{'sinv'}:$reduced_by_region[$index]->{'ctglen'}:$reduced_by_region[$index]->{'bits'}:$reduced_by_region[$index]->{'keepclus'}:$reduced_by_region[$index]->{'keepctg'}:$reduced_by_region[$index]->{'weak'}\n";
		$reduced_by_region[$index]->{'sbeg'} = $beg_align;
		$reduced_by_region[$index]->{'send'} = $end_align;
	    }
	}
    }
    return;
}

###################################################################################################
sub output_files
# output the newly calculated attribute file and matchtable column file
{
    foreach my $cluster (1 .. $num_clusters) {
	$match_table[$cluster] = "----------"; # initialize for missing clusters
    }
    my $wgsANIfile;
    unless (open ($wgsANIfile, ">", ($rootname . "_wgsANI.txt")) )  {
	die ("cannot open whole genome ANI file: $rootname'_wgsANI.txt'!\n");
    }
    my $geneANIfile;
    unless (open ($geneANIfile, ">", ($rootname . "_geneANI.txt")) )  {
	die ("cannot open whole genome ANI file: $rootname'_geneANI.txt'!\n");
    }
    my $rearrange;
    unless (open ($rearrange, ">", ($rootname . "_rearrange.txt")) )  {
	die ("cannot open rearrangementss file: $rootname'_rearrange.txt'!\n");
    }
    my $alledgesfile;
    if ($reannotate) {
	unless (open ($alledgesfile, ">", ($rootname . "_alledges.txt")) )  {
	    die ("cannot open all edges file: $rootname'_alledges.txt'!\n");
	}
    }
    my $attributefile;
    unless (open ($attributefile, ">", ($rootname . "_attributes.txt")) )  {
	die ("cannot open attributes file: $rootname'_attributes.txt'!\n");
    }
    my $new_attributefile;
    unless (open ($new_attributefile, ">", ($rootname . "_attributes_new.txt")) )  {
	die ("cannot open new attributes file: $rootname'_attributes_new.txt'!\n");
    }
    my $uniqclusfile;
    unless (open ($uniqclusfile, ">", ($rootname . "_uniq_clus.txt")) )  {
	die ("cannot open new unique clusters file: $rootname'_uniq_clus.txt'!\n");
    }
    my $sumANI = 0;
    my $sumANIlen = 0;
    my $paralog_index = $num_clusters + 1;
    foreach my $contig (keys %columns) {
	my $col_index = 0;
	my $first_clus_orient = "";
	my $prev_clus_orient = "";
	my $prev_all_clus_orient = "";
	my $clus_orient = "";
	my $next_clus_orient = "";
	my $core_orient = "";
	my $prev_core_orient = "";
	my $next_core_orient = "";
	my $core_edge_start = 0;
	my $core_edge_end = 0;
	my $edge_start = 0;
	my $edge_end = 0;
	my $edge_all_start = 0;
	foreach my $column (@{ $columns{$contig} }) {
	    if ($columns_status{$contig}->[$col_index] >= 0) {
		(my $clus, my $index) = split('_', $column);
		my $beg_align = $reduced_by_region[$index]->{'sbeg'};
		my $end_align = $reduced_by_region[$index]->{'send'};
		#print "COL$col_index:$column:$columns_status{$contig}->[$col_index]:$reduced_by_region[$index]->{'sinv'}\n";
		#print "BLAST($column_scores{$contig}[$col_index]->{'best_score'}) $column:$reduced_by_region[$index]->{'clus'}($cluster_size[$reduced_by_region[$index]->{'clus'}]):$reduced_by_region[$index]->{'ctg'}:$reduced_by_region[$index]->{'pid'}:$reduced_by_region[$index]->{'qbeg'}:$reduced_by_region[$index]->{'qend'}:$reduced_by_region[$index]->{'qlen'}:$reduced_by_region[$index]->{'sbeg'}:$reduced_by_region[$index]->{'send'}:$reduced_by_region[$index]->{'sinv'}:$reduced_by_region[$index]->{'ctglen'}:$reduced_by_region[$index]->{'bits'}:$reduced_by_region[$index]->{'keepclus'}:$reduced_by_region[$index]->{'keepctg'}:$reduced_by_region[$index]->{'weak'}\n";
		if ($columns_status{$contig}->[$col_index] == 1) {
		    my $ortholog = $target . "CL_" . $clus;
		    $clus_orient = $clus . '_' . ($reduced_by_region[$index]->{'sinv'} ? '3' : '5');
		    $edge_end = $beg_align - 1;
		    $next_clus_orient = $clus . '_' . ($reduced_by_region[$index]->{'sinv'} ? '5' : '3');
		    if ($columns_core{$contig}[$col_index]->{'is_core'}) {
			$core_orient = $clus . '_' . ($reduced_by_region[$index]->{'sinv'} ? '3' : '5');
			$core_edge_end = $beg_align - 1;
			$next_core_orient = $clus . '_' . ($reduced_by_region[$index]->{'sinv'} ? '5' : '3');
			if ($prev_core_orient ne "") {
			    if (!defined($pgg_core_edges{'(' . $prev_core_orient . ',' . $core_orient . ')'})) {
				print $rearrange "$target\t$contig\trearrange\t$core_edge_start\t$core_edge_end\t", (($core_edge_end - $core_edge_start) + 1), "\t", ('(' . $prev_core_orient . ',' . $core_orient . ')'). "\n";
			    }
			}
			$prev_core_orient = $next_core_orient;
			$core_edge_start = $end_align + 1;
		    }
		    if ($prev_clus_orient ne "") {
			my $tmp_len = ($edge_end - $edge_start) + 1;
			$hash_edges{'(' . $prev_clus_orient . ',' . $clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_start\t$edge_end\t$tmp_len\t";
			$hash_edges{'(' . $clus_orient . ',' . $prev_clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_end\t$edge_start\t$tmp_len\t";
			if ($reannotate) {
			    print $alledgesfile '(' . $prev_clus_orient . ',' . $clus_orient . ')', "\n";
			    print $alledgesfile '(' . $clus_orient . ',' . $prev_clus_orient . ')', "\n";
			}
		    }
		    if ($reannotate && ($prev_all_clus_orient ne $prev_clus_orient)) {
			my $tmp_len = ($edge_end - $edge_all_start) + 1;
			$hash_edges{'(' . $prev_all_clus_orient . ',' . $clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_all_start\t$edge_end\t$tmp_len\t";
			$hash_edges{'(' . $clus_orient . ',' . $prev_all_clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_end\t$edge_all_start\t$tmp_len\t";
		    }
		    $prev_clus_orient = $next_clus_orient;
		    $edge_start = $end_align + 1;
		    $prev_all_clus_orient = $next_clus_orient;
		    $edge_all_start = $end_align + 1;
		    if ($first_clus_orient eq "") {
			$first_clus_orient = $clus_orient;
		    }
		    $match_table[$clus] = $ortholog;
		    my $tmp_len = ($end_align - $beg_align) + 1;
		    print $attributefile "$contig\t$ortholog\t";
		    if ($reduced_by_region[$index]->{'sinv'}) {
			print $attributefile "$end_align\t$beg_align\t";
			if ($reduced_by_region[$index]->{'pid'} < 90) {
			    print $geneANIfile "$target\t$contig\tgeneANI\t$end_align\t$beg_align\t$tmp_len\t$ortholog\n";
			}
		    } else {
			print $attributefile "$beg_align\t$end_align\t";
			if ($reduced_by_region[$index]->{'pid'} < 90) {
			    print $geneANIfile "$target\t$contig\tgeneANI\t$beg_align\t$end_align\t$tmp_len\t$ortholog\n";
			}
		    }
		    print $attributefile "$medoids_anno[$clus]\t$target\n";
		    my $matchlen = ($reduced_by_region[$index]->{'qend'} - $reduced_by_region[$index]->{'qbeg'}) + 1;
		    $sumANI += $reduced_by_region[$index]->{'pid'} * $matchlen;
		    $sumANIlen += $matchlen;
		} else {
		    my $paralog = $target . "PL_" . $paralog_index . "_" . $clus;
		    $clus_orient = $paralog_index . '_' . ($reduced_by_region[$index]->{'sinv'} ? '3' : '5');
		    $edge_end = $beg_align - 1;
		    $next_clus_orient = $paralog_index . '_' . ($reduced_by_region[$index]->{'sinv'} ? '5' : '3');
		    $paralog_index++;
		    if ($reannotate) {
			if ($prev_all_clus_orient ne "") {
			    my $tmp_len = ($edge_end - $edge_all_start) + 1;
			    $hash_edges{'(' . $prev_all_clus_orient . ',' . $clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_all_start\t$edge_end\t$tmp_len\t";
			    $hash_edges{'(' . $clus_orient . ',' . $prev_all_clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_end\t$edge_all_start\t$tmp_len\t";
			}
			$prev_all_clus_orient = $next_clus_orient;
			$edge_all_start = $end_align + 1;
		    } else {
			if ($prev_clus_orient ne "") {
			    my $tmp_len = ($edge_end - $edge_start) + 1;
			    $hash_edges{'(' . $prev_clus_orient . ',' . $clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_start\t$edge_end\t$tmp_len\t";
			    $hash_edges{'(' . $clus_orient . ',' . $prev_clus_orient . ')'} = "$target\t$contig\tuniq_edge\t$edge_end\t$edge_start\t$tmp_len\t";
			}
			$prev_clus_orient = $clus_orient;
			$edge_start = $end_align + 1;
			if ($first_clus_orient eq "") {
			    $first_clus_orient = $clus_orient;
			}
		    }
		    push @new_match_table, $paralog;
		    my $tmp_len = ($end_align - $beg_align) + 1;
		    print $new_attributefile "$contig\t$paralog\t";
		    if ($reduced_by_region[$index]->{'sinv'}) {
			print $new_attributefile "$end_align\t$beg_align\t";
			print $uniqclusfile "$target\t$contig\tuniq_clus\t$end_align\t$beg_align\t$tmp_len\t$paralog\n";
		    } else {
			print $new_attributefile "$beg_align\t$end_align\t";
			print $uniqclusfile "$target\t$contig\tuniq_clus\t$beg_align\t$end_align\t$tmp_len\t$paralog\n";
		    }
		    print $new_attributefile "$medoids_anno[$clus]\t$target\n";
		}
	    }
	    $col_index++;
	}
	if (($is_circular{$contig}) && ($first_clus_orient ne "") && ($prev_clus_orient ne "")) {
	    $hash_edges{'(' . $prev_clus_orient . ',' . $first_clus_orient . ')'} = 0;
	    $hash_edges{'(' . $first_clus_orient . ',' . $prev_clus_orient . ')'} = 0;
	}
	
    }
    printf $wgsANIfile "%5.2f\n", ($sumANI / $sumANIlen);
    close ($wgsANIfile);
    close ($geneANIfile);
    close ($rearrange);
    if ($reannotate) {
	close ($alledgesfile);
    }
    close ($attributefile);
    close ($new_attributefile);
    close ($uniqclusfile);
    my $matchfile;
    unless (open ($matchfile, ">", ($rootname . "_match.col")) )  {
	die ("cannot open matchtable column file: $rootname'_match.col'!\n");
    }
    #print "Outputting column for matchtable to $rootname'_match.col'\n";
    foreach my $cluster (1 .. $num_clusters) {
	print $matchfile $match_table[$cluster], "\n";
    }
    close ($matchfile);
    my $new_matchfile;
    unless (open ($new_matchfile, ">", ($rootname . "_match_new.col")) )  {
	die ("cannot open matchtable column file: $rootname'_match_new.col'!\n");
    }
    my $paralog_index = $num_clusters + 1;
    foreach my $cluster (@new_match_table) {
	print $new_matchfile "$paralog_index\t$cluster\n";
	$paralog_index++;
    }
    close ($new_matchfile);
    my $edgefile;
    unless (open ($edgefile, ">", ($rootname . "_pgg.col")) )  {
	die ("cannot open edges column file: $rootname'_pgg.col'!\n");
    }
    foreach my $edge (@edges) {
	if (defined $hash_edges{$edge}) {
	    $hash_edges{$edge} = "NOT_UNIQUE";
	    print $edgefile "1\n";
	} else {
	    print $edgefile "0\n";
	}
    }
    close ($edgefile);
    my $new_edgefile;
    unless (open ($new_edgefile, ">", ($rootname . "_pgg_new.col")) )  {
	die ("cannot open edges column file: $rootname'_pgg_new.col'!\n");
    }
    my $uniqedgefile;
    unless (open ($uniqedgefile, ">", ($rootname . "_uniq_edge.txt")) )  {
	die ("cannot open unique edge file: $rootname'_uniq_edge.txt'!\n");
    }
    foreach my $edge (sort (keys %hash_edges)) {
	if ($hash_edges{$edge} ne "NOT_UNIQUE") {
	    print $new_edgefile "$edge\t1\n";
	    my $cluster1;
	    my $cluster2;
	    if ($edge =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
		$cluster1 = $1;
		$cluster2 = $3;
	    } else {
		die ("ERROR: Bad edge formatting $edge in hash_edges\n");
	    }
	    if ($cluster1 < $cluster2) { #only capture one representation of the symmetric edge
		print $uniqedgefile "$hash_edges{$edge}$edge\n";
	    }
	}
    }
    close ($new_edgefile);
    close ($uniqedgefile);
    return;
}

###################################################################################################
sub assign_paralogs
# output the newly calculated attribute file and matchtable column file
{
    foreach my $contig (keys %columns) {
	my $col_index = 0;
	foreach my $column (@{ $columns{$contig} }) {
	    (my $clus, my $index) = split('_', $column);
	    #print "COL$col_index:$column:$columns_status{$contig}->[$col_index]\n";
	    #print "BLAST($column_scores{$contig}[$col_index]->{'best_score'}) $column:$reduced_by_region[$index]->{'clus'}($cluster_size[$reduced_by_region[$index]->{'clus'}]):$reduced_by_region[$index]->{'ctg'}:$reduced_by_region[$index]->{'pid'}:$reduced_by_region[$index]->{'qbeg'}:$reduced_by_region[$index]->{'qend'}:$reduced_by_region[$index]->{'qlen'}:$reduced_by_region[$index]->{'sbeg'}:$reduced_by_region[$index]->{'send'}:$reduced_by_region[$index]->{'sinv'}:$reduced_by_region[$index]->{'ctglen'}:$reduced_by_region[$index]->{'bits'}:$reduced_by_region[$index]->{'keepclus'}:$reduced_by_region[$index]->{'keepctg'}:$reduced_by_region[$index]->{'weak'}\n";
	    if ($columns_status{$contig}->[$col_index] == 1) {
		#print "$contig $col_index != $cluster_colindex[$clus]\n";
		if (($contig . "_" . $col_index) ne $cluster_colindex[$clus]) {
		    $columns_status{$contig}->[$col_index] = 2;
		}
	    }
	    $col_index++;
	}
    }
    return;
}

###################################################################################################
sub check_overlaps
# check to see if columns overlap too much and if so if one match is weak enough to disregard
{
    foreach my $contig (keys %columns) {
	foreach my $i (0 .. $#{ $columns{$contig} }) {
	    (my $clus1, my $index1) = split('_', $columns{$contig}[$i]);
	    my $beg1 = $reduced_by_region[$index1]->{'sbeg'};
	    my $end1 = $reduced_by_region[$index1]->{'send'};
	    my $bits1 = $reduced_by_region[$index1]->{'bits'};
	    my $len1 = ($end1 - $beg1) + 1;
	    foreach my $j (($i + 1) .. $#{ $columns{$contig} }) {
		(my $clus2, my $index2) = split('_', $columns{$contig}[$j]);
		my $beg2 = $reduced_by_region[$index2]->{'sbeg'};
		my $end2 = $reduced_by_region[$index2]->{'send'};
		my $bits2 = $reduced_by_region[$index2]->{'bits'};
		my $len2 = ($end2 - $beg2) + 1;
		my $overlap;
		if ($beg2 >= $end1) {
		    last;
		}
		if ($end2 < $end1) {
		    $overlap = $len2;
		} else {
		    $overlap = ($end1 - $beg2) + 1;
		}
		my $cov1 = $overlap / $len1;
		my $cov2 = $overlap / $len2;
		my $maxcov = $cov1 > $cov2 ? $cov1 : $cov2;
		if (($maxcov > 0.5) || ($overlap >= 50)) {
		    if (($columns_status{$contig}->[$i] <= 0) && ($columns_status{$contig}->[$j] > 0)) {
			$columns_status{$contig}->[$i] = -1;
		    } elsif (($columns_status{$contig}->[$j] <= 0) && ($columns_status{$contig}->[$i] > 0)) {
			$columns_status{$contig}->[$j] = -1;
		    } elsif (($columns_status{$contig}->[$i] > 0) && ($columns_status{$contig}->[$j] > 0)) {
			if ($maxcov > 0.5) {
			    if ($columns_status{$contig}->[$i] > $columns_status{$contig}->[$j]) {
				$columns_status{$contig}->[$i] = -1;
			    } elsif ($columns_status{$contig}->[$i] < $columns_status{$contig}->[$j]) {
				$columns_status{$contig}->[$j] = -1;
			    } elsif ($column_scores{$contig}[$i]->{'best_score'} < $column_scores{$contig}[$j]->{'best_score'}) {
				$columns_status{$contig}->[$i] = -1;
			    } elsif ($column_scores{$contig}[$i]->{'best_score'} > $column_scores{$contig}[$j]->{'best_score'}) {
				$columns_status{$contig}->[$j] = -1;
			    } elsif ($bits1 > $bits2) {
				$columns_status{$contig}->[$j] = -1;
			    } elsif ($bits2 > $bits1) {
				$columns_status{$contig}->[$i] = -1;
			    } elsif ($len1 > $len2) {
				$columns_status{$contig}->[$j] = -1;
			    } else {
				$columns_status{$contig}->[$i] = -1;
			    }
			}
		    } elsif (($columns_status{$contig}->[$i] == 0) && ($columns_status{$contig}->[$j] == 0)) {
			if ($column_scores{$contig}[$i]->{'best_score'} < $column_scores{$contig}[$j]->{'best_score'}) {
			    $columns_status{$contig}->[$i] = -1;
			} else {
			    $columns_status{$contig}->[$j] = -1;
			}
		    } else {
		    }
		}
	    }
	}
    }
    return;
}

###################################################################################################
{#main
    print "Reading PGG\n";
    &read_pgg;
    print "Reading genomes\n"
    &read_genome;
    print "Reading medoids\n";
    &read_medoids;
    print "Reading clusters sizes\n";
    &read_cluster_sizes;
    print "Reading single copy ocres\n";
    &read_core_list;
    print "Reading Blast matches\n";
    &read_blast;
    print "Processing Blast by clusters\n";
    &process_blast_by_cluster;
    print "Processing Blast by region\n";
    &process_blast_by_region;
    print "Reading PGG neighbors\n";
    &read_neighbors;
    print "Determining genome neighbors\n";
    &determine_contig_nearest_neighbors;
    print "Scoring columns\n";
    &score_from_neighbors(0, 0, 0.8);
    print "Calculating median score\n";
    &calc_median;
    print "Resetting columns\n";
    &reset_columns;
    print "Iterating scoring\n";
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.8);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.85);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.85);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.9);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.9);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 0, 0.9);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 1, 0.95);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(0, 1, 0.95);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(1, 1, 0.95);
    &reset_columns;
    &determine_contig_nearest_neighbors;
    &score_from_neighbors(1, 1, 0.95);
    print "Assigning paralogs\n";
    &assign_paralogs;
    print "Checking overlapping columns\n";
    &check_overlaps;
    print "Outputting files\n";
    &refine_alignments;
    &check_overlaps;
    &output_files;
}