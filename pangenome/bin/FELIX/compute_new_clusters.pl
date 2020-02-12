#!/usr/bin/env perl
#Copyright (C) 2011-2015 The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D. and Granger Sutton, Ph.D.

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

my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use Cwd;
use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
getopts ('hDc:n:b:i:M:m:g:s:B:L:');
our ($opt_h,$opt_D,$opt_c,$opt_n,$opt_i,$opt_M,$opt_m,$opt_g,$opt_s,$opt_B,$opt_L);

my ($cluster_file,$start_cluster,$index_file,$match_file,$medoid_file,$gtag_file,$seq_file);
my $blast_directory = "";
my $ld_load_directory = "";
#set defaults

## use boolean logic:  TRUE = 1, FALSE = 0
my $DEBUG;
my $version = "1.0";
my $cwd = getcwd;
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_i) {$index_file = $opt_i;} else { $index_file = "gene_cluster_index.txt"; } # Read in index_file name or set default.
if ($opt_M) {$match_file = $opt_M;} else { $match_file = "new_matchtable.txt"; } # Read in match_file name or set default.
if ($opt_m) {$medoid_file = $opt_m;} else { $medoid_file = "new_medoids.fasta"; } # Read in medoid_file name or set default.
if (($opt_c) && (-s $opt_c)) {$cluster_file = $opt_c;} else { print STDERR "Error with -c $opt_c\n"; &option_help; } # if no value for option c (cluster_file), quit with help menu
if (($opt_g) && (-s $opt_g)) {$gtag_file = $opt_g;} else { print STDERR "Error with -g $opt_g\n"; &option_help; } # if no value for option g (gtab_file), quit with help menu
if (($opt_s) && (-s $opt_s)) {$seq_file = $opt_s;} else { print STDERR "Error with -s $opt_s\n"; &option_help; } # if no value for option s (seq_file), quit with help menu
if ($opt_B) {
    if (-d $opt_B) {
	$blast_directory = $opt_B;
	if (substr($blast_directory, -1, 1) ne "/") {
	    $blast_directory .= "/";
	}
	if (substr($blast_directory, 0, 1) ne "/") {
	    $blast_directory = $cwd . "/$blast_directory";
	}
    } else {
	print STDERR "Error with -B $opt_B\n";
	&option_help;
    }
}
if ($opt_L) {
    if (-d $opt_L) {
	$ld_load_directory = $opt_L;
	if (substr($ld_load_directory, -1, 1) ne "/") {
	    $ld_load_directory .= "/";
	}
	if (substr($ld_load_directory, 0, 1) ne "/") {
	    $ld_load_directory = $cwd . "/$ld_load_directory";
	}
	$blast_directory = 'export LD_LIBRARY_PATH=' . $ld_load_directory . ':$LD_LIBRARY_PATH; ' . $blast_directory;
    } else {
	print STDERR "Error with -L $opt_L\n";
	&option_help;
    }
}
if ($opt_n) {
    $start_cluster = $opt_n;
    if (($start_cluster !~ /^\d+$/) || ($start_cluster <= 0)){
	die ("ERROR: $start_cluster is not a positive integer for -n!\n");
    }
} else { print STDERR "Error with -n (no starting cluster number given\n"; &option_help; } # if no value for option n (start_cluster), quit with help menu

my %feat_hash = ();         # Key = feat_name value = genome tag
my %BestTagMatch = ();      # Key1 = Query protein, Key2 = subject genome tag, Value = bit score of best match for this genome
my %SecondBestTagMatch = ();# Key1 = Query protein, Key2 = subject genome tag, Value = bit score of best match for this genome
my %relationship_hash = (); # Key1 = Query protein, Key2 = Subject protein, Key3 = struct members with their values
my %Qbytaghash = ();        # Key1 = Query protein, Key2 = subject genome tag, Key3 = Subject protein, Key4 = same as relationship_hash Key3 (actually is the same pointer)
my %Tagbyfeatnamehash = (); # Key1 = genome tag, Key2 = feat_name (query)
my %TagByPointer = ();      # Key = feat_name-tag, value = array position location
my @tag_array = ();         # array of genome tags
my %FeatnameLookupTag_hash = ();# Key = feat_name, Value = tag.  Needed for input lacking tagged feat_names
my %TagIndex = ();    # Key = genome tag, Value = index in tag_array of associated genome tag
my %clusters = ();          # Key = feat_name, Value = array of (hashes of feat_name and genome_tag) that are in a cluster with the feat_name key
my %cluster_number = ();    # Key = feat_name, Value = cluster number that feat_name is in
my %cluster_size = ();      # Key = cluster number, Value = size of the cluster
my @cluster_for_number =(); # Index = cluster number, Value = cluster structure (array)
my @medoids = ();          # Index = cluster number, Value = feature id for the medoid of the cluster
my $genome_number;          # number of genomes defined in tagfile

#%relationship_hash stores information for the blast matches
###key1 = query
###key2 = subject
###value -> id = percent identity
###         score = BLAST bits score
###         best = 1 if best blast match in the subject genome, otherwise 0
###         bibest = 1 if best bidirectional blast match in the subject genome, otherwise 0
###         neighboirscore = score based on shared neighbors

sub get_tags {  # obtain list of genomes to compare
   
    my %temp_hash = ();
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", $gtag_file) || die ("ERROR: can't open file $gtag_file\n");
    while (<$infile>)  {
	chomp;
	if (length($_) > 0) {
	    my $name = $_;
	    $name =~ s/\s+$//;

            if (defined $temp_hash{$name})  {
               die ("ERROR:  You have more than one occurance of $name in /$gtag_file!\n");
            }
	    $temp_hash{$name} = 1;
	    push (@tag_array, $name); # populate the tag_array in the order of the tagfile (1st tag is the reference tag)
	    $genome_number++;
	    print STDERR "$name\n" if ($DEBUG);
	}  
    }
    close($infile);

    my $index = 0;
    foreach my $tag (@tag_array) {
	$TagIndex{$tag} = $index++
    }
}

sub get_seqs { #read in the fasta sequences which were compared with blast

  my @line = ();
  my $id;
  my $title = "";
  my $sequence = "";
  my $length = "";

  unless (open (PEPFILE, "<", $seq_file) )  {
    die ("can't open file $seq_file.\n");
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
    $feat_hash{$id}->{'length'} = $length;
    $feat_hash{$id}->{'sequence'} = $sequence;
    #print STDERR "$id ($feat_hash{$id}->{'header'}) = $feat_hash{$id}->{'length'}\n";
    $title = ""; # clear the title for the next round.
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (PEPFILE);
  return;
}

sub get_gene_att {

    my $tag = "";
    my $feat_name = "";
    my $failed = 0;

    unless (open (ATTFILE, "<", $cluster_file) )  {
	die ("ERROR: can not open file $cluster_file.\n");
    }
    while (<ATTFILE>) {
	my @att_line = ();
	chomp;
	@att_line = split(/\t/, $_,6);  # split the scalar $line on tab
	$feat_name = $att_line[1];
	if (defined $feat_hash{$feat_name}->{'near5'}) {
	    print STDERR "ERROR: $feat_name appears more than once in the cluster file!\n";
	    $failed = 1;
	}
	$tag = $att_line[0];
	if (!defined $TagIndex{$tag}) {
	    print STDERR "ERROR: $tag is a genome tag in the cluster file but not in the genome tag file!\n";
	    $failed = 1;
	}
	$FeatnameLookupTag_hash{$feat_name} = $tag;
	$feat_hash{$feat_name}->{'near5'} = $att_line[2];
	$feat_hash{$feat_name}->{'near3'} = $att_line[3];
	$feat_hash{$feat_name}->{'core5'} = $att_line[4];
	$feat_hash{$feat_name}->{'core3'} = $att_line[5];
	print STDERR "$feat_name $feat_hash{$feat_name}->{'near5'} $feat_hash{$feat_name}->{'near3'} $feat_hash{$feat_name}->{'core5'} $feat_hash{$feat_name}->{'core3'} $tag\n" if ($DEBUG);
    }
    close (ATTFILE);
    if ($failed) {
	exit(1);
    }

    return;
}

sub select_max_scores_from_btab { 

    my $blastin = shift (@_);
    my @btab_line = (); # array variable to store split btab lines
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $qtag; # query genome tag
    my $stag; # subject genome tag
    my $score = ""; # BLAST bit score

    open (my $infile, "<", "$blastin") || die ("ERROR: can't open file $blastin\n");

    while (<$infile>)  {
	chomp;
	@btab_line = split(/\t/);
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
	$score = $btab_line[10];
	$qid = $btab_line[0];
	$sid = $btab_line[1];
	if (!defined($qid) || !defined $feat_hash{$qid}) {print STDERR "WARNING!!! $qid is a feature identifier in the btab file ($blastin) but not in the new cluster file ($cluster_file) skipping this entry!\n"; next;}
        if (!defined($sid) || !defined $feat_hash{$sid}) {print STDERR "WARNING!!! $sid is a feature identifier in the btab file ($blastin) but not in the new cluster file ($cluster_file) skipping this entry!\n"; next;}
	$qtag = $FeatnameLookupTag_hash{$qid};
	$stag = $FeatnameLookupTag_hash{$sid};
	if (!defined($qtag) || !defined $TagIndex{$qtag}) {print STDERR "WARNING!!! $qtag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n"; next;}
	if (!defined($stag) || !defined $TagIndex{$stag}) {print STDERR "WARNING!!! $stag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n"; next;}
	print STDERR "$qtag $qid $stag $sid $score\n" if ($DEBUG);
	if (!defined $relationship_hash{$qid}{$sid}) {
	    $relationship_hash{$qid}{$sid} = {};
	    if (!defined $BestTagMatch{$qid}{$stag}) { # this assumes that the best Blast match for a query in a given subject genome appears first in the file
		$BestTagMatch{$qid}->{$stag} = $score;
	    } else {
		if (!defined $SecondBestTagMatch{$qid}{$stag}) { # this assumes that the second best Blast match for a query in a given subject genome appears second in the file
		    $SecondBestTagMatch{$qid}->{$stag} = $score;
		}
		if ($score > $BestTagMatch{$qid}->{$stag}) { # correct Best scores if sort assumption is violated
		    $SecondBestTagMatch{$qid}->{$stag} = $BestTagMatch{$qid}->{$stag};
		    $BestTagMatch{$qid}->{$stag} = $score;
		} elsif ($score > $SecondBestTagMatch{$qid}->{$stag}) {
		    $SecondBestTagMatch{$qid}->{$stag} =  $score;
		}
	    }
	}
	if (!defined $relationship_hash{$sid}{$qid}) {
	    $relationship_hash{$sid}{$qid} = {};
	    if (!defined $BestTagMatch{$sid}{$qtag}) { # unfortunately need to keep this symmetric since blast matches amy not be complete for some reason
		$BestTagMatch{$sid}->{$qtag} = $score;
	    } else {
		if (!defined $SecondBestTagMatch{$sid}{$qtag}) { # can no longer assume highest values cam first since we are storing both ways
		    $SecondBestTagMatch{$sid}->{$qtag} = $score;
		}
		if ($score > $BestTagMatch{$sid}->{$qtag}) { # correct Best scores if sort assumption is violated
		    $SecondBestTagMatch{$sid}->{$qtag} = $BestTagMatch{$sid}->{$qtag};
		    $BestTagMatch{$sid}->{$qtag} = $score;
		} elsif ($score > $SecondBestTagMatch{$sid}->{$qtag}) {
		    $SecondBestTagMatch{$sid}->{$qtag} =  $score;
		}
	    }
	}
    }
    %relationship_hash = (); #clear relationship_hash
    close ($infile);
}

sub select_data_from_btab {

    my $blastin = shift (@_);
    my @btab_line = (); # array variable to store split btab lines
    my $qmatch_length = ""; # stores the query length
    my $smatch_length = ""; # stores the subject length
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $qtag; # query genome tag
    my $stag; # subject genome tag
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $pid = ""; # percent identity
    my $evlu = ""; # e-value
    my $score = ""; # BLAST bit score
    my $qlength = ""; # size of query protein sequence
    my $slength = ""; # size of subject (database match) protein sequence

    open (my $infile, "<", "$blastin") || die ("ERROR: can't open file $blastin\n");

    while (<$infile>)  {
	chomp;
	@btab_line = split(/\t/);
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
	$qbegin = $btab_line[3];
	$qend = $btab_line[4];
	$sbegin = $btab_line[6];
	$send = $btab_line[7];
	$qid = $btab_line[0];
	$sid = $btab_line[1];
	$pid = $btab_line[2];
	$evlu = $btab_line[9];
	$score = $btab_line[10];
	print STDERR "$qid $sid $score\n" if ($DEBUG);
	if (!defined($qid) || !defined $feat_hash{$qid}) {print STDERR "WARNING!!! $qid is a feature identifier in the btab file ($blastin) but not in the new cluster file ($cluster_file) skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($sid) || !defined $feat_hash{$sid}) {print STDERR "WARNING!!! $sid is a feature identifier in the btab file ($blastin) but not in the new cluster file ($cluster_file) skipping this entry!\n" if ($DEBUG); next;}
	$qtag = $FeatnameLookupTag_hash{$qid};
	$stag = $FeatnameLookupTag_hash{$sid};
	if (!defined($qtag) || !defined $TagIndex{$qtag}) {print STDERR "WARNING!!! $qtag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	if (!defined($stag) || !defined $TagIndex{$stag}) {print STDERR "WARNING!!! $stag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	$qlength = $feat_hash{$qid}->{'length'};
	$slength = $feat_hash{$sid}->{'length'};
	$qmatch_length = ((abs($qbegin - $qend) + 1)/$qlength)*100; 
	$smatch_length = ((abs($sbegin - $send) + 1)/$slength)*100;
	if (!defined $BestTagMatch{$qid}{$stag}) {
	    print STDERR "ERROR: $qid and genome $stag had no score the first time through!\n";
	    exit(1);
	}
	if ((0.8 * $BestTagMatch{$qid}->{$stag}) > $score) {
	    print STDERR "Skipping score $score < 80% of max $BestTagMatch{$qid}->{$stag} Query: $qid X Subject $sid = $pid\n" if ($DEBUG);
	    next; #do not keep matches scoring less than 80% of the best score
	}
	if (defined $relationship_hash{$qid} && defined $relationship_hash{$qid}{$sid}) {
	    if ($relationship_hash{$qid}{$sid}->{'score'} >= $score) { # ignore lower scores for the same two proteins
		next;
	    }
	} else {
	    $Qbytaghash{$qid}{$stag}->{$sid} = $relationship_hash{$qid}{$sid} = {}; #have Qbytaghash and relationship_hash reference the same hash
	}
	    # make sure $qbegin is smaller than $qend and $sbegin is smaller than $send
	my $temp_val;
	if ($qbegin > $qend) { #swap
	    $temp_val = $qend;
	    $qend = $qbegin;
	    $qbegin = $temp_val;
	}
	if ($sbegin > $send) { #swap
	    $temp_val = $send;
	    $send = $sbegin;
	    $sbegin = $temp_val;
	}
	$relationship_hash{$qid}{$sid}->{'id'} = $pid;
	$relationship_hash{$qid}{$sid}->{'eval'} = $evlu;
	$relationship_hash{$qid}{$sid}->{'score'} = $score;
	$relationship_hash{$qid}{$sid}->{'min_query'} = $qbegin;
	$relationship_hash{$qid}{$sid}->{'max_query'} = $qend;
	$relationship_hash{$qid}{$sid}->{'min_sub'} = $sbegin;
	$relationship_hash{$qid}{$sid}->{'max_sub'} = $send;
	$relationship_hash{$qid}{$sid}->{'best'} = 0;
	$relationship_hash{$qid}{$sid}->{'bibest'} = 0;
	$relationship_hash{$qid}{$sid}->{'synbest'} = 0;
	$relationship_hash{$qid}{$sid}->{'synbibest'} = 0;
	$relationship_hash{$qid}{$sid}->{'CGN_bibest'} = 0;
	print STDERR "Query: $qid X Subject $sid = $relationship_hash{$qid}{$sid}->{'id'}\n" if ($DEBUG);
	if (!defined $Tagbyfeatnamehash{$qtag}{$qid})  {
	    $Tagbyfeatnamehash{$qtag}{$qid} = 1;
	} elsif ($Tagbyfeatnamehash{$qtag}{$qid} == 0) { #if qid had only been seen as a target (sid) before change to seen as a query (qid)
		$Tagbyfeatnamehash{$qtag}{$qid} = 1;
	}
	if (!defined $Tagbyfeatnamehash{$stag}{$sid})  {
	    $Tagbyfeatnamehash{$stag}{$sid} = 0;
	}
    }
    close ($infile);

    open (OUTMISSING, ">", "missing_blast_results.txt") || die ("ERROR: can't open file missing_blast_results.txt\n");
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	my $qtag = $FeatnameLookupTag_hash{$qid};
	if (!defined $Tagbyfeatnamehash{$qtag}{$qid})  {#output feat_names in the attribute file but missing from the blast results, and cleanup hashes to ignore these proteins
	    print OUTMISSING "$qtag:$qid\n";
	    delete $feat_hash{$qid};
	    delete $FeatnameLookupTag_hash{$qid};
	    delete $TagByPointer{$qid};
	} elsif ($Tagbyfeatnamehash{$qtag}{$qid} == 0) { #if qid had only been seen as a target
	    print OUTMISSING "$qtag:$qid only search result\n";
	}
    }
    close (OUTMISSING);
    print STDERR "check that TagByPointer is consistent\n" if ($DEBUG);
#check that BestMatchTag and SecondBestMatchTag are consistent
#also determine if self-hit exists and create self-hit record if not, no score should be above the self-hit so reset them
    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	my $num_no_selfhits = 0;
	foreach my $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    my $max_score = 0;
	    my $selfhit = 1;
	    my $selfhit_score;
	    if (defined $relationship_hash{$qid} && defined $relationship_hash{$qid}{$qid}) {
		$selfhit_score = $relationship_hash{$qid}{$qid}->{'score'};
		if (!defined $BestTagMatch{$qid}{$qtag}) {
		    die ("ERROR: selfhit was found previously but not here!\n");
		}
	    } else {
		print STDERR "WARNING: no self hit Blast score for $qtag : $qid\n" if ($DEBUG);
		$selfhit = 0;
		$num_no_selfhits++;
	    }
	    if (!defined $BestTagMatch{$qid}) {
		die ("ERROR: no  Blast scores for $qtag : $qid when some were found previously\n");
	    }
	    # Find maximum blast score to use for selfhit score if there is no selfhit score
	    if (!$selfhit) {
		foreach my $stag (@tag_array) { # loop through all genomes
		    if (defined $BestTagMatch{$qid}{$stag} && ($BestTagMatch{$qid}{$stag} > $max_score)) {
			$max_score = $BestTagMatch{$qid}{$stag};
		    }
		}
		$selfhit_score = $max_score;
	    }
	    foreach my $stag (@tag_array) { # loop through all genomes
		print STDERR "$qtag:$qid:$stag - $BestTagMatch{$qid}{$stag}:$SecondBestTagMatch{$qid}{$stag}\n" if ($DEBUG);
		if (defined $SecondBestTagMatch{$qid}{$stag} && !defined $BestTagMatch{$qid}{$stag}) {
		    die ("ERROR: Second Best Blast score for $qtag : $qid to genome $stag exists but not Best!\n");
		}
		if (!defined $BestTagMatch{$qid}{$stag}) {
		    next;
		}
		if ($BestTagMatch{$qid}{$stag} > $selfhit_score) {
		    $BestTagMatch{$qid}{$stag} = $selfhit_score;
		    if (defined $SecondBestTagMatch{$qid}{$stag} && ($SecondBestTagMatch{$qid}{$stag} > $selfhit_score)) {
			$SecondBestTagMatch{$qid}{$stag} = $selfhit_score;
		    }
		}
	    }
	    if (!$selfhit) {
		$BestTagMatch{$qid}{$qtag} = $max_score;
		if ((defined $Qbytaghash{$qid}{$qtag}->{$qid}) || (defined $relationship_hash{$qid}{$qid})) {
		    die ("ERROR: selfhit defined here ($qtag $qid) but not above\n");
		}
		$Qbytaghash{$qid}{$qtag}->{$qid} = $relationship_hash{$qid}{$qid} = {}; #have Qbytaghash and relationship_hash reference the same hash
		$relationship_hash{$qid}{$qid}->{'score'} = $max_score;
		$relationship_hash{$qid}{$qid}->{'id'} = 100;
		$relationship_hash{$qid}{$qid}->{'eval'} = 0;
		$relationship_hash{$qid}{$qid}->{'min_query'} = 1;
		$relationship_hash{$qid}{$qid}->{'max_query'} = $feat_hash{$qid}->{'length'};
		$relationship_hash{$qid}{$qid}->{'min_sub'} = 1;
		$relationship_hash{$qid}{$qid}->{'max_sub'} = $feat_hash{$qid}->{'length'};
		$relationship_hash{$qid}{$qid}->{'best'} = 0;
		$relationship_hash{$qid}{$qid}->{'bibest'} = 0;
		$relationship_hash{$qid}{$qid}->{'synbest'} = 0;
		$relationship_hash{$qid}{$qid}->{'synbibest'} = 0;
		$relationship_hash{$qid}{$qid}->{'CGN_bibest'} = 0;
	    }
	}
	if ($num_no_selfhits > 0) {
	    print STDERR "WARNING: Number of missing Blast selfhits for genome $qtag is $num_no_selfhits\n";
	}
    }

    print STDERR "reduce any scores above self score to self score\n" if ($DEBUG);
#reduce any scores above self score to self score
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	if (!defined $relationship_hash{$qid}) {
	    print STDERR "WARNING: no matches to any other gene for $qid!\n";
	    next;
	}
	if (!defined $relationship_hash{$qid}{$qid}) {
	    die ("ERROR: $qid is missing a selfhit after forcing one!\n");
	}
	foreach my $sid (keys %{ $relationship_hash{$qid} })  { # go through all featnames
	    if ($relationship_hash{$qid}{$sid}->{'score'} > $relationship_hash{$qid}{$qid}->{'score'}) {
		print STDERR "WARNING ($qid $sid $relationship_hash{$qid}{$sid}->{'score'}) > selfhit $relationship_hash{$qid}{$qid}->{'score'}\n" if ($DEBUG);
		$relationship_hash{$qid}{$sid}->{'score'} = $relationship_hash{$qid}{$qid}->{'score'};
	    }
	}
    }
    return;
}

sub calc_BSR  { # subroutine to generate the average (both directions) Blast Score Ratio (BSR) and populate Qbytaghash and relationship_hash with BSR

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    foreach $stag (@tag_array)  {
		if (defined $Qbytaghash{$qid}{$stag}) { # if query protein matches anything in subject genome, lets drill through each match
		    foreach $sid (keys %{ $Qbytaghash{$qid}{$stag} })  { # need to check value of $sid here (see comment in previous subroutine)
			# calculate the Blast Score Ratio (BSR) of each relationship
			if (!defined $relationship_hash{$sid}{$qid}) {
			    $Qbytaghash{$sid}{$qtag}->{$qid} = $relationship_hash{$sid}{$qid} = {}; #have Qbytaghash and relationship_hash reference the same hash
			    #setting these values to make a symmetric entry when the match in one direction was discarded
			    if ($relationship_hash{$qid}{$sid}->{'score'} > $relationship_hash{$sid}{$sid}->{'score'}) { # prevents norm_score being greater than 1
				$relationship_hash{$sid}{$qid}->{'score'} = $relationship_hash{$sid}{$sid}->{'score'};
			    } else {
				$relationship_hash{$sid}{$qid}->{'score'} = $relationship_hash{$qid}{$sid}->{'score'};
			    }
			    $relationship_hash{$sid}{$qid}->{'id'} = $relationship_hash{$qid}{$sid}->{'id'};
			    $relationship_hash{$sid}{$qid}->{'eval'} = $relationship_hash{$qid}{$sid}->{'eval'};
			    $relationship_hash{$sid}{$qid}->{'min_query'} = $relationship_hash{$qid}{$sid}->{'min_sub'};
			    $relationship_hash{$sid}{$qid}->{'max_query'} = $relationship_hash{$qid}{$sid}->{'max_sub'};
			    $relationship_hash{$sid}{$qid}->{'min_sub'} = $relationship_hash{$qid}{$sid}->{'min_query'};
			    $relationship_hash{$sid}{$qid}->{'max_sub'} = $relationship_hash{$qid}{$sid}->{'max_query'};
			    $relationship_hash{$sid}{$qid}->{'best'} = 0;
			    $relationship_hash{$sid}{$qid}->{'bibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'CGN_bibest'} = 0;
			}
			if (!defined $relationship_hash{$qid}{$qid}) {
			    die ("ERROR: somehow we got to calc_BSR without a selfhit for $qid!\n");
			}
			my $norm_score;
			$relationship_hash{$qid}{$sid}->{'norm_score'} = $norm_score = $relationship_hash{$qid}{$sid}->{'score'} / $relationship_hash{$qid}{$qid}->{'score'};
			if (($norm_score > 1) || ($norm_score < 0)) {
			    die ("ERROR: normalized BSR score for $qtag:$qid,$stag:$sid is not between 0-1 ($norm_score)\n");
			}
			if (!defined $relationship_hash{$sid}{$sid}) {
			    die ("ERROR: somehow we got to calc_BSR without a selfhit for $sid!\n");
			}
			$relationship_hash{$sid}{$qid}->{'norm_score'} = $norm_score = $relationship_hash{$sid}{$qid}->{'score'} / $relationship_hash{$sid}{$sid}->{'score'};
			if (($norm_score > 1) || ($norm_score < 0)) {
			    die ("ERROR: normalized BSR score for $stag:$sid;$qtag:$qid is not between 0-1 ($norm_score)\n");
			}
			#BSR (Blast Score Ratio) is mean of two normalized scores
			$Qbytaghash{$qid}{$stag}{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = ($relationship_hash{$qid}{$sid}->{'norm_score'} + $relationship_hash{$sid}{$qid}->{'norm_score'})/2 ; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#alternate choice of BSR to be max of two normalized scores - this did not work but we do use the max of the norm scores in a few places
			#if ($relationship_hash{$qid}{$sid}->{'norm_score'} > $relationship_hash{$sid}{$qid}->{'norm_score'}) {
			#    $Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = $relationship_hash{$qid}{$sid}->{'norm_score'}; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#} else {
			#    $Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = $relationship_hash{$sid}{$qid}->{'norm_score'}; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#}
		    }
		}
	    }
	}
    }
			
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if ((!defined $relationship_hash{$qid}{$sid}->{'BSR'}) || (!defined $relationship_hash{$sid}{$qid}->{'BSR'})) {
		print STDERR "ERROR!!! BSR undefined for $qid $sid in relationship_hash\n";
	    }
	}
    }
}

sub calc_bibest  { # subroutine to find bidrectional best blast matches and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best blast hit per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		my $best_score = 0;
		my $second_score = 0;
		my $found_best = 0;
		my $best_sid;
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'BSR'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'BSR'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$found_best) {
			$best_score = $relationship_hash{$qid}{$sid}->{'BSR'};
			$found_best = 1;
			$best_sid = $sid;
			next;
		    }
		    $second_score = $relationship_hash{$qid}{$sid}->{'BSR'};
		    last; #quit after finding second best blast match
		}
		if ($found_best && ((0.95 * $best_score) > $second_score)) {
		    $relationship_hash{$qid}{$best_sid}->{'best'} = 1;
		}
	    }
	}
    }

    #mark bidirectionally best blast hit per genome per qid if it exists
    foreach $qid (keys %relationship_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'best'} && $relationship_hash{$sid}{$qid}->{'best'}) {
		    $relationship_hash{$qid}{$sid}->{'bibest'}  = 1;
		}
	}
    }
			
}

sub calc_synbibest  { # subroutine to find bidrectional best synteny scores and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best synteny score per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		my $best_score = 0;
		my $second_score = 0;
		my $found_best = 0;
		my $best_sid;
#		my $second_sid;
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'synteny'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'synteny'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$found_best) {
			$best_score = $relationship_hash{$qid}{$sid}->{'synteny'};
			$found_best = 1;
			$best_sid = $sid;
			next;
		    }
		    $second_score = $relationship_hash{$qid}{$sid}->{'synteny'};
#		    $second_sid = $sid;
		    last; #quit after finding second best synteny score
		}
#		if ($found_best) {
#		    $relationship_hash{$qid}{$best_sid}->{'synbest'} = 1;
#		    if ((0.5 * $best_score) < $second_score) {
#			$relationship_hash{$qid}{$second_sid}->{'synbest'} = 1;# allow more than one synbest for realtively good second scores
#		    }
#		}
		if ($found_best && ((0.95 * $best_score) >= $second_score)) {
		    $relationship_hash{$qid}{$best_sid}->{'synbest'} = 1;
		}
	    }
	}
    }

    #mark bidirectionally best synteny scores per genome per qid if it exists
    foreach $qid (keys %relationship_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'synbest'} && $relationship_hash{$sid}{$qid}->{'synbest'}) {
		    $relationship_hash{$qid}{$sid}->{'synbibest'}  = 1;
		}
	}
    }
			
}

sub calc_synteny  { # subroutine to calculate synteny scores and populate Qbytaghash and relationship_hash with them

   #call synteny score for each relationship
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if ($qid eq $sid) {
		next; #skip matches to same sequence
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    if (!defined $qtag) {
		print STDERR "$qid has no defined tag in FeatnameLookupTag_hash\n";
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if (!defined $stag) {
		print STDERR "$sid has no defined tag in FeatnameLookupTag_hash\n";
	    }
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    print STDERR "$qid $sid $relationship_hash{$qid}{$sid}->{'id'} $relationship_hash{$qid}{$sid}->{'eval'} $relationship_hash{$qid}{$sid}->{'score'} $relationship_hash{$qid}{$sid}->{'best'} $relationship_hash{$qid}{$sid}->{'bibest'}\n" if($DEBUG);
	    $relationship_hash{$qid}{$sid}->{'synteny'} = &synteny($qtag, $qid, $stag, $sid);
	}
    }
			
}

sub merge_clusters {#&merge_clusters(qid, sid) is used to test if two clusters can be merged and if so to merge them
    my ($qid, $sid) = @_;
    my $merged = [];
    my $failed_merge = 0;
    
    #assumes the cluster arrays are sorted by genome tag index
    my $index = 0; #need to keep track of where we are in the subject array
    foreach my $cur_qid_tag ( @{ $clusters{$qid} } ) {
	while (($index <= $#{ $clusters{$sid} }) && (($clusters{$sid}->[$index])->{'tag'} < $cur_qid_tag->{'tag'})) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	if ($index <= $#{ $clusters{$sid} }) {
	    if (($clusters{$sid}->[$index])->{'tag'} == $cur_qid_tag->{'tag'}) {
		if (($clusters{$sid}->[$index])->{'id'} eq $cur_qid_tag->{'id'}) {
		    die ("ERROR: two different clusters should never contain the same protein\n");
		} else {
		    $failed_merge = 1; #cannot have two different proteins from the same genome
		    last;
		}
	    }
	}
	push(@{ $merged }, $cur_qid_tag);
	print STDERR Dumper($merged) if ($DEBUG);
    }
    if (!$failed_merge) {
	while ($index <= $#{ $clusters{$sid} }) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	foreach my $cur_id_tag ( @{ $clusters{$qid} }, @{ $clusters{$sid} } ) {
	    $clusters{$cur_id_tag->{'id'}} = $merged;
	}
    } else {
	print STDERR "failed merge\n" if ($DEBUG);
    }
    
}

sub calc_clusters { # greedily compute clusters by starting with largest relationship score to merge existing clusters (starting with every protein as a singleton) and constrainging a cluster to only contain one protein from each genome

    # start by creating a sorted array of relationships from the relationship_hash
    my @relationship_array = ();
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $qTagIndex = $TagIndex{$qtag};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    my $sTagIndex = $TagIndex{$stag};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if ((!$relationship_hash{$qid}{$sid}->{'bibest'}) && (!$relationship_hash{$qid}{$sid}->{'synbibest'})) {
		next; # ignore matches which are not bidirectionally best matches or bidirectionally best synteny matches
	    }
	    if ($relationship_hash{$qid}{$sid}->{'synteny'} < 2) {
		next; # ignore matches which do not meet this minimum
	    }
	    my $score = $relationship_hash{$qid}{$sid}->{'synteny'};
	    if (defined $relationship_hash{$sid}{$qid}) {
		if ($qid gt $sid) {
		    next; # only enter qid,sid pairs once not symmetrically since they indicate the same cluster join
		}
		if ($relationship_hash{$sid}{$qid}->{'synteny'} > $score) {# use the larger synteny score
		    $score = $relationship_hash{$sid}{$qid}->{'synteny'};
		}
	    }
	    my $hash_ref = {};
	    $hash_ref->{'query'} = $qid;
	    $hash_ref->{'subject'} = $sid;
	    $hash_ref->{'score'} = $score;
	    push (@relationship_array, $hash_ref);
	}
    }
    @relationship_array = sort {$b->{'score'} <=> $a->{'score'}} ( @relationship_array );
    print STDERR Dumper(\@relationship_array) if ($DEBUG);

    foreach my $qid (keys %relationship_hash)  { # go through all featnames and initialize clusters as singletons
	if (!defined $FeatnameLookupTag_hash{$qid}) {
	    die("ERROR!!! $qid in relationship_hash but not in FeatnameLookupTag_hash\n");
	}
	$clusters{$qid} = [{ 'id' => $qid, 'tag' => $TagIndex{$FeatnameLookupTag_hash{$qid}} }];
    }
    print STDERR Dumper(\%clusters) if ($DEBUG);
    foreach my $match ( @relationship_array ) {
	if ($clusters{$match->{'query'}} == $clusters{$match->{'subject'}}) {
	    next; #already in the same cluster
	}
	print STDERR "merge ($match->{'query'}, $match->{'subject'}): $match->{'score'}\n" if ($DEBUG);
	&merge_clusters($match->{'query'}, $match->{'subject'});
    }
    print STDERR "Done clustering\n" if ($DEBUG);
    print STDERR Dumper(\%clusters) if ($DEBUG);
			
}

sub calc_cluster_numbers { # order the cluster numbers arbitrarily

############ Assign cluster numbers to clusters in the same order they will be output in the match tables and compute cluster_size ###############
    my $cluster_num = $start_cluster;
    
    print STDERR "Calculating cluster numbers ...\n";
    foreach my $query_featname (keys %feat_hash)  { # go through all featnames
	foreach my $cluster_member ( @{ $clusters{$query_featname} } ) {
	    if (defined $cluster_number{$cluster_member->{'id'}}) { # we have seen this cluster before so skip and reset cluster number
		$cluster_num--;
		last;
	    } else {
		$cluster_number{$cluster_member->{'id'}} = $cluster_num;
		$cluster_for_number[$cluster_num] = $clusters{$query_featname};
		$cluster_size{$cluster_num}++;
	    }
	}
	$cluster_num++;
    }
    return;
}

sub calc_cluster_index {

############ Output a two column tab delimited file with feat_name and cluster number ###############
    
    open (my $indexfile, ">", $index_file) || die ("ERROR: can't open file $index_file\n");
    foreach my $query_featname (keys %cluster_number)  { # go through all featnames
	print $indexfile "$query_featname\t$cluster_number{$query_featname}\n";
    }
    close ($indexfile);
    return;
}

sub calc_cluster_medoids { # determine the medoid of each cluster based on BSR scores

    
    open (my $medoidfile, ">", $medoid_file) || die ("ERROR: can't open file $medoid_file\n");
    foreach my $cluster_id (sort {$a <=> $b} keys %cluster_size)  {
	my $max_BSR_sum = -1;
	my $medoid = undef;
	print STDERR "cluster $cluster_id\n" if ($DEBUG);
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $qid = $cluster_member->{'id'};
	    my $BSR_sum = 0;
	    print STDERR "$qid\n" if ($DEBUG);
	    foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
		print STDERR "\t$sid\t" if ($DEBUG);
		if ($cluster_id != $cluster_number{$sid}) {
		    print STDERR "in $cluster_number{$sid} $relationship_hash{$qid}{$sid}->{'BSR'}\n" if ($DEBUG);
		    next; #skip matches outside of the cluster
		}
		print STDERR "$relationship_hash{$qid}{$sid}->{'BSR'}\n" if ($DEBUG);
		$BSR_sum += $relationship_hash{$qid}{$sid}->{'BSR'};
	    }
	    print STDERR "$BSR_sum\n" if ($DEBUG);
	    if ($BSR_sum > $max_BSR_sum) {
		$max_BSR_sum = $BSR_sum;
		$medoid = $qid;
	    }
	}
	print STDERR "$max_BSR_sum\n" if ($DEBUG);
	if (!defined $medoid) {
	    die "medoid could not be determined for cluster $cluster_id !\n";
	}
	$medoids[$cluster_id] = $medoid;
	print $medoidfile ">centroid_$cluster_id $medoid $feat_hash{$medoid}->{'header'}\n";
	my $seq_len = $feat_hash{$medoid}->{'length'};
	my $sequence = $feat_hash{$medoid}->{'sequence'};
	for ( my $pos = 0 ; $pos < $seq_len ; $pos += 60 ) {
	    print $medoidfile substr($sequence, $pos, 60), "\n";
	}

    }
    close ($medoidfile);

}

sub synteny {#&synteny($genome_tag, $query_featname, $subject_tag, $subject_featname)

    my ($query_tag, $query_featname, $subject_tag, $subject_featname) = @_;
    my $total_score = 0;
    
    print STDERR "query $query_featname   subject $subject_featname\n" if ($DEBUG);
    my @q_near5 = split(/:/, $feat_hash{$query_featname}->{'near5'});
    my @q_near3 = split(/:/, $feat_hash{$query_featname}->{'near3'});
    my @q_core5 = split(/:/, $feat_hash{$query_featname}->{'core5'});
    my @q_core3 = split(/:/, $feat_hash{$query_featname}->{'core3'});
    my @s_near5 = split(/:/, $feat_hash{$subject_featname}->{'near5'});
    my @s_near3 = split(/:/, $feat_hash{$subject_featname}->{'near3'});
    my @s_core5 = split(/:/, $feat_hash{$subject_featname}->{'core5'});
    my @s_core3 = split(/:/, $feat_hash{$subject_featname}->{'core3'});
    foreach my $q_near (@q_near5) {
	(my $q_clus, my $q_inv) = split(/_/, $q_near);
	foreach my $s_near (@s_near5) {
	    (my $s_clus, my $s_inv) = split(/_/, $s_near);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
		if ($q_inv == $s_inv) {
		    $total_score += 0.5;
		}
	    }
	}
    }
    foreach my $q_near (@q_near3) {
	(my $q_clus, my $q_inv) = split(/_/, $q_near);
	foreach my $s_near (@s_near3) {
	    (my $s_clus, my $s_inv) = split(/_/, $s_near);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
		if ($q_inv == $s_inv) {
		    $total_score += 0.5;
		}
	    }
	}
    }
    foreach my $q_near (@q_near5) {
	(my $q_clus, my $q_inv) = split(/_/, $q_near);
	foreach my $s_near (@s_near3) {
	    (my $s_clus, my $s_inv) = split(/_/, $s_near);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
	    }
	}
    }
    foreach my $q_near (@q_near3) {
	(my $q_clus, my $q_inv) = split(/_/, $q_near);
	foreach my $s_near (@s_near5) {
	    (my $s_clus, my $s_inv) = split(/_/, $s_near);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
	    }
	}
    }
    foreach my $q_core (@q_core5) {
	(my $q_clus, my $q_inv, my $q_dist) = split(/_/, $q_core);
	foreach my $s_core (@s_core5) {
	    (my $s_clus, my $s_inv, my $s_dist) = split(/_/, $s_core);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
		if ($q_inv == $s_inv) {
		    $total_score += 0.5;
		    my $score = abs($q_dist - $s_dist);
		    if ($score > 10) {$score = 10;}
		    $score = (10 - $score) / 10;
		    $total_score += $score;
		}
	    }
	}
    }
    foreach my $q_core (@q_core3) {
	(my $q_clus, my $q_inv, my $q_dist) = split(/_/, $q_core);
	foreach my $s_core (@s_core3) {
	    (my $s_clus, my $s_inv, my $s_dist) = split(/_/, $s_core);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
		if ($q_inv == $s_inv) {
		    $total_score += 0.5;
		    my $score = abs($q_dist - $s_dist);
		    if ($score > 10) {$score = 10;}
		    $score = (10 - $score) / 10;
		    $total_score += $score;
		}
	    }
	}
    }
    foreach my $q_core (@q_core5) {
	(my $q_clus, my $q_inv, my $q_dist) = split(/_/, $q_core);
	foreach my $s_core (@s_core3) {
	    (my $s_clus, my $s_inv, my $s_dist) = split(/_/, $s_core);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
	    }
	}
    }
    foreach my $q_core (@q_core3) {
	(my $q_clus, my $q_inv, my $q_dist) = split(/_/, $q_core);
	foreach my $s_core (@s_core5) {
	    (my $s_clus, my $s_inv, my $s_dist) = split(/_/, $s_core);
	    if ($q_clus == $s_clus) {
		$total_score += 0.5;
	    }
	}
    }
    print STDERR "     total_score $total_score\n" if ($DEBUG);
    return($total_score);

}

sub write_files  {

    my ($choice, $query_column, $subject_column, $query_featname, $subject_tag, $subject_featname, $cluster_number) = @_;

    if ($subject_column == 1) {
	    print OUTVENN "$cluster_number\t";
    }
    if ($choice == "1") {
	if ($query_column == 1) {
	    print OUTVENN "$query_featname";
	} else {
	    print OUTVENN "\t$query_featname";
	}
    } elsif ($choice == "2") {
        if ($subject_column != 1) {
	    print OUTVENN "\t";
	}
	print OUTVENN "$subject_featname";
    } elsif ($choice == "3") {
	if ($subject_column != 1) {
	    print OUTVENN "\t"; 
	}
	print OUTVENN "----------";  # if no match print the default delimiter ----------
    } elsif ($choice == "4")  {
	print OUTVENN "\n";
    }
}

sub make_tables {

############ Generate a match table of bidirectional best hits ###############
    
    open (OUTVENN, ">", $match_file) || die ("ERROR: can't open file $match_file\n");
    foreach my $cluster_id (sort {$a <=> $b} keys %cluster_size)  { # output in previously computed order of the clusters
	my $cluster_member_first = ${ $cluster_for_number[$cluster_id] }[0];
	my $genome_tag_index = $cluster_member_first->{'tag'};
	my $query_featname = $cluster_member_first->{'id'};
	my $i = $genome_tag_index + 1; #used to know how many columns to skip for some tables
	my $subject_tag_index;
	my $subject_featname = "";
	for ($subject_tag_index = 0; $subject_tag_index < $genome_tag_index; $subject_tag_index++) {
	    &write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	}
	# We are at $genome_tag now so we are looking at self comparisons, generate output and move on
	&write_files("1", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);#$subject_featname has not been defined yet here!!
	$subject_tag_index++; #go past $genome_tag
	# Now iterate through sorted (by subject_tag_index) cluster array to print out members of cluster
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $cluster_member_tag_index = $cluster_member->{'tag'};
	    if ($cluster_member_tag_index < $subject_tag_index) {
		if ($cluster_member->{'id'} ne $query_featname) {
		    die ("ERROR!!! cluster for $query_featname should have been output sooner!\n");
		}
		next; #skip $genome_tag member we already output for
	    }
	    while ($cluster_member_tag_index > $subject_tag_index) {
		#output blanks for this genome
		&write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
		$subject_tag_index++; #go to next genome
	    }
	    $subject_featname = $cluster_member->{'id'};
	    &write_files("2", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	    $subject_tag_index++; #go to next genome
	}
	while ($subject_tag_index <= $#tag_array) {
	    &write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	    $subject_tag_index++; #go to next genome
	}
	&write_files("4", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id); #add newlines to files
    }
    close (OUTVENN);
}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog:

           A program to cluster new clusters using a PGG based annotation giving nearest neighbors plus Blast results

Copyright (C) 2019  The J. Craig Venter Institute (JCVI).  All rights reserved

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

Citation:  Derrick E. Fouts, Lauren Brinkac, Erin Beck, Jason Inman, and Granger Sutton 
           (2012) "PanOCT: Automated Clustering of Orthologs using Conserved Gene 
           Neighborhood for Pan-Genomic Analysis of Bacterial Strains and Closely Related 
           Species" Nucleic Acids Res. 40(22):e172.

  Usage: $prog <options>
Example: $prog -i index_file -M match_file -m medoid_file -c cluster_neighbors_file -g genome_list_file -s gene_seqs_file -n start_cluster_number
Version: $version
Options:
     -h: print this help page
     -c: file with nearest neighbors information for each gene [REQUIRED]
     -g: file containing unique genome identifier tags [REQUIRED]
     -s: file containing multifasta for gene sequences used to output medoids [REQUIRED]
     -n: the number of the first additional cluster to be output (next number after current set of clusters) [REQUIRED]
     -i: file name for outputting the index of gene to cluster [REQUIRED]
     -m: file name for outputting the medoids of the new clusters [REQUIRED]
     -M: file name for outputting the new clusters [REQUIRED]
     -B: directory name for where blast executables are located - default is not to use a directory
     -L: directory name for where blast libraries are located - default is not to use a directory
     -D: no argument, DEBUG MODE (DEFAULT = off)

 Authors: Granger Sutton, Ph.D.
 Date: April 22, 2019; last revised April 22, 2019

_EOB_
    exit;
}

########################################  M A I N  #####################################################
print STDERR "Fetching genome tags from $gtag_file\n";
&get_tags;
print STDERR "Gathering sequence information from $seq_file\n";
&get_seqs;
print STDERR "Gathering information from new cluster file: $cluster_file\n";
&get_gene_att;
print STDERR "blasting all versus all of $seq_file";
`mkdir M_BLAST_TMP`;
`cp $seq_file M_BLAST_TMP/temp_fasta.ftmp`;
my $makeblastdb = $blast_directory . "makeblastdb";
`$makeblastdb -in M_BLAST_TMP/temp_fasta.ftmp -dbtype nucl -out M_BLAST_TMP/temp_fasta.ftmp`;
my $blastn = $blast_directory . "blastn";
`$blastn -query $seq_file -db M_BLAST_TMP/temp_fasta.ftmp -out M_BLAST_TMP/temp_results.ftmp -task blastn -evalue 0.000001 -outfmt \"6 qseqid sseqid pident qstart qend qlen sstart send slen evalue bitscore\"`;
print STDERR "Getting data from blast .btab file: M_BLAST_TMP/temp_results.ftmp\n";
&select_max_scores_from_btab("M_BLAST_TMP/temp_results.ftmp");
&select_data_from_btab("M_BLAST_TMP/temp_results.ftmp");
`rm -r M_BLAST_TMP`;
print STDERR "Calculating the Blast Score Ratio (BSR) ...\n";
&calc_BSR;
print STDERR "Calculating reciprocal best hits (RBHs) ...\n";
&calc_bibest;
print STDERR "Calculating conserved gene neighborhood (CGN) scores ...\n";
&calc_synteny;
print STDERR "Calculating recipricol best CGN scores ...\n";
&calc_synbibest;
print STDERR "Calculating clusters ...\n";
&calc_clusters;
&calc_cluster_numbers;
print STDERR "Calculating cluster medoids ...\n";
&calc_cluster_medoids;
print STDERR "Calculating cluster index ...\n";
&calc_cluster_index;
print STDERR "generating match table!\n";
&make_tables;
exit(0);
