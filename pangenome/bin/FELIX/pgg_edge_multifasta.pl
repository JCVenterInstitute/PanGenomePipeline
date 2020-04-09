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
print STDERR "$commandline\n";
my $prog = $0;
$prog =~ s/.*\///;

use Cwd;
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
getopts ('j:I:RSALlDhb:B:m:p:P:a:g:t:M:s:T:Vk:FfC:Q:N:X');
our ($opt_j, $opt_I,$opt_S,$opt_A,$opt_L,$opt_l,$opt_D,$opt_h,$opt_b,$opt_m,$opt_p,$opt_P,$opt_a,$opt_g,$opt_t,$opt_M,$opt_r,$opt_R,$opt_B,$opt_s,$opt_T,$opt_V,$opt_k,$opt_F,$opt_f,$opt_C,$opt_Q,$opt_N,$opt_X);

## use boolean logic:  TRUE = 1, FALSE = 0

my $qsub_job_num = 0;
my $bin_directory = "/usr/local/projdata/8520/projects/PANGENOME/pangenome_bin/";
my $cwd = getcwd;
my $muscle_path = "/usr/local/bin/muscle";
my $rscript_path = "/usr/local/bin/Rscript";
my $version = "ver1.0";
my $project = "8520";
my $Gapped_Context = 100; #number of basepairs around a cluster/edge to check for gaps so as to ignore problematic assembly regions
my $keep_divergent_alignments = "";
my $basedir;
my $multifastadir;
my $use_multifasta = 0;
my $write_multifasta = 0;
my $no_stats = 0;
my $matchtable_file;
my $att_file;
my $pgg_file;
my $DEBUG;
my $genomes_file_name;
my @genome_array = ();
my $genome_number;
my $max_grid_jobs = 50;
my $compute_all = 0;
my $align_all = 0;
my $align_new = 0;
my $suppress = 0;
my $target_id = "";
my $ignore_id = "";
my $ignore_index = -1;
my $remake_files;
my $medoids_path;
my $single_cores;
my $topology_file;
my $strip_version = 0;
my $qsub_queue = "himem";
if ($opt_j) { # should really check that this a positive integer
    if (($opt_j =~ /^\d+$/) && ($opt_j > 0)) {
	$max_grid_jobs = $opt_j;
    } else {
	print STDERR "Error $opt_j for maximum number of grid jobs is not a positive integer\n";
	&option_help;
    }
}
if ($opt_P) {$project = $opt_P;} else {$project = "8520";}
if ($opt_Q) {$qsub_queue = $opt_Q;} else {$qsub_queue = "himem";}
if ($opt_M) {$medoids_path = $opt_M;} else {$medoids_path = "";}
if ($opt_s) {$single_cores = $opt_s;} else {$single_cores = "";}
if ($opt_X) {$no_stats = 1;} else {$no_stats = 0;}
if ($opt_F) {$use_multifasta = 1;} else {$use_multifasta = 0;}
if ($opt_f) {$write_multifasta = 1;} else {$write_multifasta = 0;}
if ($opt_R) {$remake_files = 1;} else {$remake_files = 0;}
if ($opt_V) {$strip_version = 1;} else {$strip_version = 0;}
if ($opt_A) {$compute_all = 1;} else {$compute_all = 0;} # flag to compute statistics for all
if ($opt_S) {$suppress = 1;} else {$suppress = 0;} # flag to suppress outputting multifasta files
if ($opt_L) {
    if ($suppress) {
	print STDERR "Error cannot have both -S (suppress multifasta output) and -L (generate multiple sequence alignments)\n";
	&option_help;
    }
    $align_all = 1;
} else {$align_all = 0;} # flag to compute multiple sequence alignment for all
if ($opt_l) {
    if ($align_all) {
	print STDERR "Error cannot have both -l (compare to multiple sequence alignments) and -L (generate multiple sequence alignments)\n";
	&option_help;
    }
    $align_new = 1;
} else {$align_new = 0;} # flag to compare to multiple sequence alignments
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_t) {$target_id = $opt_t;} #set target genome id for statistics
if ($opt_I) {$ignore_id = $opt_I;} #set a genome id to ignore for statistics (used when the target genome is one of the PGG genomes)
if ($opt_k) {
    $keep_divergent_alignments = $opt_k;
    if ((-e $keep_divergent_alignments) && !(-d $keep_divergent_alignments)) {
	print STDERR "Error with -k $keep_divergent_alignments - file exists but is not a directory\n";
	&option_help;
    } elsif (!(-e $keep_divergent_alignments)) {
	mkdir($keep_divergent_alignments) or die "Could not create directory $keep_divergent_alignments\n";
    }
}
if ($opt_N) {
    $bin_directory = $opt_N;
    if ((-e $bin_directory) && !(-d $bin_directory)) {
	print STDERR "Error with -N $bin_directory - file exists but is not a directory\n";
	&option_help;
    } elsif (!(-e $bin_directory)) {
	print STDERR "Error with -N $bin_directory - does not exist\n";
	&option_help;
    }
}
if ($opt_b) {
    $multifastadir = $opt_b;
    if ((-e $multifastadir) && !(-d $multifastadir)) {
	print STDERR "Error with -b $multifastadir - file exists but is not a directory\n";
	&option_help;
    } elsif (!(-e $multifastadir)) {
	mkdir($multifastadir) or die "Could not create directory $multifastadir\n";
    }
} else { $multifastadir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if ($opt_B) {
    $basedir = $opt_B;
    if ((-e $basedir) && !(-d $basedir)) {
	print STDERR "Error with -B $opt_B - file exists but is not a directory\n";
	&option_help;
    } elsif (!(-e $basedir)) {
	mkdir($basedir) or die "Could not create directory $basedir\n";
    }
} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if (($opt_p) && (-s "$opt_p")) {$pgg_file = $opt_p;} else { print STDERR "Error with -p $opt_p\n"; &option_help; } # if no value for option p (pan-genome graph input file), quit with help menu
if (($opt_m) && (-s "$opt_m")) {$matchtable_file = $opt_m;} else { print STDERR "Error with -m $opt_m\n"; &option_help; } # if no value for option m (matchtable input file), quit with help menu
if (($opt_a) && (-s "$opt_a")) {$att_file = $opt_a;} else { print STDERR "Error with -a $opt_a\n"; &option_help; } # if no value for option a (attribute input file), quit with help menu
if (($opt_g) && (-s "$opt_g")) {$genomes_file_name = $opt_g;} else { print STDERR "Error with -g $opt_g\n"; &option_help; } # if no value for option g (genome tags and contig file names input file), quit with help menu
if (($opt_T) && (-s "$opt_T")) {$topology_file = $opt_T;} else { print STDERR "Error with -T $opt_T\n"; &option_help; } # if no value for option T (topology input file), quit with help menu
if ($opt_C) {
    if (-x "$opt_C") {
	$muscle_path = $opt_C;
    } else {
	print STDERR "Error with -C $opt_C\n"; &option_help;
    }
} else {
    $muscle_path = "/usr/local/bin/muscle";
}
if ($opt_r) {
    if (-x "$opt_r") {
	$rscript_path = $opt_r;
    } else {
	print STDERR "Error with -r $opt_r\n"; &option_help;
    }
} else {
    $rscript_path = "/usr/local/bin/Rscript";
}

my $num_size_one_clus = 0;
my $num_shared_clus = 0;
my $num_core_clus = 0;
my $num_reduced_clus = 0;
my $num_size_one_edge = 0;
my $num_shared_edge = 0;
my $num_core_edge = 0;
my $num_reduced_edge = 0;
my $cpu_name = "$target_id" . "_pem_cpu_separate_stats";
my %edge_hash = ();            # Key1 = edge ID Key2 = struct members with their values (5p,3p,gtag, contig)
my %single_copy_core = ();     # key = cluster number, value is defined if single copy core otherwise undefiend
my %is_circular = ();          # key1 = genome ID, key2 = contig_name, value = 1 if circular 0 otherwise
my %feat_hash = ();            # Key1 = feat_name Key2 = struct members with their values (5p,3p,anno,gtag, contig)
my %cluster_to_feat_hash = (); # Key1 = genome tag Key2 = cluster_id Value = feat_name
my %genseq_hash = ();          # Key1 = genome tag Key2 = contig_name Value = contig sequence
my %genseq_len = ();           # Key1 = genome tag Key2 = contig_name Value = length of contig sequence
my %total_clus = ();           # key = genome ID, value = total number of nonsingleton clusters
my %total_edge = ();           # key = genome ID, value = total number of nonsingleton edges
my $total_clus_pgg = 0;        # number of nonsingleton clusters
my $total_edge_pgg = 0;        # number of nonsingleton edges - only count one direction
my $total_clus_alle_pgg = 0;   # number of alleles in nonsingleton clusters - only count one direction
my $total_edge_alle_pgg = 0;   # number of alleles in nonsingleton edges
my %distant_clus_alle = ();    # key = genome ID, value = number of distant cluster alleles
my %distant_edge_alle = ();    # key = genome ID, value = number of distant edge alleles
my %column_clus_alle = ();     # key = genome ID, value = number of distant (by unique column characters) cluster alleles
my %column_edge_alle = ();     # key = genome ID, value = number of distant (by unique column characters) edge alleles
my %uniq_clus = ();            # key = genome ID, value = number of unique/singleton clusters
my %uniq_edge = ();            # key = genome ID, value = number of unique/singleton edges
my %short_clus = ();           # key = genome ID, value = number of short clusters
my %short_edge = ();           # key = genome ID, value = number of short edges
my %long_clus = ();            # key = genome ID, value = number of long clusters
my %long_edge = ();            # key = genome ID, value = number of long edges
my %very_short_clus = ();      # key = genome ID, value = number of very short clusters
my %very_short_edge = ();      # key = genome ID, value = number of very short edges
my %very_long_clus = ();       # key = genome ID, value = number of very long clusters
my %very_long_edge = ();       # key = genome ID, value = number of very long edges
my %gapped_clus = ();          # key = genome ID, value = number of gapped (really too many ambiguous characters) clusters
my %gapped_edge = ();          # key = genome ID, value = number of gapped (really too many ambiguous characters) edges
my %frameshift = ();           # key = genome ID, value = number of probable frameshifts
my %miss_sing_core = ();       # key = genome ID, value = number of clusters missing  for single copy core clusters
my %miss_sing_edge = ();       # key = genome ID, value = number of edges missing  for edges between single copy core clusters
my %missing_75c = ();          # key = genome ID, value = number of clusters missing  for clusters in 75-100% of genomes
my %missing_75e = ();          # key = genome ID, value = number of edges missing  for edges in 75-100% of genomes
my %uniq_clus_alle_75_100 = ();# key = genome ID, value = number of unique alleles for clusters in 75-100% of genomes
my %uniq_clus_alle_25_75 = (); # key = genome ID, value = number of unique alleles for clusters in 25-75% of genomes
my %uniq_clus_alle_0_25 = ();  # key = genome ID, value = number of unique alleles for clusters in 0-25% of genomes
my %uniq_edge_alle_75_100 = ();# key = genome ID, value = number of unique alleles for edges in 75-100% of genomes
my %uniq_edge_alle_25_75 = (); # key = genome ID, value = number of unique alleles for edges in 25-75% of genomes
my %uniq_edge_alle_0_25 = ();  # key = genome ID, value = number of unique alleles for edges in 0-25% of genomes
my @renumber = ();             # maps old cluster numbers to new cluster numbers
my @mf_files = ();             # a list of multifasta files which are created can can be aligned (not zero length and not singleton)
 
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
	if ($ignore_id eq $tag) {
	    next; #skip over the genome to be ignored
	}
	if (($tag eq "") || ($asmbl_id eq "") || ($type eq "")) {
	    die ("ERROR: genome id, assembly id/contig id, and type  must not be empty/null in the contig topology file $topology_file.\nLine:\n$_\n");
	}
	if ($strip_version) {
	    $asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	if ($type eq "circular") {
	    $is_circular{$tag}->{$asmbl_id} = 1;
	} elsif ($type eq "linear") {
	    $is_circular{$tag}->{$asmbl_id} = 0;
	} else {
	    die ("ERROR: type $type must be either circular or linear in the  contig topology file $topology_file.\nLine:\n$_\n");
	}
    }
    close (CIRCFILE);
    return;
}

sub get_genomes {  # obtain list of genomes - must be in the same order as the matchtable columns - and the multifasta contigs file for the genomes
   
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: cannot open file $genomes_file_name\n");
    print  STDERR "Order of genomes in $genomes_file_name with array index\n" if ($DEBUG);
    my $target_found = 0;
    while (my $line1 = <$infile>)  {
	chomp $line1;
	(my $name, my $contig_file) = split(/\t/, $line1);  # split the scalar $line on tab

	if ($ignore_id eq $name) {
	    $ignore_index = $genome_number;
	    print  STDERR "Ignoring genome $ignore_id with index $ignore_index\n";
	    next; #skip over the genome to be ignored
	}
	if (defined $genseq_hash{$name})  {
	    die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
	} else  {
	    push (@genome_array, $name); # populate the genome_array in the order of the genome file
	    print  STDERR "$name\t$genome_number\n" if ($DEBUG);
	    $genome_number++;
	    if ($target_id ne "") {
		if ($target_id eq $name) {
		    $target_found = 1;
		}
	    }
	}
	if (!$use_multifasta || ($target_id eq $name)) { #only read in the target genome - will use multiple fasta files for edges and clusters instead
	    my $contigfile;
	    unless (open ($contigfile, "<", $contig_file) )  {
		die ("ERROR: cannot open file $contig_file.\n");
	    }
	    my ($save_input_separator) = $/;
	    $/="\n>";
	    while (my $line2 = <$contigfile>) {
		(my $title, my $sequence) = split(/\n/, $line2, 2); # split the header line and sequence (very cool)
		my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
		my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
		$id =~ s/>\s*//; # remove leading > and spaces
		if ($strip_version) {
		    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
		}
		$sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
		#print STDERR "genome: $name contig: $id\n";
		if (!defined $is_circular{$name}) {
		    die ("ERROR: $name is a genome in the genomes fasta list file for $genomes_file_name but not in the contig topology file $topology_file!\n");
		}
		if (!defined $is_circular{$name}->{$id}) {
		    die ("ERROR: $id is a contig in the genome fasta file for $name genome for $genomes_file_name but not in the contig topology file $topology_file!\n");
		}
		$genseq_hash{$name}->{$id} = $sequence;
		$genseq_len{$name}->{$id} = length($sequence);
		$title = ""; # clear the title for the next contig
		$sequence = ""; #clear out the sequence for the next contig
	    }
	    $/ = $save_input_separator; # restore the input separator
	    close ($contigfile);
	} else {
		$genseq_hash{$name}->{"Placeholder"} = "EMPTY"; #this is just here to allow checking for duplicate genome names
	}
    }
    if (($target_id ne "") && !$target_found) {
	die ("ERROR: Did not find target genome: $target_id in genome list file: $genomes_file_name.\n");
    }
    close($infile);
    print  STDERR "$genome_number genomes\n\n";
}

sub get_genome_names {  # obtain list of genome names - must be in the same order as the matchtable columns
   
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: cannot open file $genomes_file_name\n");
    print  STDERR "Order of genomes in $genomes_file_name with array index\n" if ($DEBUG);
    my $target_found = 0;
    while (my $line1 = <$infile>)  {
	chomp $line1;
	(my $name, my $contig_file) = split(/\t/, $line1);  # split the scalar $line on tab

	if ($ignore_id eq $name) {
	    $ignore_index = $genome_number;
	    print  STDERR "Ignoring genome $ignore_id with index $ignore_index\n";
	    next; #skip over the genome to be ignored
	}
	if (defined $genseq_hash{$name})  {
	    die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
	} else  {
	    $genseq_hash{$name}->{"Placeholder"} = "EMPTY"; #this is just here to allow checking for duplicate genome names
	    push (@genome_array, $name); # populate the genome_array in the order of the genome file
	    print  STDERR "$name\t$genome_number\n" if ($DEBUG);
	    $genome_number++;
	    if ($target_id ne "") {
		if ($target_id eq $name) {
		    $target_found = 1;
		}
	    }
	}
   }
    if (($target_id ne "") && !$target_found) {
	die ("ERROR: Did not find target genome: $target_id in genome list file: $genomes_file_name.\n");
    }
    close($infile);
    print  STDERR "$genome_number genomes\n\n";
}

sub output_multifasta {  # obtain list of genomes - must be in the same order as the matchtable columns - and the multifasta contigs file for the genomes - then write multifasta for clusters and edges
   
    $genome_number = 0;     # total number of genomes to be processed
    if ($target_id ne "") {
	die ("ERROR: The -f option for writing a multifasta file should not be used with the -t option for specifyig a target genome!\n");
    }
    if ($ignore_id ne "") {
	die ("ERROR: The -f option for writing a multifasta file should not be used with the -I option for specifyig a genome to ignore!\n");
    }

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: cannot open file $genomes_file_name\n");
    print  STDERR "Order of genomes in $genomes_file_name with array index\n" if ($DEBUG);
    while (my $line1 = <$infile>)  {
	chomp $line1;
	(my $name, my $contig_file) = split(/\t/, $line1);  # split the scalar $line on tab

	if (defined $genseq_hash{$name})  {
	    die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
	} else  {
	    push (@genome_array, $name); # populate the genome_array in the order of the genome file
	    print  STDERR "$name\t$genome_number\n" if ($DEBUG);
	    $genome_number++;
	}
	my $contigfile;
	unless (open ($contigfile, "<", $contig_file) )  {
	    die ("ERROR: cannot open file $contig_file.\n");
	}
	my ($save_input_separator) = $/;
	$/="\n>";
	while (my $line2 = <$contigfile>) {
	    (my $title, my $sequence) = split(/\n/, $line2, 2); # split the header line and sequence (very cool)
	    my @fields = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
	    my $id = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
	    $id =~ s/>\s*//; # remove leading > and spaces
	    if ($strip_version) {
		$id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    }
	    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	    #print STDERR "genome: $name contig: $id\n";
	    if (!defined $is_circular{$name}) {
		die ("ERROR: $name is a genome in the genomes fasta list file for $genomes_file_name but not in the contig topology file $topology_file!\n");
	    }
	    if (!defined $is_circular{$name}->{$id}) {
		die ("ERROR: $id is a contig in the genome fasta file for $name genome for $genomes_file_name but not in the contig topology file $topology_file!\n");
	    }
	    $genseq_hash{$name}->{$id} = $sequence;
	    $genseq_len{$name}->{$id} = length($sequence);
	    $title = ""; # clear the title for the next contig
	    $sequence = ""; #clear out the sequence for the next contig
	}
	$/ = $save_input_separator; # restore the input separator
	close ($contigfile);
	unless (open (TABLEFILE, "<", "$matchtable_file") )  {
	    die ("ERROR: cannot open file $matchtable_file.\n");
	}
	my $cluster_num = 1;
	while (my $line = <TABLEFILE>) {
	    chomp $line;
	    my @feat_names = split(/\t/, $line);  # split the scalar $line on tab
	    my $cluster_id = shift @feat_names;
	    if ($cluster_num != $cluster_id) {
		die ("ERROR: clusters are not sequentially ordered starting from 1: expecting $cluster_num but got $cluster_id\n");
	    }
	    $cluster_num++;
	    my $genome_tag = $name;
	    my $index = $genome_number - 1;
	    my $feat_name = $feat_names[$index];
	    if (($feat_name eq "----------") || ($feat_name eq "")) { #this is a placeholder and can be skipped
		next;
	    }
	    my $seq_len;
	    my $sequence;
	    if (!defined $feat_hash{$feat_name}) { # should not happen
		die ("ERROR: output_multifasta: gene identifier $feat_name in $matchtable_file is not in $att_file!\n");
	    }
	    if ($genome_tag ne $feat_hash{$feat_name}->{'gtag'}) {
		die ("ERROR: genome tag in $att_file ($feat_hash{$feat_name}->{'gtag'}) not the same as expected column in $matchtable_file ($genome_tag)");
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
	    my $contig_len = $genseq_len{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}};
	    print STDERR "$feat_name $feat_hash{$feat_name}->{'anno'} $feat_hash{$feat_name}->{'gtag'} $genome_tag\n" if ($DEBUG);
	    if (($fivep < 1) || ($threep < 1) || ($fivep > $contig_len) || ($threep > $contig_len)) {
		if ($is_circular{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}) {
		    if ($fivep <= $threep) {
			$seq_len = ($threep - $fivep) + 1;
			$sequence = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - 1));
			if (($fivep < 1) && ($threep >= 1)){
			    $sequence .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, $threep);
			} elsif (($threep > $contig_len) && ($fivep <= $contig_len)) {
			    $sequence .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, ($threep - $contig_len));
			} else {
			    die "ERROR: feature coordinates falling outside of contig boudaries (1:$contig_len) are not as expected: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
			}
		    } else {
			$seq_len = ($fivep - $threep) + 1;
			my $tmp_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - 1));
			if (($threep < 1) && ($fivep >= 1)){
			    $tmp_seq .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, $fivep);
			} elsif (($fivep > $contig_len) && ($threep <= $contig_len)) {
			    $tmp_seq .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, ($fivep - $contig_len));
			} else {
			    die "ERROR: feature coordinates falling outside of contig boudaries (1:$contig_len) are not as expected: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
			}
			$sequence = reverse($tmp_seq);
			$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		} else {
		    die "ERROR: feature coordinate falls outside of contig boudaries (1:$contig_len) and contig is not indicated to be circular: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
		}
	    } elsif ($fivep <= $threep) {
		$seq_len = ($threep - $fivep) + 1;
		$sequence = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - 1), $seq_len);
	    } else {
		$seq_len = ($fivep - $threep) + 1;
		my $tmp_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - 1), $seq_len);
		$sequence = reverse($tmp_seq);
		$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
	    }
	    $feat_hash{$feat_name}->{'len'} = $seq_len;
	    if ($seq_len <= 0) { #should not happen
		die ("ERROR: coordinates on contig sequence resulted in negative seq_len $seq_len for $feat_name!\n");
	    }
	    my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
	    if ($bad_count > 0) {
		die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
	    }
	    unless (open (OUTFILE, ">>$multifastadir/cluster_full_$cluster_id.fasta") )  {
		die ("ERROR: cannot open file $multifastadir/cluster_full_$cluster_id.fasta\n");
	    }
	    print OUTFILE ">$genome_tag\t$feat_name\n";
	    my $pos;
	    my $tmp_seq_len = $seq_len;
	    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
		print OUTFILE substr($sequence, $pos, 60), "\n";
		$tmp_seq_len -= 60;
	    }
	    print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
	    close (OUTFILE);
	}
	close (TABLEFILE);
	unless (open (PGGFILE, "<", "$pgg_file") )  {
	    die ("ERROR: cannot open file $pgg_file.\n");
	}
	while (my $line = <PGGFILE>) {
	    chomp $line;
	    my $cluster1;
	    my $cluster2;
	    my $whichend1;
	    my $whichend2;
	    my @edge_values = split(/\t/, $line);  # split the scalar $line on tab
	    my $edge_id = shift @edge_values;
	    my $edge_name = $edge_id;
	    if ($edge_id =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
		$edge_id = "edge".$1."_".$2."to".$3."_".$4;
		$cluster1 = $1;
		$cluster2 = $3;
		$whichend1 = $2;
		$whichend2 = $4;
	    } else {
		die ("ERROR: Bad edge formatting $edge_id in file $pgg_file.\n");
	    }
	    
	    my $index = $genome_number - 1;
	    my $edge_value = $edge_values[$index];
	    my $genome_tag = $name;
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
		die ("ERROR: output:multifasta cluster to feat_name mapping is in conflict $genome_tag $cluster1 $cluster2!\n");
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
	    my $seq_len;
	    my $sequence;
	    my $edge_5p;
	    my $edge_3p;
	    my $contig_len = $genseq_len{$genome_tag}->{$contig1};
	    if ($start1 < $end1) { # forward strand
		if ((($start1 > $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && !$is_circular{$genome_tag}->{$contig1}) {
		    die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		}
		if ((($start1 > $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && $is_circular{$genome_tag}->{$contig1}) { # this is a circular contig
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if (($start1 <= 0) || ($end1 <= 0)) {
			$start1 += $contig_len; # normalize edge coordinates to be > 0
			$end1 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if (($start2 <= 0) || ($end2 <= 0)) {
			$start2 += $contig_len; # normalize edge coordinates to be > 0
			$end2 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if ($start1 > $end2) { # the contigs on either side of the edge are on opposite ends of the contig
			if ($end1 <= $contig_len) {
			    $seq_len = ($contig_len - $end1) + ($start2 - 1);
			    $edge_5p = ($end1 < $contig_len) ? ($end1 + 1) : 1;
			    $edge_3p = ($end1 < $contig_len) ? (($start2 - 1) + $contig_len) : ($start2 - 1);
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			} else {
			    my $extra = $end1 - $contig_len;
			    $seq_len = ($start2 - 1) - $extra;
			    $edge_5p = $extra + 1;
			    $edge_3p = $start2 - 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			}
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    if ($end1 < $contig_len) {
				$sequence = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - $contig_len));
				$sequence .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, ($start2 - 1));
			    } else {
				$sequence = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - $contig_len), $seq_len);
			    }
			}
		    } else { #contigs on either end of the edge are on the same end of the contig
			$seq_len = ($start2 - $end1) - 1;
			$edge_5p = $end1 + 1;
			$edge_3p = $start2 - 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $seq_len);
			}
		    }
		} else { # contig is not circular
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if ($start1 > $end2) { # this should never happen for linear contigs
			die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		    } else { #contigs on either end of the edge are on the same end of the contig as expected
			$seq_len = ($start2 - $end1) - 1;
			$edge_5p = $end1 + 1;
			$edge_3p = $start2 - 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $seq_len);
			}
		    }
		}
	    } else { # reverse strand
		if ((($start1 < $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && !$is_circular{$genome_tag}->{$contig1}) {
		    die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		}
		if ((($start1 < $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && $is_circular{$genome_tag}->{$contig1}) { # this is a circular contig
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if (($start1 <= 0) || ($end1 <= 0)) {
			$start1 += $contig_len; # normalize edge coordinates to be > 0
			$end1 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if (($start2 <= 0) || ($end2 <= 0)) {
			$start2 += $contig_len; # normalize edge coordinates to be > 0
			$end2 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if ($start1 < $end2) { # the contigs on either side of the edge are on opposite ends of the contig
			if ($start2 <= $contig_len) {
			    $seq_len = ($contig_len - $start2) + ($end1 - 1);
			    $edge_5p = ($start2 < $contig_len) ? (($end1 - 1) + $contig_len) : ($end1 - 1);
			    $edge_3p = ($start2 < $contig_len) ? ($start2 + 1) : 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			} else {
			    my $extra = $start2 - $contig_len;
			    $seq_len = ($end1 - 1) - $extra;
			    $edge_5p = $end1 - 1;
			    $edge_3p = $extra + 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			}
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = "";
			    if ($start2 < $contig_len) {
				$tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - $contig_len));
				$tmp_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, ($end1 - 1));
			    } else {
				$tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - $contig_len), $seq_len);
			    }
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    } else { #contigs on either end of the edge are on the same end of the contig
			$seq_len = ($end1 - $start2) - 1;
			$edge_5p = $end1 - 1;
			$edge_3p = $start2 + 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $seq_len);
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    }
		} else { # contig is not circular
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if ($start1 < $end2) { # this should never happen for linear contigs
			die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		    } else { #contigs on either end of the edge are on the same end of the contig as expected
			$seq_len = ($end1 - $start2) - 1;
			$edge_5p = $end1 - 1;
			$edge_3p = $start2 + 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $seq_len);
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    }
		}
	    }
	    if ($cluster1 <= $cluster2) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
		if ($seq_len < 0) { #should not happen
		    die ("ERROR: coordinates on contig sequence reulted in negative seq_len $seq_len for $edge_value $feat_name1 $start1 $end1 $feat_name2 $start2 $end2!\n");
		}
		my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
		if ($bad_count > 0) {
		    die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
		}
		if ($sequence eq "") {
		    $sequence = "EMPTY";
		    $seq_len = 5;
		}
		unless (open (OUTFILE, ">>$multifastadir/full_$edge_id.fasta") )  {
		    die ("ERROR: cannot open file $multifastadir/full_$edge_id.fasta\n");
		}
		print OUTFILE ">$genome_tag\t$edge_5p\t$edge_3p\n";
		my $pos;
		my $tmp_seq_len = $seq_len;
		for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
		    print OUTFILE substr($sequence, $pos, 60), "\n";
		    $tmp_seq_len -= 60;
		}
		print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
		close (OUTFILE);
	    }
	}
	close (PGGFILE);
	$genseq_hash{$name} = (); #free up memory
	$genseq_hash{$name}->{"Placeholder"} = "EMPTY"; #this is just here to allow checking for duplicate genome names
    }
    close($infile);
    print  STDERR "$genome_number genomes\n\n";
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
	die ("ERROR: cannot open file $att_file.\n");
    }
    while (<ATTFILE>) {
	my @att_line = ();
	chomp;
	@att_line = split(/\t/, $_);  # split the scalar $line on tab
	$asmbl_id = $att_line[0];
	$feat_name = $att_line[1];
	$end5 = $att_line[2];
	$end3 = $att_line[3];
	$anno = $att_line[4];
	$tag = $att_line[5];
	if ($ignore_id eq $tag) {
	    #print STDERR "$ignore_id:$ignore_index:$tag:$feat_name\n" if ($DEBUG);
	    next; #skip over the genome to be ignored
	}
	#if ($use_multifasta && ($target_id ne $tag)) { #only read in the target genome - will use multiple fasta files for edges and clusters instead
	#    next; #skip over the genome to be ignored
	#}
	if ($asmbl_id eq "") {
	    print STDERR "ERROR: assembly id/contig id must not be empty/null in the gene attribute file\n$_\n";
	    $failed = 1;
	}
	if ($strip_version) {
	    $asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	}
	if (defined $feat_hash{$feat_name}) {
	    print STDERR "ERROR: $feat_name appears more than once in the gene attribute file $att_file!\n";
	    $failed = 1;
	}
	#if (!defined $genseq_hash{$tag}) {
	#    print STDERR "ERROR: $tag is a genome tag in the gene attribute file $att_file but not in the genome tag file $genomes_file_name!\n";
	#    $failed = 1;
	#}
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

    if ($remake_files) {
	unless (open (OUTMATCHFILE, ">", "$basedir/matchtable.txt") )  {
	    die ("ERROR: cannot open file $basedir/matchtable.txt!\n");
	}
	unless (open (SIZEFILE, ">", "$basedir/cluster_sizes.txt") ) {
	    die ("ERROR: cannot open file $basedir/cluster_sizes.txt!\n");
	}
    }
    unless (open (TABLEFILE, "<", "$matchtable_file") )  {
	die ("ERROR: cannot open file $matchtable_file.\n");
    }
    my $cluster_num = 1;
    my $reduced_cluster_num = 1;
    while (my $line = <TABLEFILE>) {
	chomp $line;
	my @feat_names = split(/\t/, $line);  # split the scalar $line on tab
	my $cluster_id = shift @feat_names;
	my $out_line = join("\t", @feat_names);
	if ($cluster_num != $cluster_id) {
	    die ("ERROR: clusters are not sequentially ordered starting from 1: expecting $cluster_num but got $cluster_id\n");
	}
	my %feat_pres = (); # key = sequence of feature, value = number of features with this sequence
	my $target_sequence = ""; # sequence for the target genome if specified
	my @tmp_array = @genome_array;
	my $genome_tag;
	my $gene_count = 0;
	my $single_genome = "";
	my $index = 0;
	my @genome_seqs = ();
	my $min = 10000000000;
	my $max = 0;
	my $sum = 0;
	my $sumsquared = 0;
	my @sizes = ();
	my $div_by_three = 0;
	my $seen = 0;
	if ($use_multifasta && !$no_stats) {
	    if (open (CLUSTERFILE, "<$multifastadir/cluster_full_$cluster_id.fasta") )  { #if the file isn't there it's because the cluster was empty and going away on the next iteration
		my ($save_input_separator) = $/;
		$/="\n>";
		while (my $line2 = <CLUSTERFILE>) {
		    (my $title, my $sequence) = split(/\n/, $line2, 2); # split the header line and sequence (very cool)
		    my @fields = split(/\t/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
		    $genome_tag = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
		    $genome_tag =~ s/>\s*//; # remove leading > and spaces
		    if ($strip_version) {
			$genome_tag =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
		    }
		    my $feat_name = $fields[1];
		    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
		    my $seq_len = length($sequence);
		    if ($genome_tag eq $ignore_id) {
			next; # need to get target genome info from actual fasta file not from multifasta file
		    }
		    if ($genome_tag eq $target_id) {
			die ("ERROR: $target_id is the target genome and should not be in the multifasta files\n"); # need to get target genome info from actual fasta file not from multifasta file
		    }
		    my $tag = shift @tmp_array;
		    while ($tag ne $genome_tag) {
			$index++;
			$tag = shift @tmp_array;
			if (!defined $tag) {
			    die ("ERROR: shifted off the end of the genome tag array while looking for $genome_tag with featname $feat_name for cluster $cluster_id and target $target_id\n");
			}
		    }
		    if (!defined $feat_hash{$feat_name}) { # should not happen
			die ("ERROR: process_amtchtable:cluster gene identifier $feat_name in $matchtable_file is not in $att_file!\n");
		    }
		    if ($genome_tag ne $feat_hash{$feat_name}->{'gtag'}) {
			die ("ERROR: genome tag in $att_file ($feat_hash{$feat_name}->{'gtag'}) not the same as expected column in $matchtable_file ($genome_tag)");
		    }
		    if (!defined $feat_hash{$feat_name}->{'contig'}) { # should not happen
			die ("ERROR: contig identifier was not assigned for $feat_name in $matchtable_file should have come from $att_file!\n");
		    }
		    if (!defined $genseq_hash{$genome_tag}) { # should not happen
			die ("ERROR: genome tag identifier was not assigned for $feat_name in $matchtable_file should have come from $att_file!\n");
		    }
		    $cluster_to_feat_hash{$genome_tag}->{$cluster_id} = $feat_name;
		    $feat_hash{$feat_name}->{'len'} = $seq_len;
		    if ($seq_len <= 0) { #should not happen
			die ("ERROR: from multifasta file seq_len $seq_len for $feat_name!\n");
		    }
		    my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
		    if ($bad_count > 0) {
			die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
		    }
		    if (($seq_len % 3) == 0) {
			$div_by_three++;
		    }
		    $genome_seqs[$index] = $sequence;
		    if ($compute_all || ($target_id ne "")) {
			if (!$gene_count) {
			    $single_genome = $genome_tag;
			}
			if ($genome_tag eq $target_id) {
			    $target_sequence = $sequence;
			}
		    }
		    if (($genome_tag ne $target_id) && ($align_all || $remake_files || $compute_all || ($target_id ne ""))) {
			if ($seq_len > $max) {
			    $max = $seq_len;
			}
			if ($seq_len < $min) {
			    $min = $seq_len;
			}
			$sum += $seq_len;
			$sumsquared += $seq_len * $seq_len;
			push @sizes, $seq_len;
		    }
		    if ($genome_tag ne $target_id) {
			$gene_count++;
		    }
		    $index++;
		    $title = ""; # clear the title for the next contig
		    $sequence = ""; #clear out the sequence for the next contig
		}		$/ = $save_input_separator; # restore the input separator
		close (CLUSTERFILE);
	    } else {
		print STDERR "cluster_full_$cluster_id.fasta not found\n";
	    }

	}
	@tmp_array = @genome_array;
	$index = 0;
	foreach my $feat_name (@feat_names) {
	    if (!$seen && ($ignore_index == $index)) {
		#print STDERR "$ignore_id:$ignore_index:$index:$feat_name\n" if ($DEBUG);
		$seen =1;
		next; # ignore the column corresponding to the genome to be ignored
	    }
	    $genome_tag = shift @tmp_array;
	    if (($feat_name eq "----------") || ($feat_name eq "")) { #this is a placeholder and can be skipped
		$index++;
		next;
	    }
	    if ($no_stats) {
		$gene_count++;
		$index++;
		next;
	    }
	    if ($use_multifasta && ($genome_tag ne $target_id)) {
		$index++;
		next; # processed these from multifasta file
	    }
	    my $seq_len;
	    my $sequence;
	    if (!defined $feat_hash{$feat_name}) { # should not happen
		die ("ERROR: process_matchtable:genome gene identifier $feat_name in $matchtable_file is not in $att_file!\n");
	    }
	    if ($genome_tag ne $feat_hash{$feat_name}->{'gtag'}) {
		die ("ERROR: genome tag in $att_file ($feat_hash{$feat_name}->{'gtag'}) not the same as expected column in $matchtable_file ($genome_tag)");
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
	    my $fivep_seq = "";
	    my $threep_seq = "";
	    my $contig_len = $genseq_len{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}};
	    print STDERR "$feat_name $feat_hash{$feat_name}->{'anno'} $feat_hash{$feat_name}->{'gtag'} $genome_tag\n" if ($DEBUG);
	    if (($fivep < 1) || ($threep < 1) || ($fivep > $contig_len) || ($threep > $contig_len)) {
		if ($is_circular{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}) {
		    if ($fivep <= $threep) {
			$seq_len = ($threep - $fivep) + 1;
			$sequence = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - 1));
			if (($fivep < 1) && ($threep >= 1)){
			    $sequence .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, $threep);
			    my $fivep_start = ($contig_len + $fivep) - $Gapped_Context;
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep_start - 1), $Gapped_Context);
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, $threep, $Gapped_Context);
			} elsif (($threep > $contig_len) && ($fivep <= $contig_len)) {
			    $sequence .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, ($threep - $contig_len));
			    my $fivep_start = $fivep - $Gapped_Context;
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep_start - 1), $Gapped_Context);
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - $contig_len), $Gapped_Context);
			} else {
			    die "ERROR: feature coordinates falling outside of contig boudaries (1:$contig_len) are not as expected: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
			}
		    } else {
			$seq_len = ($fivep - $threep) + 1;
			my $tmp_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - 1));
			if (($threep < 1) && ($fivep >= 1)){
			    $tmp_seq .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, $fivep);
			    my $threep_start = ($contig_len + $threep) - $Gapped_Context;
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep_start - 1), $Gapped_Context);
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, $fivep, $Gapped_Context);
			} elsif (($fivep > $contig_len) && ($threep <= $contig_len)) {
			    $tmp_seq .= substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, 0, ($fivep - $contig_len));
			    my $threep_start = $threep - $Gapped_Context;
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep_start - 1), $Gapped_Context);
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - $contig_len), $Gapped_Context);
			} else {
			    die "ERROR: feature coordinates falling outside of contig boudaries (1:$contig_len) are not as expected: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
			}
			$sequence = reverse($tmp_seq);
			$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		} else {
		    die "ERROR: feature coordinate falls outside of contig boudaries (1:$contig_len) and contig is not indicated to be circular: $genome_tag $feat_hash{$feat_name}->{'contig'} $fivep $threep\n";
		}
	    } elsif ($fivep <= $threep) {
		$seq_len = ($threep - $fivep) + 1;
		$sequence = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep - 1), $seq_len);
		if ($fivep <= $Gapped_Context) {
		    $feat_hash{$feat_name}->{'gapped'} = 1;
		} else {
		    my $fivep_start = $fivep - $Gapped_Context;
		    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($fivep_start - 1), $Gapped_Context);
		}
		if (($contig_len - $threep) < $Gapped_Context) {
		    $feat_hash{$feat_name}->{'gapped'} = 1;
		} else {
		    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, $threep, $Gapped_Context);
		}
	    } else {
		$seq_len = ($fivep - $threep) + 1;
		my $tmp_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep - 1), $seq_len);
		$sequence = reverse($tmp_seq);
		$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		if ($threep <= $Gapped_Context) {
		    $feat_hash{$feat_name}->{'gapped'} = 1;
		} else {
		    my $threep_start = $threep - $Gapped_Context;
		    $threep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, ($threep_start - 1), $Gapped_Context);
		}
		if (($contig_len - $fivep) < $Gapped_Context) {
		    $feat_hash{$feat_name}->{'gapped'} = 1;
		} else {
		    $fivep_seq = substr($genseq_hash{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}}, $fivep, $Gapped_Context);
		}
	    }
	    if (($fivep_seq =~ /NNNNN/) || ($threep_seq =~ /NNNNN/)) {
		$feat_hash{$feat_name}->{'gapped'} = 1;
	    }
	    $feat_hash{$feat_name}->{'len'} = $seq_len;
	    if ($seq_len <= 0) { #should not happen
		die ("ERROR: coordinates on contig sequence resulted in negative seq_len $seq_len for $feat_name!\n");
	    }
	    my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
	    if ($bad_count > 0) {
		die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
	    }
	    if (($seq_len % 3) == 0) {
		$div_by_three++;
	    }
	    $genome_seqs[$index] = $sequence;
	    if ($compute_all || ($target_id ne "")) {
		if (!$gene_count) {
		    $single_genome = $genome_tag;
		}
		if ($genome_tag eq $target_id) {
		    $target_sequence = $sequence;
		}
	    }
	    if (($genome_tag ne $target_id) && ($align_all || $remake_files || $compute_all || ($target_id ne ""))) {
		if ($seq_len > $max) {
		    $max = $seq_len;
		}
		if ($seq_len < $min) {
		    $min = $seq_len;
		}
		$sum += $seq_len;
		$sumsquared += $seq_len * $seq_len;
		push @sizes, $seq_len;
	    }
	    if ($genome_tag ne $target_id) {
		$gene_count++;
	    }
	    $index++;
	}
	my $mean;
	my $median;
	my $median_25;
	my $median_75;
	my $stddev;
	if ($align_all || $remake_files || $compute_all || ($target_id ne "")) {
	    if ($gene_count > 0) {
		if (((100 * $gene_count) / $genome_number) >= 95) {
		    $num_core_clus++;
		} elsif ($gene_count > 1) {
		    $num_shared_clus++;
		} else {
		    $num_size_one_clus++;
		}
		if ($no_stats) {
		    $max = $min = $mean = $median = $median_25 = $median_75 = $stddev = 0; # the cluster_sizes file will have less useful information in this case
		} else {
		    $mean = $sum / $gene_count;
		    @sizes = sort {$a <=> $b} @sizes;
		    $median = ($gene_count % 2) ? $sizes[($gene_count / 2)] : (($sizes[(($gene_count / 2) - 1)] + $sizes[($gene_count / 2)]) / 2);
		    if ($gene_count < 4) {
			$median_25 = $median_75 = $median;
		    } else {
			$median_25 = $sizes[int(($gene_count + 1) / 4)];
			$median_75 = $sizes[($gene_count - 1) - int(($gene_count + 1) / 4)];
		    }
		    $stddev = sqrt(($sumsquared - ($mean * $mean * $gene_count)) / (($gene_count > 1) ? ($gene_count - 1) : 1));
		}
		if ($remake_files) {
		    print SIZEFILE "$reduced_cluster_num\t$gene_count\t\t$min\t$max\t$median\t$mean\t$stddev\t\t$median_25\t$median_75\n"; # do not have "connectivity" or average %identity at this point
		    print OUTMATCHFILE "$reduced_cluster_num\t$out_line\n";
		}
	    } else {
		if ($target_sequence ne "") {
		    $max = $min = $mean = $median = $median_25 = $median_75 = length($target_sequence);
		    $stddev = 0;
		    $num_size_one_clus++;
		} else {
		    $num_reduced_clus++;
		}
	    }
	}
	if (!$suppress) {
	    unless (open (OUTFILE, ">$multifastadir/cluster_$cluster_id.fasta") )  {
		die ("ERROR: cannot open file $multifastadir/cluster_$cluster_id.fasta\n");
	    }
	    my $index = 0;
	    my $nr_allele_num = 0;
	    my $empty = 0;
	    foreach my $genome_tag (@genome_array) {
		if (defined $genome_seqs[$index]) {
		    my $sequence = $genome_seqs[$index];
		    my $seq_len = length($sequence);
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
			my $ambig_count = $sequence =~ tr/ACGTacgt//c;
			if (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20)) {
			    # do not include sequences with a lot of ambiguous base calls
			} elsif (($seq_len < ($median_25 - (0.1 * $median))) || ($seq_len > ($median_75 + (0.1 * $median)))) {
			    # do not include length outliers in the multifasta file
			} elsif ($sequence eq "EMPTY") {
			    $empty = 1;
			} else {
			    $nr_allele_num++;
			    print OUTFILE ">$genome_tag\n";
			    my $pos;
			    my $tmp_seq_len = $seq_len;
			    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
				print OUTFILE substr($sequence, $pos, 60), "\n";
				$tmp_seq_len -= 60;
			    }
			    print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
			}
		    }
		}
		$index++;
	    }
	    push (@mf_files, "$multifastadir/cluster_$cluster_id.fasta\t$nr_allele_num\t$empty");
	    close (OUTFILE);
	} else {
	    my $index = 0;
	    foreach my $genome_tag (@genome_array) {
		if (defined $genome_seqs[$index]) {
		    my $sequence = $genome_seqs[$index];
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
		    }
		}
		$index++;
	    }
	}
	if ($gene_count > 1) {
	    $total_clus_pgg++;
	    $total_clus_alle_pgg += $gene_count;
	}
	if ($no_stats) {
	    if ($gene_count > 0) {
		$renumber[$cluster_num] = $reduced_cluster_num;
		$reduced_cluster_num++;
	    } else {
		$renumber[$cluster_num] = 0;
	    }
	    $cluster_num++;
	    next; # do not generate the anomalies file or multiple sequence alignments
	}
	if ($compute_all || ($target_id ne "")) {
	    if ($compute_all && ($gene_count > 0)) {
		$index = 0;
		foreach my $genome_tag (@genome_array) {
		    if (($gene_count == 1) && ($single_genome eq $genome_tag)) {
			$uniq_clus{$genome_tag}++;
		    } elsif (defined $genome_seqs[$index]) {
			my $feat_name = $cluster_to_feat_hash{$genome_tag}->{$cluster_id};
			if (!defined $feat_pres{$genome_seqs[$index]}) {
			    print STDERR "$genome_tag:$index:$feat_name:\n$genome_seqs[$index]\n" if ($DEBUG);
			}
			my $seq_len = $feat_hash{$feat_name}->{'len'};
			my $ambig_count = $genome_seqs[$index] =~ tr/ACGTacgt//c;
			if (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20)) {
			    # do not include sequences with a lot of ambiguous base calls
			    $gapped_clus{$genome_tag}++;
			    print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tgapped_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			} else {
			    if (defined $feat_hash{$feat_name}->{'gapped'}) {
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tgapped_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    }
			    if ($feat_pres{$genome_seqs[$index]} == 1) {
				if (((100 * $gene_count) / $genome_number) >= 75) {
				    $uniq_clus_alle_75_100{$genome_tag}++;
				} elsif (((100 * $gene_count) / $genome_number) <= 25) {
				    $uniq_clus_alle_0_25{$genome_tag}++;
				} else {
				    $uniq_clus_alle_25_75{$genome_tag}++;
				}
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tuniq_clus_allele\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    } else {
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tidentical_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    }
			    if (($seq_len < $min) && ($seq_len < ($median - (0.02 * $median)))) {
				$very_short_clus{$genome_tag}++;
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tvery_short_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    } elsif ($seq_len < ($median_25 - (0.1 * $median))) {
				$short_clus{$genome_tag}++;
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tshort_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    } elsif (($seq_len > $max) && ($seq_len > ($median + (0.02 * $median)))) {
				$very_long_clus{$genome_tag}++;
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tvery_long_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    } elsif ($seq_len > ($median_75 + (0.1 * $median))) {
				$long_clus{$genome_tag}++;
				print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tlong_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    }
			    if ($div_by_three > (0.7 * $gene_count)) { # best approximation for frameshift
				if (($seq_len % 3) != 0) {
				    $frameshift{$genome_tag}++;
				    print DIFFFILE "$genome_tag\t$feat_hash{$feat_name}->{'contig'}\tframeshift_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
				}
			    }
			}
			$total_clus{$genome_tag}++;
		    } else {
			if (((100 * $gene_count) / $genome_number) >= 75) {
			    if (defined $single_copy_core{$cluster_num}) {
				$miss_sing_core{$genome_tag}++;
			    } else {
				$missing_75c{$genome_tag}++;
			    }
			}
		    }
		    $index++;
		}
	    } else {
		if (($gene_count == 0) && ($target_sequence ne "")) {
		    $uniq_clus{$target_id}++;
		} elsif ($target_sequence ne "") {
		    my $feat_name = $cluster_to_feat_hash{$target_id}->{$cluster_id};
		    my $seq_len = $feat_hash{$feat_name}->{'len'};
		    my $ambig_count = $target_sequence =~ tr/ACGTacgt//c;
		    if (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20)) {
			# do not include sequences with a lot of ambiguous base calls
			$gapped_clus{$target_id}++;
			print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tgapped_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
		    } else {
			if (defined $feat_hash{$feat_name}->{'gapped'}) {
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tgapped_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			}
			if ($feat_pres{$target_sequence} == 1) {
			    if (((100 * $gene_count) / $genome_number) >= 75) {
				$uniq_clus_alle_75_100{$target_id}++;
			    } elsif (((100 * $gene_count) / $genome_number) <= 25) {
				$uniq_clus_alle_0_25{$target_id}++;
			    } else {
				$uniq_clus_alle_25_75{$target_id}++;
			    }
			    if ($align_new) {
				print STDERR "$target_id\t$cluster_id\t$seq_len\t$median_25\t$median\t$median_75\n" if ($DEBUG);
				if (($seq_len >= ($median_25 - (0.1 * $median))) && ($seq_len <= ($median_75 + (0.1 * $median)))) {
				    my $mf_file = "$multifastadir/cluster_$cluster_id.fasta";
				    my $msa_file = "$multifastadir/cluster_$cluster_id.afa";
				    my $stats_file = "$multifastadir/cluster_$cluster_id.stats";
				    my $empty_file = "$multifastadir/cluster_$cluster_id.empty";
				    my $max_all;
				    my $median_target;
				    my $max_cols_all;
				    my $cols_target;
				    print STDERR "Good length\t$stats_file\t$msa_file\t$mf_file\n" if ($DEBUG);
				    if ((-s $stats_file) && ((-s $msa_file) || (-s $mf_file))){
					unless (open (STATSFILE, "<", $stats_file) )  {
					    die ("ERROR: cannot open file $stats_file\n");
					}
					while (my $line = <STATSFILE>)  {
					    chomp $line;
					    (my @fields) = split(/\t/, $line);  # split the scalar $line on tab
					    if ($fields[0] eq "All") {
						$max_all = $fields[4];
					    } elsif ($fields[0] eq "UniqueAlleleCount") {
						$max_cols_all = $fields[4];
					    }
					}
					close(STATSFILE);
					my $seq_file = "$multifastadir/cluster_TMP_" . $$ . "_$target_id.$cluster_id.fasta";
					my $combined_file = "$multifastadir/cluster_TMP_" . $$ . "_$target_id.$cluster_id.afa";
					my $com_stats_file = "$multifastadir/cluster_TMP_" . $$ . "_$target_id.$cluster_id.stats";
					unless (open (OUTFILE, ">", $seq_file) )  {
					    die ("ERROR: cannot open file $seq_file\n");
					}
					print OUTFILE ">$target_id\n";
					my $pos;
					my $tmp_seq_len = $seq_len;
					for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
					    print OUTFILE substr($target_sequence, $pos, 60), "\n";
					    $tmp_seq_len -= 60;
					}
					print OUTFILE substr($target_sequence, $pos, $tmp_seq_len), "\n";
					close (OUTFILE);
					if (-s $msa_file) {
					    my $muscle_args = " -profile -in1 $msa_file -in2 $seq_file -out $combined_file -diags -quiet -verbose";
					    my $muscle_exec = $muscle_path . $muscle_args;
					    `/usr/bin/time -o tmp_cpu_stats -v $muscle_exec`;
					    `echo "***$muscle_exec***length($seq_len)" >> $cpu_name`;
					    `cat tmp_cpu_stats >> $cpu_name`;
					    `rm tmp_cpu_stats`;
					    &bash_error_check($muscle_exec, $?, $!);
					} elsif (-s $mf_file) {
					    `cat $mf_file >> $seq_file`;
					    my $muscle_args = " -in $seq_file -out $combined_file -diags -quiet -verbose";
					    my $muscle_exec = $muscle_path . $muscle_args;
					    `/usr/bin/time -o tmp_cpu_stats -v $muscle_exec`;
					    `echo "***$muscle_exec***length($seq_len)" >> $cpu_name`;
					    `cat tmp_cpu_stats >> $cpu_name`;
					    `rm tmp_cpu_stats`;
					    &bash_error_check($muscle_exec, $?, $!);
					}
					my $stats_path = "$bin_directory/summarize_alignment.R";
					my $stats_args = " $combined_file $target_id";
					my $stats_exec = $rscript_path . " " . $stats_path . $stats_args;
					`/usr/bin/time -o tmp_cpu_stats -v $stats_exec > $com_stats_file`;
					`echo "***$stats_exec***" >> $cpu_name`;
					`cat tmp_cpu_stats >> $cpu_name`;
					`rm tmp_cpu_stats`;
					&bash_error_check("$stats_exec > $com_stats_file", $?, $!);
					unless (open (STATSFILE, "<", $com_stats_file) )  {
					    die ("ERROR: cannot open file $com_stats_file\n");
					}
					while (my $line = <STATSFILE>)  {
					    print STDERR $line if ($DEBUG);
					    chomp $line;
					    (my @fields) = split(/\t/, $line);  # split the scalar $line on tab
					    if ($fields[0] eq $target_id) {
						$median_target = $fields[2];
						$cols_target = $fields[27];
					    }
					}
					close(STATSFILE);
					if ($keep_divergent_alignments) {
					    `mv $combined_file $keep_divergent_alignments/cluster_$target_id.$cluster_id.afa`;
					    `rm $seq_file $com_stats_file`;
					} else {
					    `rm $seq_file $combined_file $com_stats_file`;
					}
					print STDERR "$target_id\t$cluster_id\t$median_target\t$max_all\t$cols_target\t$max_cols_all\n" if ($DEBUG);
					if ($median_target > $max_all) {
					    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tdivergent_clus_allele\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
					    $distant_clus_alle{$target_id}++;
					}
					if ($cols_target > $max_cols_all) {
					    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tdivergent_columns_clus_allele\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
					    $column_clus_alle{$target_id}++;
					}
				    }
				}
			    }
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tuniq_clus_allele\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			} else {
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tidentical_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			}
			if (($seq_len < $min) && ($seq_len < ($median - (0.02 * $median)))) {
			    $very_short_clus{$target_id}++;
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tvery_short_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			} elsif ($seq_len < ($median_25 - (0.1 * $median))) {
			    $short_clus{$target_id}++;
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tshort_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			} elsif (($seq_len > $max) && ($seq_len > ($median + (0.02 * $median)))) {
			    $very_long_clus{$target_id}++;
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tvery_long_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			} elsif ($seq_len > ($median_75 + (0.1 * $median))) {
			    $long_clus{$target_id}++;
			    print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tlong_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			}
			if ($div_by_three > (0.7 * $gene_count)) { # best approximation for frameshift
			    if (($seq_len % 3) != 0) {
				$frameshift{$target_id}++;
				print DIFFFILE "$target_id\t$feat_hash{$feat_name}->{'contig'}\tframeshift_clus\t$feat_hash{$feat_name}->{'5p'}\t$feat_hash{$feat_name}->{'3p'}\t$feat_hash{$feat_name}->{'len'}\t$feat_name\n";
			    }
			}
		    }
		    $total_clus{$target_id}++;
		} else {
		    if (((100 * $gene_count) / $genome_number) >= 75) {
			if (defined $single_copy_core{$cluster_num}) {
			    $miss_sing_core{$target_id}++;
			} else {
			    $missing_75c{$target_id}++;
			}
		    }
		}
	    }
	}
	if ($gene_count > 0) {
	    $renumber[$cluster_num] = $reduced_cluster_num;
	    $reduced_cluster_num++;
	} else {
	    $renumber[$cluster_num] = 0;
	}
	$cluster_num++;
    }
    if ($remake_files) {
	close (OUTMATCHFILE);
	close(SIZEFILE);
    }
    close (TABLEFILE);
    return;
}

sub process_pgg {

    if ($remake_files) {
	unless (open (OUTPGGFILE, ">", "$basedir/pgg.txt") )  {
	    die ("ERROR: cannot open file $basedir/pgg.txt!\n");
	}
	unless (open (SIZEFILE, ">", "$basedir/edge_sizes.txt") ) {
	    die ("ERROR: cannot open file $basedir/edge_sizes.txt!\n");
	}
    }
    unless (open (PGGFILE, "<", "$pgg_file") )  {
	die ("ERROR: cannot open file $pgg_file.\n");
    }
    while (my $line = <PGGFILE>) {
	chomp $line;
	my $cluster1;
	my $cluster2;
	my $whichend1;
	my $whichend2;
	my @edge_values = split(/\t/, $line);  # split the scalar $line on tab
	my $edge_id = shift @edge_values;
	my $alt_edge_id;
	my $edge_name = $edge_id;
	my $out_line = join("\t", @edge_values);
	my $out_cluster1;
	my $out_cluster2;
	if ($edge_id =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
	    $edge_id = "edge".$1."_".$2."to".$3."_".$4;
	    $alt_edge_id = "edge".$3."_".$4."to".$1."_".$2;
	    $cluster1 = $1;
	    $cluster2 = $3;
	    $whichend1 = $2;
	    $whichend2 = $4;
	    $out_cluster1 = $renumber[$cluster1];
	    $out_cluster2 = $renumber[$cluster2];
	} else {
	    die ("ERROR: Bad edge formatting $edge_id in file $pgg_file.\n");
	}
	#print STDERR "$edge_name:$edge_id\n";
	my %feat_pres = (); # key = sequence of feature, value = number of features with this sequence
	my $target_sequence = ""; # sequence for the target genome if specified
	my $gene_count = 0;
	my $single_genome = "";
	my $index = 0;
	my @genome_seqs = ();
	my @tmp_array = @genome_array;
	my $genome_tag;
	my $min = 10000000000;
	my $max = 0;
	my $sum = 0;
	my $sumsquared = 0;
	my @sizes = ();
	my $seen = 0;
	if ($use_multifasta && !$no_stats) {
	    my $edge_file = ($cluster1 < $cluster2) ? "full_$edge_id.fasta" : "full_$alt_edge_id.fasta";
	    if (open (EDGEFILE, "<$multifastadir/$edge_file") )  { #if the file isn't there it's because the edge was empty and going away on the next iteration
		my ($save_input_separator) = $/;
		$/="\n>";
		while (my $line2 = <EDGEFILE>) {
		    (my $title, my $sequence) = split(/\n/, $line2, 2); # split the header line and sequence (very cool)
		    my @fields = split(/\t/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @line
		    $genome_tag = $fields[0]; # unique orf identifier is in column 0, com_name is in rest
		    $genome_tag =~ s/>\s*//; # remove leading > and spaces
		    if ($strip_version) {
			$genome_tag =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
		    }
		    my $edge_5p = $fields[1];
		    my $edge_3p = $fields[2];
		    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
		    my $seq_len = length($sequence);
		    if ($genome_tag eq $ignore_id) {
			next; # need to get target genome info from actual fasta file not from multifasta file
		    }
		    if ($genome_tag eq $target_id) {
			die ("ERROR: $target_id is the target genome and should not be in the multifasta files\n"); # need to get target genome info from actual fasta file not from multifasta file
		    }
		    my $tag = shift @tmp_array;
		    while ($tag ne $genome_tag) {
			$index++;
			$tag = shift @tmp_array;
			if (!defined $tag) {
			    die ("ERROR: shifted off the end of the genome tag array while looking for $genome_tag with edge $edge_name $edge_5p $edge_3p and target $target_id\n");
			}
		    }
		    my $feat_name1 = $cluster_to_feat_hash{$genome_tag}->{$cluster1};
		    my $feat_name2 = $cluster_to_feat_hash{$genome_tag}->{$cluster2};
		    my $contig1 = $feat_hash{$feat_name1}->{'contig'};
		    my $contig2 = $feat_hash{$feat_name2}->{'contig'};
		    if ((!defined $feat_name1) || (!defined $feat_name2) || (!defined $contig1) || (!defined $contig2)) { # should not happen
			die ("ERROR: process_pgg:edge cluster to feat_name mapping is in conflict $genome_tag $cluster1 $cluster2!\n");
		    }
		    if ($contig1 ne $contig2) { # should not happen
			die ("ERROR: for edge $edge_id cluster features $feat_name1:$feat_name2 are not on the same contig $contig1:$contig2");
		    }
		    if (($genome_tag ne $feat_hash{$feat_name1}->{'gtag'}) || ($genome_tag ne $feat_hash{$feat_name2}->{'gtag'})) { # should not happen
			die ("ERROR: Inconsistency in genome tag between $att_file, $matchtable_file, and $pgg_file for $feat_name1 and $feat_name2");
		    }
		    if (!defined $genseq_hash{$genome_tag}) { # should not happen
			die ("ERROR: genome tag identifier was not assigned for $pgg_file should have come from $att_file!\n");
		    }
		    $edge_hash{$genome_tag . $edge_id} = {};
		    $edge_hash{$genome_tag . $edge_id}->{'gtag'} = $genome_tag;
		    $edge_hash{$genome_tag . $edge_id}->{'contig'} = $contig1;
		    if ($sequence eq "EMPTY") {
			$seq_len = 0;
		    }
		    $edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
		    $edge_hash{$genome_tag . $edge_id}->{'5p'} = $edge_5p;
		    $edge_hash{$genome_tag . $edge_id}->{'3p'} = $edge_3p;
		    if (($genome_tag ne $target_id) && ($align_all || $remake_files || $compute_all || ($target_id ne ""))) {
			if ($seq_len > $max) {
			    $max = $seq_len;
			}
			if ($seq_len < $min) {
			    $min = $seq_len;
			}
			$sum += $seq_len;
			$sumsquared += $seq_len * $seq_len;
			push @sizes, $seq_len;
		    }
		    if ($cluster1 <= $cluster2) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
			if ($seq_len < 0) { #should not happen
			    die ("ERROR: zero seq_len $seq_len for $edge_id $feat_name1 $edge_5p $feat_name2 $edge_3p!\n");
			}
			if ($sequence ne "EMPTY") {
			    my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
			    if ($bad_count > 0) {
				die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
			    }
			}
			$genome_seqs[$index] = $sequence;
			if ($compute_all || ($target_id ne "")) {
			    if (!$gene_count) {
				$single_genome = $genome_tag;
			    }
			    if ($genome_tag eq $target_id) {
				$target_sequence = $sequence;
			    }
			}
		    }
		    if ($genome_tag ne $target_id) {
			$gene_count++;
		    }
		    $index++;
		    $title = ""; # clear the title for the next contig
		    $sequence = ""; #clear out the sequence for the next contig
		}		$/ = $save_input_separator; # restore the input separator
		close (EDGEFILE);
	    } else {
		print STDERR "$edge_file not found\n";
	    }

	}
	@tmp_array = @genome_array;
	$index = 0;
	foreach my $edge_value (@edge_values) {
	    if (!$seen && ($ignore_index == $index)) {
		print STDERR "$ignore_id:$ignore_index:$index:$edge_value\n" if ($DEBUG);
		$seen = 1;
		next; # ignore the column corresponding to the genome to be ignored
	    }
	    $genome_tag = shift @tmp_array;
	    if ($edge_value == 0) { #this is a placeholder and can be skipped
		$index++;
		next;
	    }
	    if ($no_stats) {
		$gene_count++;
		$index++;
		next;
	    }
	    if ($use_multifasta && ($genome_tag ne $target_id)) {
		$index++;
		next; # processed these from multifasta file
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
		die ("ERROR: process_pgg:pgg cluster to feat_name mapping is in conflict $genome_tag $cluster1 $cluster2!\n");
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
	    $edge_hash{$genome_tag . $edge_id} = {};
	    $edge_hash{$genome_tag . $edge_id}->{'gtag'} = $genome_tag;
	    $edge_hash{$genome_tag . $edge_id}->{'contig'} = $contig1;
	    my $fivep_seq = "";
	    my $threep_seq = "";
	    my $seq_len;
	    my $sequence;
	    my $contig_len = $genseq_len{$genome_tag}->{$contig1};
	    if ($start1 < $end1) { # forward strand
		if ((($start1 > $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && !$is_circular{$genome_tag}->{$contig1}) {
		    die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		}
		if ((($start1 > $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && $is_circular{$genome_tag}->{$contig1}) { # this is a circular contig
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if (($start1 <= 0) || ($end1 <= 0)) {
			$start1 += $contig_len; # normalize edge coordinates to be > 0
			$end1 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if (($start2 <= 0) || ($end2 <= 0)) {
			$start2 += $contig_len; # normalize edge coordinates to be > 0
			$end2 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if ($start1 > $end2) { # the contigs on either side of the edge are on opposite ends of the contig
			if ($end1 <= $contig_len) {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1), $Gapped_Context);
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - $Gapped_Context), $Gapped_Context);
			    $seq_len = ($contig_len - $end1) + ($start2 - 1);
			    $edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			    $edge_hash{$genome_tag . $edge_id}->{'5p'} = ($end1 < $contig_len) ? ($end1 + 1) : 1;
			    $edge_hash{$genome_tag . $edge_id}->{'3p'} = ($end1 < $contig_len) ? (($start2 - 1) + $contig_len) : ($start2 - 1);
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			} else {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 -1), $Gapped_Context);
			    my $extra = $end1 - $contig_len;
			    if ($extra >= $Gapped_Context) {
				$fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($extra - $Gapped_Context), $Gapped_Context);
			    } else {
				$fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($extra - $Gapped_Context));
				$fivep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $extra);
			    }
			    $seq_len = ($start2 - 1) - $extra;
			    $edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			    $edge_hash{$genome_tag . $edge_id}->{'5p'} = $extra + 1;
			    $edge_hash{$genome_tag . $edge_id}->{'3p'} = $start2 - 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			}
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    if ($end1 < $contig_len) {
				$sequence = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - $contig_len));
				$sequence .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, ($start2 - 1));
			    } else {
				$sequence = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - $contig_len), $seq_len);
			    }
			}
		    } else { #contigs on either end of the edge are on the same end of the contig
			if ((($start2 - 1) + $Gapped_Context) > $contig_len) {
			    my $extra = (($start2 - 1) + $Gapped_Context) - $contig_len;
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1));
			    $threep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $extra);
			} else {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1), $Gapped_Context);
			}
			if ($end1 >= $Gapped_Context) {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $end1 - $Gapped_Context, $Gapped_Context);
			} else {
			    my $extra = $end1 - $Gapped_Context;
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $extra);
			    $fivep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $end1);
			}
			$seq_len = ($start2 - $end1) - 1;
			$edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			$edge_hash{$genome_tag . $edge_id}->{'5p'} = $end1 + 1;
			$edge_hash{$genome_tag . $edge_id}->{'3p'} = $start2 - 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $seq_len);
			}
		    }
		} else { # contig is not circular
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if ($start1 > $end2) { # this should never happen for linear contigs
			die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		    } else { #contigs on either end of the edge are on the same end of the contig as expected
			if ((($start2 - 1) + $Gapped_Context) > $contig_len) {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1));
			} else {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1), $Gapped_Context);
			}
			if ($end1 >= $Gapped_Context) {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $end1 - $Gapped_Context, $Gapped_Context);
			} else {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, 0, $end1);
			}
			$seq_len = ($start2 - $end1) - 1;
			$edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			$edge_hash{$genome_tag . $edge_id}->{'5p'} = $end1 + 1;
			$edge_hash{$genome_tag . $edge_id}->{'3p'} = $start2 - 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $seq_len);
			}
		    }
		}
	    } else { # reverse strand
		if ((($start1 < $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && !$is_circular{$genome_tag}->{$contig1}) {
		    die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		}
		if ((($start1 < $end2) || ($start2 <= 0) || ($end1 <= 0) || ($start2 > $contig_len) || ($end1 > $contig_len)) && $is_circular{$genome_tag}->{$contig1}) { # this is a circular contig
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if (($start1 <= 0) || ($end1 <= 0)) {
			$start1 += $contig_len; # normalize edge coordinates to be > 0
			$end1 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if (($start2 <= 0) || ($end2 <= 0)) {
			$start2 += $contig_len; # normalize edge coordinates to be > 0
			$end2 += $contig_len; # normalize edge coordinates to be > 0
		    }
		    if ($start1 < $end2) { # the contigs on either side of the edge are on opposite ends of the contig
			if ($start2 <= $contig_len) {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - $Gapped_Context), $Gapped_Context);
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - 1), $Gapped_Context);
			    $seq_len = ($contig_len - $start2) + ($end1 - 1);
			    $edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			    $edge_hash{$genome_tag . $edge_id}->{'5p'} = ($start2 < $contig_len) ? (($end1 - 1) + $contig_len) : ($end1 - 1);
			    $edge_hash{$genome_tag . $edge_id}->{'3p'} = ($start2 < $contig_len) ? ($start2 + 1) : 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			} else {
			    my $extra = $start2 - $contig_len;
			    if ($extra >= $Gapped_Context) {
				$threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($extra - $Gapped_Context), $Gapped_Context);
			    } else {
				$threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($extra - $Gapped_Context));
				$threep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $extra);
			    }
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - 1), $Gapped_Context);
			    $seq_len = ($end1 - 1) - $extra;
			    $edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			    $edge_hash{$genome_tag . $edge_id}->{'5p'} = $end1 - 1;
			    $edge_hash{$genome_tag . $edge_id}->{'3p'} = $extra + 1;
			    #print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			}
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = "";
			    if ($start2 < $contig_len) {
				$tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - $contig_len));
				$tmp_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, ($end1 - 1));
			    } else {
				$tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - $contig_len), $seq_len);
			    }
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    } else { #contigs on either end of the edge are on the same end of the contig
			if ((($start2 - 1) + $Gapped_Context) > $contig_len) {
			    my $extra = (($start2 - 1) + $Gapped_Context) - $contig_len;
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1));
			    $threep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $extra);
			} else {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($start2 - 1), $Gapped_Context);
			}
			if ($end1 >= $Gapped_Context) {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $end1 - $Gapped_Context, $Gapped_Context);
			} else {
			    my $extra = $end1 - $Gapped_Context;
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $extra);
			    $fivep_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, 0, $end1);
			}
			$seq_len = ($end1 - $start2) - 1;
			$edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			$edge_hash{$genome_tag . $edge_id}->{'5p'} = $end1 - 1;
			$edge_hash{$genome_tag . $edge_id}->{'3p'} = $start2 + 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $seq_len);
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    }
		} else { # contig is not circular
		    #print STDERR "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n" if ($DEBUG);
		    if ($start1 < $end2) { # this should never happen for linear contigs
			die ("ERROR: bad edge coordinates for noncircular contig: $genome_tag $contig1($contig_len) $edge_id $edge_value $feat_name1 $feat_name2 $start1 $end1 $start2 $end2\n");
		    } else { #contigs on either end of the edge are on the same end of the contig as expected
			if ((($end1 - 1) + $Gapped_Context) > $contig_len) {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - 1));
			} else {
			    $threep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, ($end1 - 1), $Gapped_Context);
			}
			if ($start2 >= $Gapped_Context) {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2 - $Gapped_Context, $Gapped_Context);
			} else {
			    $fivep_seq = substr($genseq_hash{$genome_tag}->{$contig1}, 0, $start2);
			}
			$seq_len = ($end1 - $start2) - 1;
			$edge_hash{$genome_tag . $edge_id}->{'len'} = $seq_len;
			$edge_hash{$genome_tag . $edge_id}->{'5p'} = $end1 - 1;
			$edge_hash{$genome_tag . $edge_id}->{'3p'} = $start2 + 1;
			#print STDERR "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n" if ($DEBUG);
			if ($seq_len <= 0) {
			    $seq_len = 0;
			    $sequence = "";
			} else {
			    my $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $seq_len);
			    $sequence = reverse($tmp_seq);
			    $sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
			}
		    }
		}
	    }
	    if (($fivep_seq =~ /NNNNN/) || ($threep_seq =~ /NNNNN/)) {
		$edge_hash{$genome_tag . $edge_id}->{'gapped'} = 1;
	    }
	    if (($genome_tag ne $target_id) && ($align_all || $remake_files || $compute_all || ($target_id ne ""))) {
		if ($seq_len > $max) {
		    $max = $seq_len;
		}
		if ($seq_len < $min) {
		    $min = $seq_len;
		}
		$sum += $seq_len;
		$sumsquared += $seq_len * $seq_len;
		push @sizes, $seq_len;
	    }
	    if ($cluster1 <= $cluster2) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
		if ($seq_len < 0) { #should not happen
		    die ("ERROR: coordinates on contig sequence reulted in negative seq_len $seq_len for $edge_value $feat_name1 $start1 $end1 $feat_name2 $start2 $end2!\n");
		}
		my $bad_count = $sequence =~ tr/AGCTYRWSKMDVHBNagctyrwskmdvhbn//c;
		if ($bad_count > 0) {
		    die ("ERROR: Unexpected character not in [AGCTYRWSKMDVHBNagctyrwskmdvhbn] found in genome fasta sequence!\n$sequence\n");
		}
		if ($sequence eq "") {
		    $sequence = "EMPTY";
		}
		$genome_seqs[$index] = $sequence;
		if ($compute_all || ($target_id ne "")) {
		    if (!$gene_count) {
			$single_genome = $genome_tag;
		    }
		    if ($genome_tag eq $target_id) {
			$target_sequence = $sequence;
		    }
		}
	    }
	    if ($genome_tag ne $target_id) {
		$gene_count++;
	    }
	    $index++;
	}
	my $mean;
	my $median;
	my $median_25;
	my $median_75;
	my $stddev;
	if ($align_all || $remake_files || $compute_all || ($target_id ne "")) {
	    if ($gene_count > 0) {
		if ($cluster1 <= $cluster2) { # only count edge in one direction
		    if (((100 * $gene_count) / $genome_number) >= 95) {
			$num_core_edge++;
		    } elsif ($gene_count > 1) {
			$num_shared_edge++;
		    } else {
			$num_size_one_edge++;
		    }
		}
		if ($no_stats) {
		    $max = $min = $mean = $median = $median_25 = $median_75 = $stddev = 0; # the cluster_sizes file will have less useful information in this case
		} else {
		    $mean = $sum / $gene_count;
		    @sizes = sort {$a <=> $b} @sizes;
		    $median = ($gene_count % 2) ? $sizes[($gene_count / 2)] : (($sizes[(($gene_count / 2) - 1)] + $sizes[($gene_count / 2)]) / 2);
		    if ($gene_count < 4) {
			$median_25 = $median_75 = $median;
		    } else {
			$median_25 = $sizes[int(($gene_count + 1) / 4)];
			$median_75 = $sizes[($gene_count - 1) - int(($gene_count + 1) / 4)];
		    }
		    $stddev = sqrt(($sumsquared - ($mean * $mean * $gene_count)) / (($gene_count > 1) ? ($gene_count - 1) : 1));
		}
		if ($remake_files) {
		    print SIZEFILE "($out_cluster1", "_$whichend1,$out_cluster2", "_$whichend2)\t$gene_count\t\t$min\t$max\t$median\t$mean\t$stddev\t\t$median_25\t$median_75\n"; #  do not have "connectivity" or average %identity at this point
		    print OUTPGGFILE "($out_cluster1", "_$whichend1,$out_cluster2", "_$whichend2)\t$out_line\n";
		}
	    } else {
		if ($target_sequence ne "") {
		    $max = $min = $mean = $median = $median_25 = $median_75 = length($target_sequence);
		    $stddev = 0;
		    if ($cluster1 <= $cluster2) { # only count edge in one direction
			$num_size_one_edge++;
		    }
		} else {
		    if ($cluster1 <= $cluster2) { # only count edge in one direction
			$num_reduced_edge++;
		    }
		}
	    }
	}
	if (!$suppress && ($cluster1 <= $cluster2)) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
	    unless (open (OUTFILE, ">$multifastadir/$edge_id.fasta") )  {
		die ("ERROR: cannot open file $multifastadir/$edge_id.fasta\n");
	    }
	    my $index = 0;
	    my $nr_allele_num = 0;
	    my $empty = 0;
	    foreach my $genome_tag (@genome_array) {
		if (defined $genome_seqs[$index]) {
		    my $sequence = $genome_seqs[$index];
		    my $seq_len = length($sequence);
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
			my $ambig_count;
			if ($sequence eq "EMPTY") {
			    $sequence = "";
			    $seq_len = 0;
			}
			$ambig_count = $sequence =~ tr/ACGTacgt//c;
			if ($sequence eq "") {
			    $empty = 1;
			} elsif (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20)) {
			    # do not include sequences with a lot of ambiguous base calls
			} elsif (($seq_len < ($median_25 - (0.1 * $median))) || ($seq_len > ($median_75 + (0.1 * $median)))) {
			    # do not include length outliers in the multifasta file
			} else {
			    $nr_allele_num++;
			    print OUTFILE ">$genome_tag\n";
			    my $pos;
			    my $tmp_seq_len = $seq_len;
			    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
				print OUTFILE substr($sequence, $pos, 60), "\n";
				$tmp_seq_len -= 60;
			    }
			    print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
			}
		    }
		}
		$index++;
	    }
	    push (@mf_files, "$multifastadir/$edge_id.fasta\t$nr_allele_num\t$empty");
	    close (OUTFILE);
	} else {
	    my $index = 0;
	    foreach my $genome_tag (@genome_array) {
		if (defined $genome_seqs[$index]) {
		    my $sequence = $genome_seqs[$index];
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
		    }
		}
		$index++;
	    }
	}
	if ($no_stats) {
	    if ($gene_count > 1) {
		if ($cluster1 <= $cluster2) { # only count edge in one direction
		    $total_edge_pgg++;
		    $total_edge_alle_pgg += $gene_count;
		}
	    }
	    next; # do not generate the anomalies file or multiple sequence alignments
	}
	if ($cluster1 <= $cluster2) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
	    if ($gene_count > 1) {
		if ($cluster1 <= $cluster2) { # only count edge in one direction
		    $total_edge_pgg++;
		    $total_edge_alle_pgg += $gene_count;
		}
	    }
	    if ($compute_all || ($target_id ne "")) {
		if ($compute_all && ($gene_count > 0)) {
		    $index = 0;
		    foreach my $genome_tag (@genome_array) {
			if (($gene_count == 1) && ($single_genome eq $genome_tag)) {
			    $uniq_edge{$genome_tag}++;
			} elsif (defined $genome_seqs[$index]) {
			    my $feat_name = $genome_tag . $edge_id;
			    if (!defined $feat_pres{$genome_seqs[$index]}) {
				print STDERR "$genome_tag:$index:$feat_name:\n$genome_seqs[$index]\n" if ($DEBUG);
			    }
			    my $seq_len = $edge_hash{$feat_name}->{'len'};
			    my $ambig_count = $genome_seqs[$index] =~ tr/ACGTacgt//c;
			    if (($genome_seqs[$index] ne "EMPTY") && (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20))) {
				# do not include sequences with a lot of ambiguous base calls
				$gapped_edge{$genome_tag}++;
				print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tgapped_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    } else {
				if (defined $edge_hash{$feat_name}->{'gapped'}) {
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tgapped_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				}
				if ($feat_pres{$genome_seqs[$index]} == 1) {
				    if (((100 * $gene_count) / $genome_number) >= 75) {
					$uniq_edge_alle_75_100{$genome_tag}++;
				    } elsif (((100 * $gene_count) / $genome_number) <= 25) {
					$uniq_edge_alle_0_25{$genome_tag}++;
				    } else {
					$uniq_edge_alle_25_75{$genome_tag}++;
				    }
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tuniq_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				} else {
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tidentical_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				}
				if ($seq_len < 0) {
				    $seq_len = 0;
				}
				#print STDERR "$genome_tag($mean):$median_25:$median:$median_75:$seq_len\n" if ($DEBUG);
				if (($seq_len < $min) && ($seq_len < ($median - (0.02 * $median)))) {
				    $very_short_edge{$genome_tag}++;
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tvery_short_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				    #print STDERR "VERY SHORT\n" if ($DEBUG);
				} elsif ($seq_len < ($median_25 - (0.1 * $median))) {
				    $short_edge{$genome_tag}++;
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tshort_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				    #print STDERR "SHORT\n" if ($DEBUG);
				} elsif (($seq_len > $max) && ($seq_len > ($median + (0.02 * $median)))) {
				    $very_long_edge{$genome_tag}++;
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tvery_long_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				    #print STDERR "VERY LONG\n" if ($DEBUG);
				} elsif ($seq_len > ($median_75 + (0.1 * $median))) {
				    $long_edge{$genome_tag}++;
				    print DIFFFILE "$genome_tag\t$edge_hash{$feat_name}->{'contig'}\tlong_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
				    #print STDERR "LONG\n" if ($DEBUG);
				}
			    }
			    $total_edge{$genome_tag}++;
			} else {
			    if (((100 * $gene_count) / $genome_number) >= 75) {
				if ((defined $single_copy_core{$cluster1}) && (defined $single_copy_core{$cluster2})){
				    $miss_sing_edge{$genome_tag}++;
				} else {
				    $missing_75e{$genome_tag}++;
				}
			    }
			}
			$index++;
		    }
		} else {
		    if (($gene_count == 0) && ($target_sequence ne "")) {
			$uniq_edge{$target_id}++;
		    } elsif ($target_sequence ne "") {
			my $feat_name = $target_id . $edge_id;
			my $seq_len = $edge_hash{$feat_name}->{'len'};
			my $ambig_count = $target_sequence =~ tr/ACGTacgt//c;
			if (($target_sequence ne "EMPTY") && (($ambig_count >= 10) || ((($ambig_count * 100) / $seq_len) > 20))) {
			    # do not include sequences with a lot of ambiguous base calls
			    $gapped_edge{$target_id}++;
			    print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tgapped_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			} else {
			    if (defined $edge_hash{$feat_name}->{'gapped'}) {
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tgapped_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    }
			    if ($feat_pres{$target_sequence} == 1) {
				if (((100 * $gene_count) / $genome_number) >= 75) {
				    $uniq_edge_alle_75_100{$target_id}++;
				} elsif (((100 * $gene_count) / $genome_number) <= 25) {
				    $uniq_edge_alle_0_25{$target_id}++;
				} else {
				    $uniq_edge_alle_25_75{$target_id}++;
				}
				if ($align_new) {
				    print STDERR "$target_id\t$edge_id\t$seq_len\t$median_25\t$median\t$median_75\n" if ($DEBUG);
				    if (($seq_len >= ($median_25 - (0.1 * $median))) && ($seq_len <= ($median_75 + (0.1 * $median)))) {
					my $mf_file = "$multifastadir/$edge_id.fasta";
					my $msa_file = "$multifastadir/$edge_id.afa";
					my $stats_file = "$multifastadir/$edge_id.stats";
					my $empty_file = "$multifastadir/$edge_id.empty";
					my $max_all;
					my $median_target;
					my $max_cols_all;
					my $cols_target;
					print STDERR "Good length\t$stats_file\t$msa_file\t$mf_file\n" if ($DEBUG);
					if ((-s $stats_file) && ((-s $msa_file) || (-s $mf_file))){
					    unless (open (STATSFILE, "<", $stats_file) )  {
						die ("ERROR: cannot open file $stats_file\n");
					    }
					    while (my $line = <STATSFILE>)  {
						chomp $line;
						(my @fields) = split(/\t/, $line);  # split the scalar $line on tab
						if ($fields[0] eq "All") {
						    $max_all = $fields[4];
						} elsif ($fields[0] eq "UniqueAlleleCount") {
						    $max_cols_all = $fields[4];
						}
					    }
					    close(STATSFILE);
					    my $seq_file = "$multifastadir/edge_TMP_" . $$ . "_$target_id.$edge_id.fasta";
					    my $combined_file = "$multifastadir/edge_TMP_" . $$ . "_$target_id.$edge_id.afa";
					    my $com_stats_file = "$multifastadir/edge_TMP_" . $$ . "_$target_id.$edge_id.stats";
					    unless (open (OUTFILE, ">", $seq_file) )  {
						die ("ERROR: cannot open file $seq_file\n");
					    }
					    print OUTFILE ">$target_id\n";
					    my $pos;
					    my $tmp_seq_len = $seq_len;
					    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
						print OUTFILE substr($target_sequence, $pos, 60), "\n";
						$tmp_seq_len -= 60;
					    }
					    print OUTFILE substr($target_sequence, $pos, $tmp_seq_len), "\n";
					    close (OUTFILE);
					    if (-s $msa_file) {
						my $muscle_args = " -profile -in1 $msa_file -in2 $seq_file -out $combined_file -diags -quiet -verbose";
						my $muscle_exec = $muscle_path . $muscle_args;
						`/usr/bin/time -o tmp_cpu_stats -v $muscle_exec`;
						`echo "***$muscle_exec***length($seq_len)" >> $cpu_name`;
						`cat tmp_cpu_stats >> $cpu_name`;
						`rm tmp_cpu_stats`;
						&bash_error_check($muscle_exec, $?, $!);
					    } elsif (-s $mf_file) {
						`cat $mf_file >> $seq_file`;
						my $muscle_args = " -in $seq_file -out $combined_file -diags -quiet -verbose";
						my $muscle_exec = $muscle_path . $muscle_args;
						`/usr/bin/time -o tmp_cpu_stats -v $muscle_exec`;
						`echo "***$muscle_exec***length($seq_len)" >> $cpu_name`;
						`cat tmp_cpu_stats >> $cpu_name`;
						`rm tmp_cpu_stats`;
						&bash_error_check($muscle_exec, $?, $!);
					    }
					    my $stats_path = "$bin_directory/summarize_alignment.R";
					    my $stats_args = " $combined_file $target_id";
					    my $stats_exec = $rscript_path . " " . $stats_path . $stats_args;
					    `/usr/bin/time -o tmp_cpu_stats -v $stats_exec > $com_stats_file`;
					    `echo "***$stats_exec***" >> $cpu_name`;
					    `cat tmp_cpu_stats >> $cpu_name`;
					    `rm tmp_cpu_stats`;
					    &bash_error_check("$stats_exec > $com_stats_file", $?, $!);
					    unless (open (STATSFILE, "<", $com_stats_file) )  {
						die ("ERROR: cannot open file $com_stats_file\n");
					    }
					    while (my $line = <STATSFILE>)  {
						print STDERR $line if ($DEBUG);
						chomp $line;
						(my @fields) = split(/\t/, $line);  # split the scalar $line on tab
						if ($fields[0] eq $target_id) {
						    $median_target = $fields[2];
						    $cols_target = $fields[27];
						}
					    }
					    close(STATSFILE);
					    if ($keep_divergent_alignments) {
						`mv $combined_file $keep_divergent_alignments/edge_$target_id.$edge_id.afa`;
						`rm $seq_file $com_stats_file`;
					    } else {
						`rm $seq_file $combined_file $com_stats_file`;
					    }
					    print STDERR "$target_id\t$edge_id\t$median_target\t$max_all\t$cols_target\t$max_cols_all\n" if ($DEBUG);
					    if ($median_target > $max_all) {
						print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tdivergent_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
						$distant_edge_alle{$target_id}++;
					    }
					    if ($cols_target > $max_cols_all) {
						print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tdivergent_columns_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
						$column_edge_alle{$target_id}++;
					    }
					}
				    }
				}
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tuniq_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    } else {
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tidentical_edge\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    }
			    if ($seq_len < 0) {
				$seq_len = 0;
			    }
			    if (($seq_len < $min) && ($seq_len < ($median - (0.02 * $median)))) {
				$very_short_edge{$target_id}++;
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tvery_short_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    } elsif ($seq_len < ($median_25 - (0.1 * $median))) {
				$short_edge{$target_id}++;
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tshort_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    } elsif (($seq_len > $max) && ($seq_len > ($median + (0.02 * $median)))) {
				$very_long_edge{$target_id}++;
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tvery_long_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    } elsif ($seq_len > ($median_75 + (0.1 * $median))) {
				$long_edge{$target_id}++;
				print DIFFFILE "$target_id\t$edge_hash{$feat_name}->{'contig'}\tlong_edge_allele\t$edge_hash{$feat_name}->{'5p'}\t$edge_hash{$feat_name}->{'3p'}\t$edge_hash{$feat_name}->{'len'}\t$edge_name\n";
			    }
			}
			$total_edge{$target_id}++;
		    } else {
			if (((100 * $gene_count) / $genome_number) >= 75) {
			    if ((defined $single_copy_core{$cluster1}) && (defined $single_copy_core{$cluster2})){
				$miss_sing_edge{$target_id}++;
			    } else {
				$missing_75e{$target_id}++;
			    }
			}
		    }
		}
	    }
	}
    }
    if ($remake_files) {
	close (OUTPGGFILE);
	close(SIZEFILE);
    }
    close (PGGFILE);
    return;
}

sub process_medoids {  # read in the medoids for the PGG

    my $medoidfile;
    unless (open ($medoidfile, "<", $medoids_path) )  {
	die ("cannot open medoid fasta file: $medoids_path!\n");
    }
    my $medoidout;
    unless (open ($medoidout, ">", "$basedir/medoids.fasta") )  {
	die ("cannot open output medoid fasta file: $basedir/medoids.fasta!\n");
    }
    my $save_input_separator = $/;
    my $line;
    my $cluster_num = 1;
    $/="\n>";
    while ($line = <$medoidfile>) {
	chomp $line;
	(my $title, my $sequence) = split(/\n/, $line, 2); # split the header line and sequence
	$title =~ s/^>//; # remove leading ">" for first medoid
	$sequence =~ s/\n$//; # remove trailing newline for last medoid
	if ($title =~ /(centroid_)(\d+)(.*)/) {
	    if ($2 != $cluster_num) {
	    die ("ERROR: Bad medoid cluster number in header expecting $cluster_num but found $2\n>$title\nin file $medoids_path!\n");
	    }
	    $title = ">$1$renumber[$cluster_num]$3\n";
	} else {
	    die ("ERROR: Bad medoid header formatting\n>$title\nin file $medoids_path!\n");
	}
	if ($renumber[$cluster_num]) {
	    print $medoidout $title, $sequence, "\n";
	}
	$cluster_num++;
    }
    $/ = $save_input_separator; # restore the input separator
    close ($medoidfile);
    close ($medoidout);
    return;
}

sub read_single_cores                                               # Read in the set of clusters which are single copy core, renumber and output
{
    unless (open(CLUSTER_CORES, "<", $single_cores)) {
	die ("cannot open cores file: $single_cores!\n");
    }
    while (my $line = <CLUSTER_CORES>) {
	chomp $line;
	$line =~ s/\s*//g; # remove all whitespace characters
	$line =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	$single_copy_core{$line} = 1;
	#print STDERR "$line:$single_copy_core{$line}\n" if ($DEBUG);
    }
    close(CLUSTER_CORES);
}

sub write_single_cores                                               # Read in the set of clusters which are single copy core, renumber and output
{
    unless (open(OUT_CORES, ">", "$basedir/single_copy_clusters.txt")) {
	die ("cannot open cores file: $basedir/single_copy_clusters.txt!\n");
    }
    foreach my $clus (sort {$a <=> $b} (keys %single_copy_core)) {
	if ($renumber[$clus] != 0) {
	    print OUT_CORES "$renumber[$clus]\n";
	} else {
	    print STDERR "WARNING: single copy core cluster $clus was not annotated in any genome and hence reduced away!\n";
	}
	#print STDERR "$clus:$renumber[$single_copy_core{$clus}]\n" if ($DEBUG);
    }
    close(OUT_CORES);
}

sub bash_error_check {
    my ($command, $error, $message) = @_;
    if (!$error) {
	return(0);
    }
    print STDERR "$command FAILED\n";
    if ($error == -1) {
	printf STDERR "failed to execute code(%d): %s\n", $error >> 8, $message;
    } elsif ($error & 127) {
	printf STDERR "child died with code %d signal %d, %s coredump\n", $error >> 8, ($error & 127),  ($error & 128) ? 'with' : 'without';
    } else {
	printf STDERR "child exited with value %d\n", $error >> 8;
    }
    return(1);
}

sub launch_grid_job {
# Given a shell script, launch it via qsub.

    my ( $name, $project_code, $working_dir, $shell_script, $stdoutdir, $stderrdir, $queue, $job_array_max ) = @_;

    my $qsub_command = "qsub -V -o $stdoutdir -e $stderrdir -r n -N $name";
    if ($queue eq "NONE") {
	$qsub_command .= " -d $working_dir";
    } else {
	$qsub_command .= " -wd $working_dir";
	$qsub_command .= " -terse";
    }
    $qsub_command .= " -P $project_code" if ($project_code && ($project_code ne "NONE"));
    $qsub_command .= " -l $queue" if ($queue && ($queue ne "NONE"));
    $qsub_command .= " -t 1-$job_array_max" if $job_array_max;

    #$qsub_command .= " $shell_script";
    $qsub_job_num++;
    my $qsub_exec = $cwd . "/TMP_" . $qsub_job_num . "_" . $name;
    unless (open(OUT_QSUB, ">", $qsub_exec)) {
	die ("cannot open qsub executable file $qsub_exec!\n");
    }
    if (substr($shell_script, 0, 1) ne "/") {
	$shell_script = $cwd . "/$shell_script";
    }
    print OUT_QSUB $shell_script;
    close(OUT_QSUB);
    `chmod +x $qsub_exec`;
    $qsub_command .= " $qsub_exec";

    my $job_id = `$qsub_command`;
    $job_id =~ s/\s*//g; # remove all whitespace characters

    if (&bash_error_check($qsub_command, $?, $!)) {
        die "Problem submitting the job!: $job_id\n$qsub_command\n$shell_script\n$qsub_exec\n";
    }

    return $job_id;

}


sub wait_for_grid_jobs {
    # Given a hash of job ids wait until hash is reduced to number of jobs specified and return number of jobs; name is the job name
    
    my ( $queue, $name, $number, $job_ids ) = @_;
    my $size = scalar( keys %{$job_ids} );

    if ($queue eq "NONE") {
	sleep 300; # need to wait to make sure qstat knows about all submitted jobs
    }
    while ( $size > $number ) {
	sleep 60;
	if ($queue eq "NONE") {
	    my $response = `qstat 2>&1`;
	    &parse_response_qstat( $response, $name, $job_ids );
	} else {
	    my $response = `qacct -j $name 2>&1`;
	    &parse_response_qacct( $response, $job_ids );
	}
	$size = scalar( keys %{$job_ids} );
    }
    return ($size);
}

sub parse_response_qstat {
# NOT INTENDED TO BE CALLED DIRECTLY.
# Given a qstat response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $job_name, $job_ids ) = @_;
    my @qstat_array = split ( /\n/, $response );
    my %running = ();
    foreach my $line (@qstat_array) {
	my @fields = split ( /\s+/, $line );
	if (($fields[1] eq $job_name) && ($fields[4] ne "C")) {
	    $running{$fields[0]} = $fields[4];
	}
    }
    foreach my $job_id (keys %{ $job_ids }) {
	if (! ( defined $running{$job_id} )) {
	    delete ( $job_ids->{$job_id} )
	}
    }
    return;
}

sub parse_response_qacct {
# NOT INTENDED TO BE CALLED DIRECTLY.
# Given a qacct response, delete a job id from the loop-control-hash when
# a statisfactory state is seen.

    my ( $response, $job_ids ) = @_;
    my @qacct_array = split ( /=+\n/, $response );
    if (scalar(@qacct_array) <= 1) {
	return; # jobs haven't hit the grid yet
    }
    shift @qacct_array; # get rid of empty record at beginning.

    for my $record ( @qacct_array ) {

        next if ( $record =~ /error: ignoring invalid entry in line/ );

        chomp $record;

        my @rec_array = split ( "\n", $record );

        for my $line (@rec_array) {

            $line =~ s/(.*\S)\s+$/$1/;
            my ( $key, $value ) = split ( /\s+/, $line, 2 );
	    if ($key eq "jobnumber") {
		if ( defined $job_ids->{$value} ) {
		    delete ( $job_ids->{$value} )
		}
	    }
	}
    }
    return;
}

sub compute_alignments
{
    my $job_name = "pggmus_" . $$; #use a common job name so that qacct can access all of them together
    my %job_ids = ();
    my $num_jobs = 0;
    my $total_jobs = 0;
    my @msa_files = ();
    foreach my $line (@mf_files) {
	(my $mf_file, my $num_alleles, my $empty) = split(/\t/, $line);  # split the scalar $line on tab
	my $msa_file = $mf_file;
	my $empty_file;
	my $stats_file;
	$msa_file =~ s/\.fasta$//;
	$empty_file = $stats_file = $msa_file;
	$msa_file .= ".afa";
	$empty_file .= ".empty";
	$stats_file .= ".stats";
	print STDERR "$mf_file:$num_alleles:$empty:$msa_file:$empty_file:$stats_file\n" if ($DEBUG);
	if ((-e $msa_file) && (-s $msa_file)) {
	    next; # skip if file already exists and is not zero size
	}
	if ($empty) {
	    unless (open(OUT_EMPTY, ">", "$empty_file")) {
		die ("cannot open empty indicator file: $empty_file!\n");
	    }
	    close(OUT_EMPTY); # creating a zero length indicator file that this edge contained a zero length entry
	}
	if ($num_alleles <= 1) {
	    unless (open(OUT_STATS, ">", "$stats_file")) {
		die ("cannot open alignment stats file: $stats_file!\n");
	    }
	    print OUT_STATS "\tMean\tMedian\tMin\tMax\tSD\t0%\t5%\t10%\t15%\t20%\t25%\t30%\t35%\t40%\t45%\t50%\t55%\t60%\t65%\t70%\t75%\t80%\t85%\t90%\t95%\t100%\nAll\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\nUniqueAlleleCount\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"; #print baseline stats when no alignment is available (mean, median, maximmum, std. dev. of the normalized edit distance - percent identity difference)
	    close(OUT_STATS);
	    next; #nothing to align
	}
	my $identifier = basename($mf_file);
	$identifier =~ s/\.fasta$//;
	my $muscle_args = " -in $mf_file -out $msa_file -diags -quiet -verbose";
	my $muscle_exec = $muscle_path . $muscle_args;
	my $stdoutfile = $cwd . "/" . $identifier . "_stdout";
	my $stderrfile = $cwd . "/" . $identifier . "_stderr";
	my $working_dir = $cwd;
	print STDERR "Launched $muscle_exec\n" if ($DEBUG);
	$job_ids{&launch_grid_job($job_name, $project, $working_dir, $muscle_exec, $stdoutfile, $stderrfile, $qsub_queue)} = 1;
	$num_jobs++;
	$total_jobs++;
	if ($num_jobs >= $max_grid_jobs) {
	    $num_jobs = &wait_for_grid_jobs($qsub_queue, $job_name, ((($max_grid_jobs - 10) > 0) ? ($max_grid_jobs - 10) : 0), \%job_ids);
	}
    }
    &wait_for_grid_jobs($qsub_queue, $job_name, 0, \%job_ids);
    $num_jobs = 0;
    foreach my $line (@mf_files) {
	(my $mf_file, my $num_alleles, my $empty) = split(/\t/, $line);  # split the scalar $line on tab
	my $msa_file = $mf_file;
	$msa_file =~ s/\.fasta$//;
	$msa_file .= ".afa";
	if ($num_alleles <= 1) {
	    next; #nothing to align
	}
	if ((-e $msa_file) && (-s $msa_file)) {
	    push (@msa_files, $msa_file);
	    next; # skip if file already exists and is not zero size
	}
	$num_jobs++;
    }
    if ($num_jobs > ($total_jobs / 10)) {
	die "Too many Muscle grid jobs failed $num_jobs\n";
    } elsif ($num_jobs > 0) {
	for (my $k=0; $k <= 2; $k++){ #try a maximum of 3 times on failed jobs
	    print STDERR "Resubmit $num_jobs jobs Iteration $k\n" if ($DEBUG);
	    %job_ids = ();
	    $num_jobs = 0;
	    @msa_files = ();
	    foreach my $line (@mf_files) {
		(my $mf_file, my $num_alleles, my $empty) = split(/\t/, $line);  # split the scalar $line on tab
		my $msa_file = $mf_file;
		$msa_file =~ s/\.fasta$//;
		$msa_file .= ".afa";
		if ($num_alleles <= 1) {
		    next; #nothing to align
		}
		if ((-e $msa_file) && (-s $msa_file)) {
		    push (@msa_files, $msa_file);
		    next; # skip if file already exists and is not zero size
		}
		my $identifier = basename($mf_file);
		$identifier =~ s/\.fasta$//;
		my $muscle_args = " -in $mf_file -out $msa_file -diags -quiet -verbose";
		my $muscle_exec = $muscle_path . $muscle_args;
		my $stdoutfile = $cwd . "/" . $identifier . "_stdout";
		my $stderrfile = $cwd . "/" . $identifier . "_stderr";
		my $working_dir = $cwd;
		print STDERR "Relaunched $muscle_exec\n" if ($DEBUG);
		$job_ids{&launch_grid_job($job_name, $project, $working_dir, $muscle_exec, $stdoutfile, $stderrfile, $qsub_queue)} = 1;
		$num_jobs++;
	    }
	    if ($num_jobs == 0) {
		last; # no failed jobs
	    }
	    print STDERR "Iteration $k: $num_jobs Muscle jobs relaunched\n" if ($DEBUG);
	    &wait_for_grid_jobs($qsub_queue, $job_name, 0, \%job_ids);
	}
    }
    print STDERR "$num_jobs failed Muscle jobs\n";
    my $stats_path = "$bin_directory/summarize_alignment.R";
    $job_name = "pggstat_" . $$; #use a common job name so that qacct can access all of them together
    %job_ids = ();
    $num_jobs = 0;
    foreach my $msa_file (@msa_files) {
	my $stats_file = $msa_file;
	$stats_file =~ s/\.afa$//;
	$stats_file .= ".stats";
	print STDERR "$msa_file:$stats_file\n" if ($DEBUG);
	if ((-e $stats_file) && (-s $stats_file)) {
	    next; # skip if file already exists and is not zero size
	}
	my $identifier = basename($msa_file);
	$identifier =~ s/\.afa$//;
	my $stats_args = " $msa_file";
	my $stats_exec = $rscript_path . " " . $stats_path . $stats_args;
	my $stdoutfile = $stats_file;
	my $stderrfile = $cwd . "/" . $identifier . "_stderr";
	my $working_dir = $cwd;
	print STDERR "Launched $stats_exec\n" if ($DEBUG);
	$job_ids{&launch_grid_job($job_name, $project, $working_dir, $stats_exec, $stdoutfile, $stderrfile, $qsub_queue)} = 1;
	$num_jobs++;
	if ($num_jobs >= $max_grid_jobs) {
	    $num_jobs = &wait_for_grid_jobs($qsub_queue, $job_name, ((($max_grid_jobs - 10) > 0) ? ($max_grid_jobs - 10) : 0), \%job_ids);
	}
    }
    &wait_for_grid_jobs($qsub_queue, $job_name, 0, \%job_ids);
    $num_jobs = 0;
    foreach my $msa_file (@msa_files) {
	my $stats_file = $msa_file;
	$stats_file =~ s/\.afa$//;
	$stats_file .= ".stats";
	if ((-e $stats_file) && (-s $stats_file)) {
	    next; # skip if file already exists and is not zero size
	}
	$num_jobs++;
    }
    if ($num_jobs > ($max_grid_jobs / 2)) {
	die "Too many Alignment Stats grid jobs failed $num_jobs\n";
    } elsif ($num_jobs > 0) {
	for (my $k=0; $k <= 2; $k++){ #try a maximum of 3 times on failed jobs
	    print STDERR "Resubmit $num_jobs jobs Iteration $k\n" if ($DEBUG);
	    %job_ids = ();
	    $num_jobs = 0;
	    foreach my $msa_file (@msa_files) {
		my $stats_file = $msa_file;
		$stats_file =~ s/\.afa$//;
		$stats_file .= ".stats";
		if ((-e $stats_file) && (-s $stats_file)) {
		    next; # skip if file already exists and is not zero size
		}
		my $identifier = basename($msa_file);
		$identifier =~ s/\.afa$//;
		my $stats_args = " $msa_file";
		my $stats_exec = $rscript_path . " " . $stats_path . $stats_args;
		my $stdoutfile = $stats_file;
		my $stderrfile = $cwd . "/" . $identifier . "_stderr";
		my $working_dir = $cwd;
		print STDERR "Relaunched $stats_exec\n";
		$job_ids{&launch_grid_job($job_name, $project, $working_dir, $stats_exec, $stdoutfile, $stderrfile, $qsub_queue)} = 1;
		$num_jobs++;
	    }
	    if ($num_jobs == 0) {
		last; # no failed jobs
	    }
	    print STDERR "Iteration $k: $num_jobs Alignment Stats jobs relaunched\n" if ($DEBUG);
	    &wait_for_grid_jobs($qsub_queue, $job_name, 0, \%job_ids);
	}
    }
    print STDERR "$num_jobs failed Summarize alignment jobs\n";
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
Example: $prog -b output_dir -g panoct_genome_tag_file -m panoct_matchtable_file -a panoct_attribute_file -p panoct_pan-genome_graph_file -T topology_file
Version: $version
 Option:
     -h: print this help page
     -b: directory path for where to put output multifasta files[DEFAULT = PWD]
     -B: directory path for where to put other output files[DEFAULT = PWD]
     -N: directory path for where to find executables for pan-genome tools
     -g: two colums tab delimited file (input): col1 genome_tag used in panoct in same order as panoct; col2 file name for the fasta contig file for this genome
     -m: panoct matchtable file (input)
     -a: combined panoct attribute file (input)
     -p: a panoct pan_genome_graph_file
     -T: a topology file (3 columns tab delimited: genome id, contig id, circular/linear
     -t: target genome id if stats only wanted for one genome
     -A: output stats for all genomes
     -L: generate multiple sequence alignment files in the multifasta output directory - cannot be combined with -S
     -l: generate tempory incremental multiple sequence alignment files in the multifasta output directory to compare new sequences
     -k: keep mutliple sequence alignments from the target genome which are divergent from the PGG (argument is the directory to keep them in)
     -S: suppress multifasta output
     -M: medoids input file name - output will be in -B directory (medoids.fasta)
     -c: output size files in -B directory computed cluster sizes(program  adds cluster_sizes.txt) and edge sizes(program  adds edge_sizes.txt) files
     -P: project code for qsub jobs
     -j: the maximum number of grid jobs to run for multiple sequence alignments
     -f: generate multifasta files with all alleles for clusters and edges
     -F: use multifasta files with all alleles for clusters and edges instead of extracting them directly from genome fasta files
     -C: path to the Muslce executable for multiple sequence alignments - default /usr/local/bin/muscle
     -r: path to the Rscript executable for Rscript scripts - default /usr/local/bin/Rscript
     -Q: name of qsub queue to use for Muscle jobs - default himem
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
If -A or -t are used a cluster_stats.txt file is created with a header line and tab delimitted cluster and edge stats.
_EOB_
    exit;
}

########################################  M A I N  #####################################################
if ($no_stats) {
    &get_genome_names;
    print STDERR "No stats just remaking files so not reading genomes\n";
} else {
    print  STDERR "Getting topology from $topology_file\n";
    &read_topology;
    print  STDERR "Gathering gene coordinates and attribute information from $att_file\n";
    &get_attributes;
    if ($write_multifasta) {
	print STDERR "Reading genomes, matchtable, and pgg then writing multifasta clusters and edges\n";
	&output_multifasta;
    } else {
	print  STDERR "Getting genome names and contig sequences from $genomes_file_name\n";
	&get_genomes;
    }
}
print  STDERR "Reading single copy core clusters from $single_cores\n";
&read_single_cores;
foreach my $genome_tag (@genome_array) {
    $total_clus{$genome_tag} = 0;           # key = genome ID, value = number of nonsingleton clusters
    $total_edge{$genome_tag} = 0;           # key = genome ID, value = number of nonsingleton edges
    $uniq_clus{$genome_tag} = 0;            # key = genome ID, value = number of unique/singleton clusters
    $uniq_edge{$genome_tag} = 0;            # key = genome ID, value = number of unique/singleton edges
    $short_clus{$genome_tag} = 0;           # key = genome ID, value = number of short clusters
    $short_edge{$genome_tag} = 0;           # key = genome ID, value = number of short edges
    $long_clus{$genome_tag} = 0;            # key = genome ID, value = number of long clusters
    $long_edge{$genome_tag} = 0;            # key = genome ID, value = number of long edges
    $very_short_clus{$genome_tag} = 0;      # key = genome ID, value = number of very short clusters
    $very_short_edge{$genome_tag} = 0;      # key = genome ID, value = number of very short edges
    $very_long_clus{$genome_tag} = 0;       # key = genome ID, value = number of very long clusters
    $very_long_edge{$genome_tag} = 0;       # key = genome ID, value = number of very long edges
    $frameshift{$genome_tag} = 0;           # key = genome ID, value = number of probable frameshifts
    $miss_sing_core{$genome_tag} = 0;       # key = genome ID, value = number of clusters missing  for single copy core clusters
    $miss_sing_edge{$genome_tag} = 0;       # key = genome ID, value = number of edges missing  for edges between single copy core clusters
    $missing_75c{$genome_tag} = 0;          # key = genome ID, value = number of clusters missing  for clusters in 75-100% of genomes
    $missing_75e{$genome_tag} = 0;          # key = genome ID, value = number of edges missing  for edges in 75-100% of genomes
    $uniq_clus_alle_75_100{$genome_tag} = 0;# key = genome ID, value = number of unique alleles for clusters in 75-100% of genomes
    $uniq_clus_alle_25_75{$genome_tag} = 0; # key = genome ID, value = number of unique alleles for clusters in 25-75% of genomes
    $uniq_clus_alle_0_25{$genome_tag} = 0;  # key = genome ID, value = number of unique alleles for clusters in 0-25% of genomes
    $uniq_edge_alle_75_100{$genome_tag} = 0;# key = genome ID, value = number of unique alleles for edges in 75-100% of genomes
    $uniq_edge_alle_25_75{$genome_tag} = 0; # key = genome ID, value = number of unique alleles for edges in 25-75% of genomes
    $uniq_edge_alle_0_25{$genome_tag} = 0;  # key = genome ID, value = number of unique alleles for edges in 0-25% of genomes
    $distant_clus_alle{$genome_tag} = 0;    # key = genome ID, value = number of distant cluster alleles
    $distant_edge_alle{$genome_tag} = 0;    # key = genome ID, value = number of distant edge alleles
    $column_clus_alle{$genome_tag} = 0;     # key = genome ID, value = number of distant (by unique column characters) cluster alleles
    $column_edge_alle{$genome_tag} = 0;     # key = genome ID, value = number of distant (by unique column characters) edge alleles
    $gapped_clus{$genome_tag} = 0;          # key = genome ID, value = number of gapped (really too many ambiguous characters) clusters
    $gapped_edge{$genome_tag} = 0;          # key = genome ID, value = number of gapped (really too many ambiguous characters) edges
}    
print  STDERR "Reading matchtable from $matchtable_file\n";
if (!$suppress) {
    print  STDERR " and outputting cluster multifasta files to $multifastadir\n";
}
if (($compute_all || ($target_id ne "")) && (!$no_stats)) {
    unless (open (DIFFFILE, ">", "$basedir/anomalies.txt") ) {
	die ("ERROR: cannot open file $basedir/anomalies.txt\n");
    }
    print DIFFFILE "Genome ID\tContig ID\tType\tStart\tEnd\tLength\tCluster-Edge ID\n";
}
&process_matchtable;
print  STDERR "Reading pan-genome graph from $pgg_file\n";
if (!$suppress) {
    print  STDERR " and outputting edge multifasta files to $multifastadir\n";
}
&process_pgg;
unless (open (SIZESFILE, ">", "$basedir/ce_sizes.txt") ) {
    die ("ERROR: cannot open file $basedir/ce_sizes.txt\n");
}
print SIZESFILE "Type\tCore\tShared\tSingle\tReduced\tShared+Core\tS+C Alleles\nCluster\t$num_core_clus\t$num_shared_clus\t$num_size_one_clus\t$num_reduced_clus\t$total_clus_pgg\t$total_clus_alle_pgg\nEdge\t$num_core_edge\t$num_shared_edge\t$num_size_one_edge\t$num_reduced_edge\t$total_edge_pgg\t$total_edge_alle_pgg\n";
close (SIZESFILE);
if (($compute_all || ($target_id ne "")) && (!$no_stats)) {
    close (DIFFFILE);
    unless (open (STATFILE, ">", "$basedir/cluster_stats.txt") ) {
	die ("ERROR: cannot open file $basedir/cluster_stats.txt\n");
    }
    print  STDERR "Writing feature stats to $basedir/cluster_stats.txt\n";
    print STATFILE "Genome\tuniq_clus\tuniq_clus_alle_0_25\tuniq_clus_alle_25_75\tuniq_clus_alle_75_100\tuniq_edge\tuniq_edge_alle_0_25\tuniq_edge_alle_25_75\tuniq_edge_alle_75_100\tshort_clus\tlong_clus\tshort_edge\tlong_edge\tvery_short_clus\tvery_long_clus\tvery_short_edge\tvery_long_edge\tframeshift\tmissing_75clus\tmissing_75edge\tmiss_sing_core\tmiss_sing_edge\tdist_clus_alle\tcol_clus_alle\tdist_edge_alle\tcol_edge_alle\tgapped_clus\tgapped_edge\ttotal_clus($total_clus_pgg,$total_clus_alle_pgg)\ttotal_edge($total_edge_pgg,$total_edge_alle_pgg)\n";
    if ($compute_all) {
	foreach my $genome_tag (@genome_array) {
	    print STATFILE "$genome_tag\t$uniq_clus{$genome_tag}\t$uniq_clus_alle_0_25{$genome_tag}\t$uniq_clus_alle_25_75{$genome_tag}\t$uniq_clus_alle_75_100{$genome_tag}\t$uniq_edge{$genome_tag}\t$uniq_edge_alle_0_25{$genome_tag}\t$uniq_edge_alle_25_75{$genome_tag}\t$uniq_edge_alle_75_100{$genome_tag}\t$short_clus{$genome_tag}\t$long_clus{$genome_tag}\t$short_edge{$genome_tag}\t$long_edge{$genome_tag}\t$very_short_clus{$genome_tag}\t$very_long_clus{$genome_tag}\t$very_short_edge{$genome_tag}\t$very_long_edge{$genome_tag}\t$frameshift{$genome_tag}\t$missing_75c{$genome_tag}\t$missing_75e{$genome_tag}\t$miss_sing_core{$genome_tag}\t$miss_sing_edge{$genome_tag}\t$distant_clus_alle{$genome_tag}\t$column_clus_alle{$genome_tag}\t$distant_edge_alle{$genome_tag}\t$column_edge_alle{$genome_tag}\t$gapped_clus{$genome_tag}\t$gapped_edge{$genome_tag}\t$total_clus{$genome_tag}\t$total_edge{$genome_tag}\n";
	}
    } else {
	print STATFILE "$target_id\t$uniq_clus{$target_id}\t$uniq_clus_alle_0_25{$target_id}\t$uniq_clus_alle_25_75{$target_id}\t$uniq_clus_alle_75_100{$target_id}\t$uniq_edge{$target_id}\t$uniq_edge_alle_0_25{$target_id}\t$uniq_edge_alle_25_75{$target_id}\t$uniq_edge_alle_75_100{$target_id}\t$short_clus{$target_id}\t$long_clus{$target_id}\t$short_edge{$target_id}\t$long_edge{$target_id}\t$very_short_clus{$target_id}\t$very_long_clus{$target_id}\t$very_short_edge{$target_id}\t$very_long_edge{$target_id}\t$frameshift{$target_id}\t$missing_75c{$target_id}\t$missing_75e{$target_id}\t$miss_sing_core{$target_id}\t$miss_sing_edge{$target_id}\t$distant_clus_alle{$target_id}\t$column_clus_alle{$target_id}\t$distant_edge_alle{$target_id}\t$column_edge_alle{$target_id}\t$gapped_clus{$target_id}\t$gapped_edge{$target_id}\t$total_clus{$target_id}\t$total_edge{$target_id}\n";
    }
    close (STATFILE);
}
print STDERR "Checking for further processing\n";
if ($remake_files) {
    print  STDERR "Remaking medoids and single core files accounting for cluster dropouts\n";
    &process_medoids;
    &write_single_cores;
}
if ((!$suppress) && ($align_all)) {
    print  STDERR "Computing multiple sequence alignments and  statistics\n";
    &compute_alignments;
}
print STDERR "Finished processing - exiting\n";
exit(0);
