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
getopts ('RSADhb:B:m:p:a:g:t:M:s:');
our ($opt_S,$opt_A,$opt_D,$opt_h,$opt_b,$opt_m,$opt_p,$opt_a,$opt_g,$opt_t,$opt_M,$opt_R,$opt_B,$opt_s);

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "ver1.0";
my $basedir;
my $multifastadir;
my $matchtable_file;
my $att_file;
my $pgg_file;
my $DEBUG;
my $genomes_file_name;
my @genome_array = ();
my $genome_number;
my $compute_all = 0;
my $suppress = 0;
my $target_id = "";
my $remake_files;
my $medoids_path;
my $single_cores;
if ($opt_M) {$medoids_path = $opt_M;} else {$medoids_path = "";}
if ($opt_s) {$single_cores = $opt_s;} else {$single_cores = "";}
if ($opt_R) {$remake_files = 1;} else {$remake_files = 0;}
if ($opt_A) {$compute_all = 1;} else {$compute_all = 0;} # flag to compute statistics for all
if ($opt_S) {$suppress = 1;} else {$suppress = 0;} # flag to suppress outputting multifasta files
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_t) {$target_id = $opt_t;} #set target genome id for statistics
if ($opt_b) {$multifastadir = $opt_b;} else { $multifastadir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if ($opt_B) {$basedir = $opt_B;} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory for output files
if (($opt_p) && (-s "$opt_p")) {$pgg_file = $opt_p;} else { print STDERR "Error with -p $opt_p\n"; &option_help; } # if no value for option p (pan-genome graph input file), quit with help menu
if (($opt_m) && (-s "$opt_m")) {$matchtable_file = $opt_m;} else { print STDERR "Error with -m $opt_m\n"; &option_help; } # if no value for option m (matchtable input file), quit with help menu
if (($opt_a) && (-s "$opt_a")) {$att_file = $opt_a;} else { print STDERR "Error with -a $opt_a\n"; &option_help; } # if no value for option a (attribute input file), quit with help menu
if (($opt_g) && (-s "$opt_g")) {$genomes_file_name = $opt_g;} else { print STDERR "Error with -g\n"; &option_help; } # if no value for option g (genome tags and contig file names input file), quit with help menu

my %is_circular = ();          # key1 = genome ID, key2 = contig_name, value = 1 if circular 0 otherwise
my %feat_hash = ();            # Key1 = feat_name Key2 = struct members with their values (5p,3p,anno,gtag)
my %cluster_to_feat_hash = (); # Key1 = genome tag Key2 = cluster_id Value = feat_name
my %genseq_hash = ();          # Key1 = genome tag Key2 = contig_name Value = contig sequence
my %genseq_len = ();           # Key1 = genome tag Key2 = contig_name Value = length of contig sequence
my %uniq_clus = ();            # key = genome ID, value = number of unique/singleton clusters
my %uniq_edge = ();            # key = genome ID, value = number of unique/singleton edges
my %short_clus = ();           # key = genome ID, value = number of short clusters
my %short_edge = ();           # key = genome ID, value = number of short edges
my %long_clus = ();            # key = genome ID, value = number of long clusters
my %long_edge = ();            # key = genome ID, value = number of long edges
my %frameshift = ();           # key = genome ID, value = number of probable frameshifts
my %missing_75c = ();          # key = genome ID, value = number of clusters missing  for clusters in 75-100% of genomes
my %missing_75e = ();          # key = genome ID, value = number of edges missing  for edges in 75-100% of genomes
my %uniq_clus_alle_75_100 = ();# key = genome ID, value = number of unique alleles for clusters in 75-100% of genomes
my %uniq_clus_alle_25_75 = (); # key = genome ID, value = number of unique alleles for clusters in 25-75% of genomes
my %uniq_clus_alle_0_25 = ();  # key = genome ID, value = number of unique alleles for clusters in 0-25% of genomes
my %uniq_edge_alle_75_100 = ();# key = genome ID, value = number of unique alleles for edges in 75-100% of genomes
my %uniq_edge_alle_25_75 = (); # key = genome ID, value = number of unique alleles for edges in 25-75% of genomes
my %uniq_edge_alle_0_25 = ();  # key = genome ID, value = number of unique alleles for edges in 0-25% of genomes
my @renumber = ();             # maps old cluster numbers to new cluster numbers

sub get_genomes {  # obtain list of genomes - must be in the same order as the matchtable columns - and the mulitfasta contigs file for the genomes
   
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$genomes_file_name") || die ("ERROR: cannot open file $genomes_file_name\n");
    print "Order of genomes in $genomes_file_name with array index\n";
    my $target_found = 0;
    while (my $line1 = <$infile>)  {
	chomp $line1;
	(my $name, my $contig_file) = split(/\t/, $line1);  # split the scalar $line on tab

	if (defined $genseq_hash{$name})  {
	    die ("ERROR:  You have more than one occurance of $name in $genomes_file_name!\n");
	} else  {
	    push (@genome_array, $name); # populate the genome_array in the order of the genome file
	    print "$name\t$genome_number\n";
	    $genome_number++;
	    if ($target_id ne "") {
		if ($target_id eq $name) {
		    $target_found = 1;
		}
	    }
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
	    $id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
	    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet characters
	    $genseq_hash{$name}->{$id} = $sequence;
	    $genseq_len{$name}->{$id} = length($sequence);
	    if ($fields[1] eq "circular") {
		$is_circular{$name}->{$id} = 1;
	    } else {
		$is_circular{$name}->{$id} = 0;
	    }
	    $title = ""; # clear the title for the next contig
	    $sequence = ""; #clear out the sequence for the next contig
	}
	$/ = $save_input_separator; # restore the input separator
	close ($contigfile);
    }
    if (($target_id ne "") && !$target_found) {
	die ("ERROR: Did not find target genome: $target_id in genome list file: $genomes_file_name.\n");
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
	die ("ERROR: cannot open file $att_file.\n");
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
	$asmbl_id =~ s/\.\d+$//; # remove trailing version number if it exists - hopefully nonversioned contig names do not have this!
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

    if ($remake_files) {
	unless (open (OUTMATCHFILE, ">", "$basedir/matchtable.txt") )  {
	    die ("ERROR: cannot open file $basedir/matchtable.txt!\n");
	}
	unless (open (SIZEFILE, ">", "$basedir/cluster_sizes.txt") ) {
	    die ("ERROR: cannot open file $basedir/cluster_sizes.txt!\n");
	}
    }
    unless (open (TABLEFILE, "<$matchtable_file") )  {
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
	if (!$suppress) {
	    unless (open (OUTFILE, ">$multifastadir/cluster_$cluster_id.fasta") )  {
		die ("ERROR: cannot open file $multifastadir/cluster_$cluster_id.fasta\n");
	    }
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
	foreach my $feat_name (@feat_names) {
	    $genome_tag = shift @tmp_array;
	    if (($feat_name eq "----------") || ($feat_name eq "")) { #this is a placeholder and can be skipped
		$index++;
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
	    my $contig_len = $genseq_len{$genome_tag}->{$feat_hash{$feat_name}->{'contig'}};
	    print STDERR "$feat_name $feat_hash{$feat_name}->{'anno'} $feat_hash{$feat_name}->{'gtag'} $genome_tag\n" if ($DEBUG);
	    if (!$suppress) {
		print OUTFILE ">$genome_tag";
		print OUTFILE " $feat_hash{$feat_name}->{'contig'} $feat_name $fivep $threep $feat_hash{$feat_name}->{'anno'}\n";
	    }
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
	    if ($seq_len <= 0) { #should not happen
		die ("ERROR: coordinates on contig sequence resulted in negative seq_len $seq_len for $feat_name!\n");
	    } else {
		if (($seq_len % 3) == 0) {
		    $div_by_three++;
		}
		if ($compute_all || ($target_id ne "")) {
		    if (!$gene_count) {
			$single_genome = $genome_tag;
		    }
		    $genome_seqs[$index] = $sequence;
		    if ($genome_tag eq $target_id) {
			$target_sequence = $sequence;
		    }
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
		    }
		}
		if (!$suppress) {
		    my $pos;
		    my $tmp_seq_len = $seq_len;
		    for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
			print OUTFILE substr($sequence, $pos, 60), "\n";
			$tmp_seq_len -= 60;
		    }
		    print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
		}
	    }
	    if ($remake_files || $compute_all || ($target_id ne "")) {
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
	    $gene_count++;
	    $index++;
	}
	if (!$suppress) {
	    close (OUTFILE);
	}
	my $mean;
	my $median;
	my $median_25;
	my $median_75;
	my $stddev;
	if ($remake_files || $compute_all || ($target_id ne "")) {
	    if ($gene_count > 0) {
		$mean = $sum / $gene_count;
		@sizes = sort {$a <=> $b} @sizes;
		$median = ($gene_count % 2) ? $sizes[($gene_count / 2)] : (($sizes[(($gene_count / 2) - 1)] + $sizes[($gene_count / 2)]) / 2);
		$median_25 = $sizes[($gene_count  / 4)];
		$median_75 = $sizes[(($gene_count * 3) / 4)];
		$stddev = sqrt(($sumsquared - ($mean * $mean * $gene_count)) / (($gene_count > 1) ? ($gene_count - 1) : 1));
		if ($remake_files) {
		    print SIZEFILE "$reduced_cluster_num\t$gene_count\t\t$min\t$max\t$median\t$mean\t$stddev\t\t$median_25\t$median_75\n"; # need to add st. dev. do not have "connectivity" or average %identity at this point
		    print OUTMATCHFILE "$reduced_cluster_num\t$out_line\n";
		}
	    }
	}
	if ($compute_all || ($target_id ne "")) {
	    if ($compute_all) {
		$index = 0;
		foreach my $genome_tag (@genome_array) {
		    if (($gene_count == 1) && ($single_genome eq $genome_tag)) {
			$uniq_clus{$genome_tag}++;
		    } elsif (defined $genome_seqs[$index]) {
			if ($feat_pres{$genome_seqs[$index]} == 1) {
			    if (((100 * $gene_count) / $genome_number) >= 75) {
				$uniq_clus_alle_75_100{$genome_tag}++;
			    } elsif (((100 * $gene_count) / $genome_number) <= 25) {
				$uniq_clus_alle_0_25{$genome_tag}++;
			    } else {
				$uniq_clus_alle_25_75{$genome_tag}++;
			    }
			}
			my $seq_len = length($genome_seqs[$index]);
			if ($seq_len < ($median_25 - (0.1 * $median))) {
			    $short_clus{$genome_tag}++;
			} elsif ($seq_len > ($median_75 + (0.1 * $median))) {
			    $long_clus{$genome_tag}++;
			}
			if ($div_by_three > (0.7 * $gene_count)) { # best approximation for frameshift
			    if (($seq_len % 3) != 0) {
				$frameshift{$genome_tag}++;
			    }
			}
		    } else {
			if (((100 * $gene_count) / $genome_number) >= 75) {
			    $missing_75c{$genome_tag}++;
			}
		    }
		    $index++;
		}
	    } else {
		if (($gene_count == 1) && ($single_genome eq $target_id)) {
		    $uniq_clus{$target_id}++;
		} elsif ($target_sequence ne "") {
		    if ($feat_pres{$target_sequence} == 1) {
			if (((100 * $gene_count) / $genome_number) >= 75) {
			    $uniq_clus_alle_75_100{$target_id}++;
			} elsif (((100 * $gene_count) / $genome_number) <= 25) {
			    $uniq_clus_alle_0_25{$target_id}++;
			} else {
			    $uniq_clus_alle_25_75{$target_id}++;
			}
		    }
		    my $seq_len = length($target_sequence);
		    if ($seq_len < ($median_25 - (0.1 * $median))) {
			$short_clus{$target_id}++;
		    } elsif ($seq_len > ($median_75 + (0.1 * $median))) {
			$long_clus{$target_id}++;
		    }
		    if ($div_by_three > (0.7 * $gene_count)) { # best approximation for frameshift
			if (($seq_len % 3) != 0) {
			    $frameshift{$target_id}++;
			}
		    }
		} else {
		    if (((100 * $gene_count) / $genome_number) >= 75) {
			$missing_75c{$target_id}++;
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
    unless (open (PGGFILE, "<$pgg_file") )  {
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
	my $out_line = join("\t", @edge_values);
	my $out_cluster1;
	my $out_cluster2;
	if ($edge_id =~ /\((\d+)_([35]),(\d+)_([35])\)/) {
	    $edge_id = "edge".$1."_".$2."to".$3."_".$4;
	    $cluster1 = $1;
	    $cluster2 = $3;
	    $whichend1 = $2;
	    $whichend2 = $4;
	    $out_cluster1 = $renumber[$cluster1];
	    $out_cluster2 = $renumber[$cluster2];
	} else {
	    die ("ERROR: Bad edge formatting $edge_id in file $pgg_file.\n");
	}
	    
	if (!$suppress) {
	    unless (open (OUTFILE, ">$multifastadir/$edge_id.fasta") )  {
		die ("ERROR: cannot open file $multifastadir/$edge_id.fasta\n");
	    }
	}
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
	foreach my $edge_value (@edge_values) {
	    $genome_tag = shift @tmp_array;
	    if ($edge_value == 0) { #this is a placeholder and can be skipped
		$index++;
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
	    if (!$suppress) {
		print OUTFILE ">$genome_tag";
		print OUTFILE " $contig1 $feat_name1 $start1 $end1 $feat_name2 $start2 $end2\n";
	    }
	    my $seq_len;
	    my $sequence;
	    my $contig_len = $genseq_len{$genome_tag}->{$contig1};
	    if ((abs($start2 - $end1) > (0.9 * $contig_len)) && $is_circular{$genome_tag}->{$contig1}) {
		#print "$genome_tag-$contig1($contig_len) $start1:$end1 - $start2:$end2\n";
		if ($start1 < $end2) {
		    my $tmp_seq = "";
		    my $beg_offset = 0;
		    my $end_offset = $contig_len - $start2;
		    if ($start2 > $contig_len) {
			$beg_offset = $start2 - $contig_len;
		    }
		    if ($end1 <= 0) {
			$end_offset += ($end1 - 1);
		    }
		    $seq_len = ($contig_len - $start2) + ($end1 - 1);
		    #print "if $seq_len:$start2:$end1:$beg_offset:$end_offset\n";
		    if ($seq_len <= 0) {
			$seq_len = 0;
			$sequence = "";
		    } else {
			if ($end_offset > 0) {
			    $tmp_seq = substr($genseq_hash{$genome_tag}->{$contig1}, $start2, $end_offset);
			}
			if ($end1 > 1) {
			    $tmp_seq .= substr($genseq_hash{$genome_tag}->{$contig1}, $beg_offset, (($end1 - 1) - $beg_offset));
			}
			$sequence = reverse($tmp_seq);
			$sequence =~ tr/AGCTYRWSKMDVHBagctyrwskmdvhb/TCGARYWSMKHBDVtcgarywsmkhbdv/;
		    }
		} else {
		    $sequence = "";
		    my $beg_offset = 0;
		    my $end_offset = $contig_len - $end1;
		    if ($end1 > $contig_len) {
			$beg_offset = $end1 - $contig_len;
		    }
		    if ($start2 <= 0) {
			$end_offset += ($start2 - 1);
		    }
		    $seq_len = ($contig_len - $end1) + ($start2 - 1);
		    #print "else $seq_len:$end1:$start2:$beg_offset:$end_offset\n";
		    if ($seq_len <= 0) {
			$seq_len = 0;
			$sequence = "";
		    } else {
			if ($end_offset > 0) {
			    $sequence = substr($genseq_hash{$genome_tag}->{$contig1}, $end1, $end_offset);
			}
			if ($start2 > 1) {
			    $sequence .= substr($genseq_hash{$genome_tag}->{$contig1}, $beg_offset, (($start2 - 1) - $beg_offset));
			}
		    }
		}
	    }elsif ((($start1 < $end2) && (($end1 + 1) >= $start2)) || (($start1 >= $end2) && (($end1 - 1) <= $start2))) {
		# cluster/gene features overlap or abut on contig so edge sequence is empty set length to 0 and sequence to empty string but output empty line for fasta
		$seq_len = 0;
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
	    if ($remake_files || $compute_all || ($target_id ne "")) {
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
		} elsif ($seq_len == 0) {
		    if (!$suppress) {
			print OUTFILE "\n";
		    }
		} else {
		    if (!$suppress) {
			my $pos;
			my $tmp_seq_len = $seq_len;
			for ( $pos = 0 ; $tmp_seq_len > 60 ; $pos += 60 ) {
			    print OUTFILE substr($sequence, $pos, 60), "\n";
			    $tmp_seq_len -= 60;
			}
			print OUTFILE substr($sequence, $pos, $tmp_seq_len), "\n";
		    }
		}
		if ($compute_all || ($target_id ne "")) {
		    if (!$gene_count) {
			$single_genome = $genome_tag;
		    }
		    $genome_seqs[$index] = $sequence;
		    if ($genome_tag eq $target_id) {
			$target_sequence = $sequence;
		    }
		    if (defined $feat_pres{$sequence}) {
			$feat_pres{$sequence}++;
		    } else {
			$feat_pres{$sequence} = 1;
		    }
		}
	    }
	    $gene_count++;
	    $index++;
	}
	if (!$suppress) {
	    close (OUTFILE);
	}
	my $mean;
	my $median;
	my $median_25;
	my $median_75;
	my $stddev;
	if ($remake_files || $compute_all || ($target_id ne "")) {
	    if ($gene_count > 0) {
		$mean = $sum / $gene_count;
		@sizes = sort {$a <=> $b} @sizes;
		$median = ($gene_count % 2) ? $sizes[($gene_count / 2)] : (($sizes[(($gene_count / 2) - 1)] + $sizes[($gene_count / 2)]) / 2);
		$median_25 = $sizes[($gene_count / 4)];
		$median_75 = $sizes[(($gene_count * 3) / 4)];
		$stddev = sqrt(($sumsquared - ($mean * $mean * $gene_count)) / (($gene_count > 1) ? ($gene_count - 1) : 1));
		if ($remake_files) {
		    print SIZEFILE "($out_cluster1", "_$whichend1,$out_cluster2", "_$whichend2)\t$gene_count\t\t$min\t$max\t$median\t$mean\t$stddev\t\t$median_25\t$median_75\n"; # need to add st. dev. do not have "connectivity" or average %identity at this point
		    print OUTPGGFILE "($out_cluster1", "_$whichend1,$out_cluster2", "_$whichend2)\t$out_line\n";
		}
	    }
	}
	if ($cluster1 <= $cluster2) { # only need to do this for one orientation of the edge - not sure if the clusters can be equal or if there are two edges in this case - do a 3' 5' test?
	    if ($compute_all || ($target_id ne "")) {
		if ($compute_all) {
		    $index = 0;
		    foreach my $genome_tag (@genome_array) {
			if (($gene_count == 1) && ($single_genome eq $genome_tag)) {
			    $uniq_edge{$genome_tag}++;
			} elsif (defined $genome_seqs[$index]) {
			    if ($feat_pres{$genome_seqs[$index]} == 1) {
				if (((100 * $gene_count) / $genome_number) >= 75) {
				    $uniq_edge_alle_75_100{$genome_tag}++;
				} elsif (((100 * $gene_count) / $genome_number) <= 25) {
				    $uniq_edge_alle_0_25{$genome_tag}++;
				} else {
				    $uniq_edge_alle_25_75{$genome_tag}++;
				}
			    }
			    my $seq_len = length($genome_seqs[$index]);
			    #print "$genome_tag($mean):$median_25:$median:$median_75:$seq_len\n";
			    if ($seq_len < ($median_25 - (0.1 * $median))) {
				$short_edge{$genome_tag}++;
				#print "SHORT\n";
			    } elsif ($seq_len > ($median_75 + (0.1 * $median))) {
				$long_edge{$genome_tag}++;
				#print "LONG\n";
			    }
			} else {
			    if (((100 * $gene_count) / $genome_number) >= 75) {
				$missing_75e{$genome_tag}++;
			    }
			}
			$index++;
		    }
		} else {
		    if (($gene_count == 1) && ($single_genome eq $target_id)) {
			$uniq_edge{$target_id}++;
		    } elsif ($target_sequence ne "") {
			if ($feat_pres{$target_sequence} == 1) {
			    if (((100 * $gene_count) / $genome_number) >= 75) {
				$uniq_edge_alle_75_100{$target_id}++;
			    } elsif (((100 * $gene_count) / $genome_number) <= 25) {
				$uniq_edge_alle_0_25{$target_id}++;
			    } else {
				$uniq_edge_alle_25_75{$target_id}++;
			    }
			}
			my $seq_len = length($target_sequence);
			if ($seq_len < ($median_25 - (0.1 * $median))) {
			    $short_edge{$target_id}++;
			} elsif ($seq_len > ($median_75 + (0.1 * $median))) {
			    $long_edge{$target_id}++;
			}
		    } else {
			if (((100 * $gene_count) / $genome_number) >= 75) {
			    $missing_75e{$target_id}++;
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

sub process_single_cores                                               # Read in the set of clusters which are single copy core, renumber and output
{
    unless (open(CLUSTER_CORES, "<", $single_cores)) {
	die ("cannot open cores file: $single_cores!\n");
    }
    unless (open(OUT_CORES, ">", "$basedir/single_copy_clusters.txt")) {
	die ("cannot open cores file: $basedir/single_copy_clusters.txt!\n");
    }
    while (my $line = <CLUSTER_CORES>) {
	chomp $line;
	$line =~ s/\s*//g; # remove all whitespace characters
	$line =~ s/^.*_//; # remove centroid_, medoid_, cluster_ or any other verbiage before the cluster number
	if ($renumber[$line]) {
	    print OUT_CORES "$renumber[$line]\n";
	}
    }
    close(CLUSTER_CORES);
    close(OUT_CORES);
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
     -b: directory path for where to put output multifasta files[DEFAULT = PWD]
     -B: directory path for where to put other output files[DEFAULT = PWD]
     -g: two colums tab delimited file (input): col1 genome_tag used in panoct in same order as panoct; col2 file name for the fasta contig file for this genome
     -m: panoct matchtable file (input)
     -a: combined panoct attribute file (input)
     -p: a panoct pan_genome_graph_file
     -t: target genome id if stats only wanted for one genome
     -A: output stats for all genomes
     -S: suppress multifasta output
     -M: medoids input file name - output will be in -B directory (medoids.fasta)
     -c: output size files in -B directory computed cluster sizes(program  adds cluster_sizes.txt) and edge sizes(program  adds edge_sizes.txt) files
     -P: output files in -B directory new matchtable(matchtable.txt) and pgg(pgg.txt) files
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
print "Getting genome names and contig sequences from $genomes_file_name\n";
&get_genomes;
print "Gathering gene coordinates and attribute information from $att_file\n";
foreach my $genome_tag (@genome_array) {
    $uniq_clus{$genome_tag} = 0;            # key = genome ID, value = number of unique/singleton clusters
    $uniq_edge{$genome_tag} = 0;            # key = genome ID, value = number of unique/singleton edges
    $short_clus{$genome_tag} = 0;           # key = genome ID, value = number of short clusters
    $short_edge{$genome_tag} = 0;           # key = genome ID, value = number of short edges
    $long_clus{$genome_tag} = 0;            # key = genome ID, value = number of long clusters
    $long_edge{$genome_tag} = 0;            # key = genome ID, value = number of long edges
    $frameshift{$genome_tag} = 0;           # key = genome ID, value = number of probable frameshifts
    $missing_75c{$genome_tag} = 0;          # key = genome ID, value = number of clusters missing  for clusters in 75-100% of genomes
    $missing_75e{$genome_tag} = 0;          # key = genome ID, value = number of edges missing  for edges in 75-100% of genomes
    $uniq_clus_alle_75_100{$genome_tag} = 0;# key = genome ID, value = number of unique alleles for clusters in 75-100% of genomes
    $uniq_clus_alle_25_75{$genome_tag} = 0; # key = genome ID, value = number of unique alleles for clusters in 25-75% of genomes
    $uniq_clus_alle_0_25{$genome_tag} = 0;  # key = genome ID, value = number of unique alleles for clusters in 0-25% of genomes
    $uniq_edge_alle_75_100{$genome_tag} = 0;# key = genome ID, value = number of unique alleles for edges in 75-100% of genomes
    $uniq_edge_alle_25_75{$genome_tag} = 0; # key = genome ID, value = number of unique alleles for edges in 25-75% of genomes
    $uniq_edge_alle_0_25{$genome_tag} = 0;  # key = genome ID, value = number of unique alleles for edges in 0-25% of genomes
}    
&get_attributes;
print "Reading matchtable from $matchtable_file and outputting cluster multifasta files to $basedir\n";
&process_matchtable;
print "Reading pan-genome graph from $pgg_file and outputting edge multifasta files to $basedir\n";
&process_pgg;
if ($compute_all || ($target_id ne "")) {
    unless (open (STATFILE, ">$basedir/cluster_stats.txt") ) {
	die ("ERROR: cannot open file $basedir/cluster_stats.txt\n");
    }
    print STATFILE "Genome\tuniq_clus\tuniq_clus_alle_0_25\tuniq_clus_alle_25_75\tuniq_clus_alle_75_100\tuniq_edge\tuniq_edge_alle_0_25\tuniq_edge_alle_25_75\tuniq_edge_alle_75_100\tshort_clus\tlong_clus\tshort_edge\tlong_edge\tframshift\tmissing_75clus\tmissing_75edge\n";
    if ($compute_all) {
	foreach my $genome_tag (@genome_array) {
	    print STATFILE "$genome_tag\t$uniq_clus{$genome_tag}\t$uniq_clus_alle_0_25{$genome_tag}\t$uniq_clus_alle_25_75{$genome_tag}\t$uniq_clus_alle_75_100{$genome_tag}\t$uniq_edge{$genome_tag}\t$uniq_edge_alle_0_25{$genome_tag}\t$uniq_edge_alle_25_75{$genome_tag}\t$uniq_edge_alle_75_100{$genome_tag}\t$short_clus{$genome_tag}\t$long_clus{$genome_tag}\t$short_edge{$genome_tag}\t$long_edge{$genome_tag}\t$frameshift{$genome_tag}\t$missing_75c{$genome_tag}\t$missing_75e{$genome_tag}\n";
	}
    } else {
	print STATFILE "$target_id\t$uniq_clus{$target_id}\t$uniq_clus_alle_0_25{$target_id}\t$uniq_clus_alle_25_75{$target_id}\t$uniq_clus_alle_75_100{$target_id}\t$uniq_edge{$target_id}\t$uniq_edge_alle_0_25{$target_id}\t$uniq_edge_alle_25_75{$target_id}\t$uniq_edge_alle_75_100{$target_id}\t$short_clus{$target_id}\t$long_clus{$target_id}\t$short_edge{$target_id}\t$long_edge{$target_id}\t$frameshift{$target_id}\t$missing_75c{$target_id}\t$missing_75e{$target_id}\n";
    }
    close (STATFILE);
}
if ($remake_files) {
    &process_medoids;
    &process_single_cores;
}
exit(0);
