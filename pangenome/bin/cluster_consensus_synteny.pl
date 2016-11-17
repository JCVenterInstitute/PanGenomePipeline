#!/usr/bin/env perl
# $Id: cluster_consensus.pl 35964 2011-12-07 14:30:16Z ekelsey $

##---------------------------------------------------------------------------##
##  File: cluster_consensus.pl
##
##  Author:
##        Malay <malay@bioinformatics.org>
##
##  Description:
##
#******************************************************************************
#* Copyright (C) 2011 Malay K Basu <malay@bioinformatics.org>
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

cluster_consensus.pl - Create consensus cluster from individual clustering method.

=head1 SYNOPSIS

cluster_consensus.pl [options] <method1.result> <method2.result> >result.txt


=head1 DESCRIPTION

Very initial implemenation of creating consensus among all methods in SaucePan.


=head1 ARGUMENTS 

Individual method results files.




=head1 OPTIONS

=item B<--ortholog|-o>

Create indexed orthologs. Guranteed to give you maximum one gene per genome.

=item B<--blast|-b blast_file>

The blast result file for all against all in -m 9 format.

=item B<--genomes| -g <directory> >

A directory containing individual genome multifasta files.

=item B<--coords|c <coordinates_file> >

A file in the same format as gene_att file of panoct containing coordinates.

=item <B--help|-h>

Print this help.


=back


=head1 SEE ALSO

=head1 COPYRIGHT

Copyright (c) 2011 Malay K Basu <malay@bioinformatics.org>

=head1 AUTHORS

Malay K Basu <malay@bioinformatics.org>

=cut

##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../lib";
use SeqToolBox;
use SeqToolBox::SeqDB;
use SeqToolBox::Tools::Set;
use Pangenome::Clusterset;
use Pangenome::Cluster;
use Pangenome::BLASTResult;
use Carp;

##---------------------------------------------------------------------------##
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
##---------------------------------------------------------------------------##

my %options = ();    # this hash will have the options
my $blast_result = Pangenome::BLASTResult->new();
my $clusters;        # $clusters->{$tag_name} = $clusterset_obj;

#my $data;            # $data->{tag_name}->{gi}->{cluster} = 1
my @genomes;         #array of genome name
my $genome_data;     # hash containing key=gi value=genome_name
my @methods;         #array of all the method names
my %all_genes;       #hash containing all genes included in all clusters
my $result;
my %centroids;
my %centroid_affinity;
my %marged_centroids;
my $bigram_data;
my $genome_map_index;
my $genome_map_array;
#my $coordinates_data;

#
# Get the supplied command line options, and set flags
#

GetOptions( \%options, 'help|h', 'blast|b=s', 'ortholog|o', 'genomes|g=s','coords|c=s' )
  || pod2usage( -verbose => 0 );

_check_params( \%options );

my @files = @ARGV;
read_blast_result( $options{blast} );    #sets $blast_result

foreach my $file (@files) {
	warn "Read data\n";
	read_data($file);                    #sets clusters and data
}

if ( $options{ortholog} ) {
	my ( $g, $gd ) = read_genome_files( $options{genomes} );
	@genomes     = @{$g};
	$genome_data = $gd;
}

@methods = keys %{$clusters};

## resolve conflict
foreach my $m (@methods) {
	warn "Resolve conflict\n";
	while (1) {
		my $retval = $clusters->{$m}->resolve_conflict($blast_result);
		last if $retval;
	}
}

find_all_genes();    # sets all_genes

my $total_gene_num = scalar( keys %all_genes );
print STDERR "Total $total_gene_num genes found in clusters\n";
#
##print STDERR join (" ",keys %all_genes),"\n";
#
##resolve conflict in cluster data
#foreach my $g ( keys %all_genes ) {
#
#	foreach my $m (@methods) {
#
#		#get all the clusters for the method where a gene belong
#		my @clusters = get_clusters( $m, $g );
#		next unless @clusters;
#
#	 # what is the true cluster if the gene belongs to more than one clusters?
#		my $true_cluster;
#
#		if ( @clusters > 1 )
#		{ #We need to resolve conflict only when a gene belongs to more than one cluster
#			print STDERR "$g belongs to @clusters in $m. Resolving...\n";
#			$true_cluster = resolve_conflict( $m, $g );
#		} else {
#			$true_cluster = $clusters[0];
#		}
#	}
#}

#find_centroids();
#warn Dumper(%centroids);
#centroid_affinity();
#
my $intersections;
my $count = 0;
my $intersection_members;
my %seen;

#print STDERR Dumper($clusters);
foreach my $g ( keys %all_genes ) {

	#	print STDERR "Second interation\n";
	next if exists $seen{$g};
	my $member_cluster;

	for ( my $i = 0 ; $i < @methods ; $i++ ) {

		#		my @cluster = get_clusters( $methods[$i], $g );
		my @cluster = $clusters->{ $methods[$i] }->get_cluster_by_gene($g);

		#		print STDERR "@cluster $methods[$i] $g\n";
		#		print STDERR "Cluster: $cluster[0]\n";

		if ( @cluster && scalar(@cluster) > 1 ) {
			die "More than one clusters found for $g for $methods[$i]";
		}
		$member_cluster->{ $methods[$i] } = $cluster[0] if @cluster;
	}
	my @intersection = find_common_gene_set($member_cluster);

	next unless @intersection;

	print STDERR "Intersection: @intersection\n";
	$count++;
	$intersections->{$count} = \@intersection;
	foreach my $i (@intersection) {
		$seen{$i} = 1;
		$intersection_members->{$i}->{$count} = 1;
	}

}

#
#Merge sequence to merge duplicate intersections


my $count1 = 0;
my %seen1;

foreach my $i ( keys %{$intersection_members} ) {
	next if exists $seen1{$i};
	my @clusters = keys %{ $intersection_members->{$i} };

# @clusters will be more than one if the gene belongs to more than one intersection
	if ( @clusters > 1 ) {
		print STDERR "$i found in more than one intersections. Merging\n";
		my @merged_members = merge_intersection(@clusters);

		#print STDERR "@merged_members\n";
		foreach my $m (@merged_members) {
			$seen1{$m} = 1;
		}
	}
}

# We need to remap intersection members for the final time
$intersection_members = undef;

foreach my $i (keys %{$intersections}) {
	my @members = @{$intersections->{$i}};
	foreach my $j (@members) {
		$intersection_members->{$j} = $i;
	}
}

my @genes_in_intersection = keys %{$intersection_members};
print STDERR "Total genes in intersection:", scalar(@genes_in_intersection), "\n";

die "Could not find intersection for all the genes" unless (scalar(@genes_in_intersection) == $total_gene_num); 

#Time to read synteny information
if ($options{coords}) {
	read_coords_file ();
}


my $all_orthologous_genes;
my $orthologous_cluster;
my $ortholog_counter = 0;

#if ( $options{ortholog} ) {

##	print STDERR Dumper($intersections);

foreach my $c ( sort { $a <=> $b } keys %{$intersections} ) {

	#
	my @members = @{ $intersections->{$c} };
	my $g_m     = get_genome_classified(@members);
	my @number_of_genomes;
	my @orthologs;

	foreach my $gen1 ( keys %{$g_m} ) {
		push @number_of_genomes, $gen1;

		#			foreach my $gen2 ( keys %{$g_m} ) {
		my @array1 = keys %{ $g_m->{$gen1} };

		if ( my $best_scoring_id
			 = get_best_scoring_hits_from_genome( $gen1, \@array1, \@members ) )
		{

			push @orthologs, $best_scoring_id;
		}

		#				my @array2 = keys %{$g_m->{$gen2}};
		#				my ($a1, $b1) = get_best_hit (\@array1, \@array2);
		#				print STDERR "==============$c==================\n";
		#				print STDERR "Genome1: $gen1, Genome2:$gen2 A1:$a1,$b1\n";
		#				my ($a2, $b2) = get_best_hit (\@array2, \@array1);
		#				print STDERR "Genome1: $gen1, Genome2:$gen2 A1:$a2,$b2\n";
		#				if ($a1 eq $b2 && $a2 eq $b1) {
		#					$orthologs{$a1} = 1;
		#					$orthologs{$a2} = 1;
		#					$orthologs{$b2} = 1;
		#					$orthologs{$b1} = 1;
		#				}
		#			}
		#
		#		}
		#
		#	}

	}

	#		foreach my $o {@orthologs} {
	#			$all_orthologus_genes{$o} = 1;
	#		}
	if ( scalar(@orthologs) > scalar(@number_of_genomes) ) {
		print STDERR "Error: @orthologs\n";
		die "Could not dtermine orthologs\n";
	}

	#		foreach my $g (@genes) {
	#			if ( exists $data->{$tag_name}->{$g} ) {
	#				print STDERR "$g found in more than one cluster in $tag_name\n";
	#				$data->{$tag_name}->{$g}->{$cluster} = 1;
	#			} else {
	#				$data->{$tag_name}->{$g}->{$cluster} = 1;
	#			}
	#		}
	my $already_seen = 0;
	next unless @orthologs;

	foreach my $o_g (@orthologs) {

		if ( exists $all_orthologous_genes->{$o_g} ) {

			#				print "ortholog $o_g exists\n";
			$all_orthologous_genes->{$o_g} = 1;
			$already_seen = 1;
		} else {
			$all_orthologous_genes->{$o_g} = 1;
		}
	}
	if ($options{ortholog}) {
		if (scalar(@orthologs) != @genomes) {
			next;
		}
	}
	unless ($already_seen) {
		print ++$ortholog_counter, "\t", join( "\t", @orthologs ), "\n";
	}

	#			$all_orthologous_genes->{$o_g} = 1;

	#			}else {
	#				$all
}

#		my @orthos = keys %orthologs;

#}

#merging step;

#foreach my $o_g (keys %all_orhologs_genes)

exit(0);

######################## S U B R O U T I N E S ############################

# for each method and clusters calculate the gene with the best score with the other
# other genes. This is considered to be the centroid of the cluster. Directly stores
# the result into centroids hash.

sub read_coords_file {
	my $file = $options{coords};
	die "Could not read coordinate information" unless (-s $file);
	open (FILE, $file) || die "Could not read $file\n";
	my $coordinates_data;
	while (my $line = <FILE>) {
		chomp $line;
		my @f = split (/\t/, $line);
		my $gene = $f[1];
		my $start = $f[2];
		my $end = $f[3];
		my $genome = $f[5];
		if ($start > $end) {
			my $tmp = $end;
			$end = $start;
			$start = $tmp;
		}
		# Sorting according to the start position. May not be ideal
		$coordinates_data->{$genome}->{$gene} = $start;
		
	}
	close (FILE);
	foreach my $genome (keys %{$coordinates_data}) {
		my @sorted_genes = sort {$coordinates_data->{$genome}->{$a} <=>$coordinates_data->{$genome}->{$b} } keys %{$coordinates_data->{$genome}};
		print STDERR join("\n", @sorted_genes),"\n";
		for (my $i = 0; $i <@sorted_genes; $i++) {
			$genome_map_index->{$genome}->{$sorted_genes[$i]} = $i;
			my $previous_gene;
			if ($i == 0) {
				$previous_gene = $sorted_genes[$#sorted_genes];
			}
			my $previous_class;
			my $this_class;
			if (exists $intersection_members->{$sorted_genes[$i]}) {
				$this_class = $intersection_members->{$sorted_genes[$i]};
			}else {
				$this_class = 'unknown';
			}
			if (exists $intersection_members->{$previous_gene}) {
				$previous_class = $intersection_members->{$previous_gene};
			}else {
				$previous_class = 'unknown';
			}
			$bigram_data->{$genome}->{$previous_class}->{$this_class}++;
			
			
		}
		$genome_map_array->{$genome} = \@sorted_genes;
	}
	
}

sub find_centroids {
	print STDERR ("Finding centroids...");
	my $counter = 0;
	foreach my $tag ( keys %{$clusters} ) {
		my $clusterset = $clusters->{$tag};
		my @clusters   = $clusterset->get_all_clusters();
		foreach my $c (@clusters) {

		}

		#		foreach my $cluster ( keys %{ $clusters->{$tag} } ) {

	 #			my $best_score = 0;
	 #			my $c_gene;
	 #			my $geneset
	 #			  = get_cluster_members_by_cluster_name( $tag, $cluster );
	 #
	 #			foreach my $gene ( keys %{$geneset} ) {
	 #				my $score = get_avg_par_score( $geneset, $gene );
	 #				if ( $score > $best_score ) {
	 #					$best_score = $score;
	 #					$c_gene     = $gene;
	 #				} elsif ( $score == $best_score ) {
	 #
	 #		  # Small probability but actually it's possible that
	 #		  # more than one gene can be centroid. We are in trouble. I need to
	 #		  # handle this situation but for the time being just kill the program
	 #		  # TODO: Handle more than on centroid edge cases.
	 #
	 #					my $existing_gene = $c_gene || " ";
	 #					print STDERR
	 #					  "Centroid: $gene has the same score as $existing_gene
	 #							in $tag cluster: $cluster, trying to resolve\n";
	 #
	 #		#					my $best_gene = find_best_centroid($tag, $existing_gene, $gene);
	 #		#					croak(
	 #		#						"$gene has the same score as $existing_gene
	 #		#							in $tag cluster: $cluster\n" );
	 #				}
	 #			}
	 #
	 #			$centroids{$tag}->{$cluster} = $c_gene;
	}

	#}
	warn "done\n";
}

# given a method and two genes the function tries first to resolve the tie by finding how
# how many methods the gene belongs to. Even after this it can't resolve the conflict then it tries resolve the conflict by the total genes in a cluster;

sub find_best_centroid {
	my ( $method, @genes ) = @_;

}

sub merge_intersection {
	my @clusters = @_;

	#print STDERR "@clusters\n";
	my %all_members;
	foreach my $c (@clusters) {
		my @members = @{ $intersections->{$c} };
		foreach my $m (@members) {
			$all_members{$m} = 1;
		}
	}
	my @merged_members = keys %all_members;
	for ( my $i = 0 ; $i < @clusters ; $i++ ) {
		if ( $i == 0 ) {
			$intersections->{ $clusters[$i] } = \@merged_members;
		} else {
			delete $intersections->{ $clusters[$i] };
		}
	}
	return @merged_members;
}

sub get_best_scoring_hits_from_genome {

	#get two array refs
	my ( $genome, $genome_members, $intersection_members ) = @_;
#	my %inter_members = map { $_ => 1 } @{$intersection_members};
	my $best_score = 0;
	my $best_scoring_id;

	foreach my $g ( @{$genome_members} ) {
#		if ( my $score = get_avg_par_score( \%inter_members, $g ) ) {
		my $blast_score = $blast_result->get_avg_distance_score($g, $intersection_members);
		my $bigram_score = get_bigram_score ($g, $genome);
#		if (my $score = $blast_result->get_avg_distance_score($g, $intersection_members))
		my $score = $blast_score + $bigram_score;
			if ( $score > $best_score ) {
				$best_score      = $score;
				$best_scoring_id = $g;
			}
#		}
	}
	return $best_scoring_id ? $best_scoring_id : undef;
}

sub get_bigram_score {
	my ($gene, $genome ) = @_;
	my $this_class = get_class($gene);
	my $previous_class = get_previous_gene_class($gene, $genome);
	my $total_bigram = 0;
	my $denom  = 0;
	foreach my $g (@genomes) {
		next if ($g eq $genome);
		if (exists $bigram_data->{$g}->{$previous_class}->{$this_class}) {
			$total_bigram += $bigram_data->{$g}->{$previous_class}->{$this_class};
		}
		foreach my $c (keys %{$bigram_data->{$g}->{$previous_class}}) {
			$denom += $bigram_data->{$g}->{$previous_class}->{$c};
		}
	}
	#Laplace correction
	my @total_class = keys (%{$intersections});
	my $num = scalar (@total_class);
	return ($total_bigram + 1)/($denom + $num);
	
}

sub get_previous_gene_class {
	my ($gene, $genome) =@_;
	my $this_gene_index = $genome_map_index->{$genome}->{$gene};
	my $previous_gene;
	if ($this_gene_index == 0) {
		my @gene_array = @{$genome_map_array->{$genome}};
		$previous_gene = $gene_array[$#gene_array]; 
	}else {
		my $prv_gene_index = $this_gene_index - 1;
		$previous_gene = $genome_map_array->{$genome}->[$prv_gene_index];
	}
	return get_class($previous_gene);
}

sub get_class {
	my $gene = shift;
	if (exists $intersection_members->{$gene}) {
		return $intersection_members->{$gene};
	}else {
		return 'unknnown';
	}
}

sub get_best_hit {
	my ( $array1, $array2 ) = @_;
	my $a;
	my $b;
	my $last_score = 0;

	for ( my $i = 0 ; $i < @{$array1} ; $i++ ) {
		for ( my $j = 0 ; $j < @{$array2} ; $j++ ) {

			my $score = get_score( $array1->[$i], $array2->[$j] );
			if ( $score > $last_score ) {
				$last_score = $score;
				$a          = $array1->[$i];
				$b          = $array2->[$j];
			}
		}
	}
	return $a, $b;

}

sub get_score {
	my ( $g1, $g2 ) = @_;
	my $score1 = 0;
	my $score2 = 0;
	$score1 = $blast_result->{$g1}->{$g2}
	  if exists $blast_result->{$g1}->{$g2};
	$score2 = $blast_result->{$g2}->{$g1}
	  if exists $blast_result->{$g2}->{$g1};

	if ( $score1 && $score2 ) {
		return ( $score1 + $score2 ) / 2;

	} elsif ( $score1 && !$score2 ) {
		return $score1;
	} elsif ( !$score1 && $score2 ) {
		return $score2;
	} else {
		print STDERR
		  "Error: Could not find blast score between $g1 and $g2\n";
		return;
	}

}

sub get_genome_classified {
	my @m = @_;

	#	print STDERR @m, "\n";
	my $return;

	foreach my $i (@m) {

		#		print STDERR "$i", "\n";
		if ( exists $genome_data->{$i} ) {

			#			print STDERR $genome_data->{$i}, "\n";
			$return->{ $genome_data->{$i} }->{$i} = 1;
		} else {
			die "Could not find genome data for ID $i\n";
		}
	}
	return $return;
}

sub find_common_gene_set {
	my $member_cluster = shift;

	#	print STDERR "Find common gene set\n";
	my $set = SeqToolBox::Tools::Set->new();

	foreach my $m ( keys %{$member_cluster} ) {

	  #		my $members
	  #		  = get_cluster_members_by_cluster_name( $m, $member_cluster->{$m} );

		#		my @array = keys %{$members};
		my @array = $member_cluster->{$m}->get_members_as_array();

		#		print STDERR "@array\n";
		$set->add_set( \@array );
	}

	return $set->get_all_intersection();
}

sub _check_params {
	my $opts = shift;
	pod2usage( -verbose => 2 ) if ( $opts->{help} || $opts->{'?'} );
	pod2usage( -verbose => 1 ) unless ( $opts->{'blast'} );

	if ( $opts->{orthologs} ) {
		die "Genome files must be given for otholog dertermination\n"
		  unless $opts->{genomes};
	}

}

# Takes a method name and a gi name to then returns only a single cluster name

#sub resolve_conflict {
#	my ( $method, $gi ) = @_;
#	my @clusters = get_clusters( $method, $gi );
#	my $highest_score_cluster = 0;
#	my $cluster;
#
#	foreach my $c (@clusters) {
#		my $members = get_cluster_members_by_cluster_name( $method, $c );
#
#		#		print STDERR Dumper($members);
#		my $score = get_avg_par_score( $members, $gi );
#		print STDERR "Cluster $c score: $score\n";
#
#		if ( $score > $highest_score_cluster ) {
#			$cluster               = $c;
#			$highest_score_cluster = $score;
#		}
#	}
#	print STDERR "$gi belongs to $cluster\n";
#
#	foreach my $c (@clusters) {
#		next if $c eq $cluster;
#		print STDERR "Removing $gi from $c\n";
#		my @array;
#
#		foreach my $member ( @{ $clusters->{$method}->{$c} } ) {
#			next if ( $member eq $gi );
#			push @array, $member;
#		}
#		$clusters->{$method}->{$c} = \@array;
#		delete( $data->{$method}->{$gi}->{$c} );
#
#	}
#	return $cluster;
#
#}

sub get_cluster_members_by_cluster_name {
	my ( $method, $c ) = @_;

	#	print STDERR "$method, $c\n";
	#	my @super;
	my %r;

	#	print STDERR Dumper($clusters);
	foreach my $c ( @{ $clusters->{$method}->{$c} } ) {

		#		print STDERR $c, "\n";
		$r{$c} = 1;
	}
	return \%r;

}

sub get_avg_par_score {
	my ( $diff_1, $gi ) = @_;
	my $count       = 0;
	my $total_score = 0;

	foreach my $i ( keys %{$diff_1} ) {
		if ( $i eq $gi ) { next }

		#		if ( exists $blast_result->{$gi}->{$i} ) {
		if ( my $score = get_score( $gi, $i ) ) {

			#			$total_score += $blast_result->{$gi}->{$i};
			$total_score += $score;
			$count++;
		}
	}
	if ( $count == 0 ) { return 0; }
	return $total_score / $count;
}

sub find_distance_in_a_cluster {

}

sub find_all_genes {

	foreach my $m (@methods) {

		#		foreach my $c ( keys %{ $clusters->{$m} } ) {
		#		print STDERR $m, "\n";
		my $c     = $clusters->{$m};
		my $genes = $c->get_all_genes();
		foreach my $g ( @{$genes} ) {
			$all_genes{$g} = 1;
		}

		#			foreach my $g ( @{ $clusters->{$m}->{$c} } ) {
		#				$all_genes{$g} = 1;
		#			}
		#		}
	}
}

# Given a method name and gene id return list of clusters
#sub get_clusters {
#	my ( $tag, $gi ) = @_;
#	return unless exists $data->{$tag}->{$gi};
#	my @clusters = keys( %{ $data->{$tag}->{$gi} } );
#
#	#	print STDERR "Inside get cluster @clusters\n";
#	return @clusters;
#}

#sub get_cluster_members {
#	my ( $tag, $gi ) = @_;
#	my @clusters = keys( %{ $data->{$tag}->{$gi} } );
#	my @super;
#	my %r;
#
#	foreach my $c (@clusters) {
#		push @super, @{ $clusters->{$tag}->{$c} };
#	}
#
#	foreach my $s (@super) {
#		$r{$s} = 1;
#	}
#	return \%r;
#}

sub read_blast_result {
	my $file = shift;
	open( my $infile, $file ) || die "Can't open $file\n";

	while ( my $line = <$infile> ) {
		next if $line =~ /^\#/;
		chomp $line;
		my @f = split( /\t/, $line );
		my $q;
		my $s;

		#		if ( $f[0] =~ /(\S+)\|/ ) {
		#			$q = $1;
		#		} else {
		$q = $f[0];

		#		}

		#		if ( $f[1] =~ /(\S+)\|/ ) {
		#			$s = $1;
		#		} else {
		$s = $f[1];

		#		}
		my $score = $f[11];
		$blast_result->add_hit( $q, $s, $score );

		#		$blast_result->{$q}->{$s} = $score;
	}
}

sub read_data {
	my $file = shift;
	my $tag  = $file;
	print STDERR $tag, "\n";
	$tag =~ s/\S*\///;
	my $tag_name;

	if ( $tag =~ /parsed\_(\S+)\.txt/ ) {
		$tag_name = $1;
	} else {
		$tag_name = $tag;
	}
	die "Could not parse the tag from filename" unless $tag_name;
	open( my $infile, $file ) || die "Can't open $file\n";
	my $clusterset = Pangenome::Clusterset->new($tag_name);

	while ( my $line = <$infile> ) {
		chomp $line;
		my ( $cluster_id, @genes ) = split( /\t/, $line );
		next unless ( @genes >= 2 );
		my $c = Pangenome::Cluster->new( $cluster_id, @genes );
		$clusterset->add_cluster($c);

		#		$clusters->{$tag_name}->{$cluster} = \@genes;
		#
		#		foreach my $g (@genes) {
		#			if ( exists $data->{$tag_name}->{$g} ) {
		#				print STDERR
		#				  "$g found in more than one cluster in $tag_name\n";
		#				$data->{$tag_name}->{$g}->{$cluster} = 1;
		#			} else {
		#				$data->{$tag_name}->{$g}->{$cluster} = 1;
		#			}
		#		}
	}

	#	print STDERR Dumper($clusters);
	close($infile);

	$clusters->{$tag_name} = $clusterset;
}

sub read_genome_files {
	my $dir = shift;
	opendir( DIR, $dir ) || die "Can't open $dir\n";
	my @files;

	while ( my $file = readdir(DIR) ) {
		if ( $file =~ /\.pep$/ || $file =~ /\.fas$/ || $file =~ /\.faa$/ ) {
			push @files, File::Spec->catfile( $dir, $file );
		}
	}
	die "Could not find multifasta files for genomes" unless (@files);
	return _parse_genome_data(@files);
}

sub _parse_genome_data {
	my @files = @_;
	my @array;
	my %data;

	foreach my $f (@files) {
		my @f = File::Spec->splitpath($f);
		my $genome_db;

		if ( $f[2] =~ /^(\w+)\./ ) {
			print STDERR $1, "\n";
			$genome_db = $1;
		} else {
			$genome_db = $f[2];
		}
		my $seqdb = SeqToolBox::SeqDB->new( -file => $f );

		while ( my $seq = $seqdb->next_seq() ) {
			my $id = $seq->get_id();

			#			$data{$id}->{$genome_db} = 1;
			$data{$id} = $genome_db;
		}
		push @array, $genome_db;
	}

	#	print STDERR Dumper (%data);
	return \@array, \%data;    # genome_list, id->genome map
}
