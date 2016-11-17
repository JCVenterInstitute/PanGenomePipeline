#!/usr/bin/env perl
# $Id: calculate_distance_matrix.pl 32957 2010-08-05 16:22:09Z mbasu $

##---------------------------------------------------------------------------##
##  File: calculate_distance_matrix.pl
##       
##  Author:
##        Malay <malay@bioinformatics.org>
##
##  Description:
##     
#******************************************************************************
#* Copyright (C) 2010 Malay K Basu <malay@bioinformatics.org> 
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

calculate_distance_matrix.pl - One line description.

=head1 SYNOPSIS

calculate_distance_matrix.pl [options] -o <option>


=head1 DESCRIPTION

Write a description of your prgram. 


=head1 ARGUMENTS 

=over 4

=item B<--option|-o>

First option.



=back

=head1 OPTIONS

Something here.


=head1 SEE ALSO

=head1 COPYRIGHT

Copyright (c) 2010 Malay K Basu <malay@bioinformatics.org>

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

my ($ref_genome, @files) = @ARGV;

my $data;
my $clusters;
my $results;

foreach my $file (@files) {
	read_data($file);
}
my @tags = keys (%{$data});
open (my $fh, $ref_genome) || die "Can't open $ref_genome\n";
while (my $line = <$fh>) {
	next unless $line =~ /^\>(\S+)\|/;
	my $gi = $1;
	die "Can't parse gi" unless $gi;
	calculate_similarity ($gi);
}
close ($fh);

print "Method1\tMethod2\tScore\tMatched_cluster\tSimilarity\n";

foreach my $i (@tags) {
	foreach my $j (@tags){
		print $i, "\t", $j, "\t", sprintf("%.2f",$results->{$i}->{$j}->{result}), "\t", $results->{$i}->{$j}->{number}, "\t",
		sprintf ("%.3f", $results->{$i}->{$j}->{result}/$results->{$i}->{$j}->{number} ), "\n";
	}
}


sub calculate_similarity {
	my $gi = shift;
	
	for (my $i = 0; $i < scalar(@tags); $i++) {
		for (my $j = 0; $j <scalar (@tags); $j++){
		my ($num_1, $num_2);
		
		if (!exists $data->{$tags[$i]}->{$gi} || !exists $data->{$tags[$j]}->{$gi}) {
			next;
		}
		my $cluster_1 = get_clusters($tags[$i], $gi);
		my $cluster_2 = get_clusters ($tags[$j], $gi);
		my $sim = get_similar_count ($cluster_1, $cluster_2);
		my $count_1 = scalar(keys %{$cluster_1});
		my $count_2 = scalar(keys %{$cluster_2});
		my $result = (($sim/$count_1) + ($sim/$count_2))/2;
		$results->{$tags[$i]}->{$tags[$j]}->{result} += $result;
		$results->{$tags[$i]}->{$tags[$j]}->{number}++;
				
	}
	}
}

sub get_similar_count {
	my ($array1, $array2) = @_;
	my $match_count = 0;
	foreach my $i (keys %{$array1}) {
		if (exists $array2->{$i}) {
			$match_count++;
		}
	}
	return $match_count;
	
}

sub get_clusters {
	my ($tag, $gi) = @_;
	my @clusters = keys (%{$data->{$tag}->{$gi}});
	my @super;
	my %r;
	foreach my $c(@clusters) {
		push @super, @{$clusters->{$tag}->{$c}};
	}
	foreach my $s (@super) {
		$r{$s} = 1;
	}
	return \%r;
}
sub read_data {
	my $file = shift;
	my $tag = $file;
	print STDERR $tag, "\n";
	$tag =~ s/\S*\///;
	$tag =~ /parsed\_(\S+)\.txt/;
	my $tag_name = $1;
	die "Could not parse the tag from filename" unless $tag_name;
	open (my $infile, $file) || die "Can't open $file\n";
	while (my $line = <$infile>) {
		chomp $line;
		my ($cluster, @genes) = split (/\t/, $line);
		$clusters->{$tag_name}->{$cluster} = \@genes;
		foreach my $g (@genes) {
			if (exists $data->{$tag_name}->{$g}) {
				print STDERR "$g found in more than one cluster in $tag_name\n";
				$data->{$tag_name}->{$g}->{$cluster} = 1;
			}else {
				$data->{$tag_name}->{$g}->{$cluster} = 1;
			}
		}
	}
	close ($infile);
}

