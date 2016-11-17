# $Id$
# Perl module for Pangenome::Cluster;
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2012 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

Pangenome::Cluster;  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package Pangenome::Cluster;
use Carp;
use vars qw(@ISA);
@ISA       = qw();
@EXPORT_OK = qw();
use strict;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	return $self;
}

# _init is where the heavy stuff will happen when new is called

sub _init {
	my ( $self, @args ) = @_;
	$self->{_id}      = shift(@args);
	$self->{_members} = \@args;

}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

sub get_members_as_array {
	my $self = shift;
	return @{ $self->{_members} };
}

sub get_members_as_hashref {
	my $self = shift;
	my %hash;
	foreach my $i ( @{ $self->{_members} } ) {
		$hash{$i} = 1;
	}
	return \%hash;
}

sub delete_member {
	my ( $self, $gi ) = @_;
	my @members = @{ $self->{_members} };
	my @new_members;
	foreach my $i (@members) {
		next if ( $gi eq $i );
		push @new_members, $i;
	}
	$self->{_members} = \@new_members;
}

sub get_id {
	my $self = shift;
	return $self->{_id};

}

sub get_avg_score_by_gene {
	my ( $self, $g, $blast ) = @_;
	my @genes       = $self->get_members_as_array();
	my $count       = 0;
	my $total_score = 0;
	foreach my $m (@genes) {
		next if ( $m eq $g );
		$count++;
		my $score = $blast->get_avg_score( $g, $m ) || 0;
		$total_score += $score;

	}
	if ( $count != ( scalar(@genes) - 1 ) ) {
		croak("$g more than once found\n");
	}elsif ($count == 0) {
		return 0;
	} 
	else {
		warn ("$g has score $total_score with count $count\n");
		return $total_score / $count;
	}
}

sub get_best_scoring_genes {
	my $self  = shift;
	my $blast = shift;
	my @genes = $self->get_members_as_array();
	my @scores;
	my @return;
	my $best_score = 0;
	foreach my $g (@genes) {
		my $score = $self->get_avg_score_by_gene( $g, $blast );
		push @scores, $score;
		$best_score = $score if ( $score > $best_score );
	}
	for ( my $i = 0 ; $i < @genes ; $i++ ) {
		if ( $scores[$i] == $best_score ) {
			push @return, $genes[$i];
		}
	}
	return @return;

}

sub set_centroid {
	my $self = shift;
	my $id   = shift;
	$self->{_centroid} = $id;
}

sub get_centroid {
	my $self = shift;
	my $id   = shift;
	return $self->{_centroid} ? $self->{_centroid} : undef;
}

sub is_inside_of {
	my $self = shift; 
	my $cluster = shift;
	my @given_genes = $cluster->get_members_as_array();
	my @self_genes = $self->get_members_as_array ();
	if (@self_genes > @given_genes) {
		return 0;
	}
	my $hash = $cluster->get_members_as_hashref();
	foreach my $gene (@self_genes) {
		return 0 unless exists $hash->{$gene};
	}
	return 1;
}

=head1 PRIVATE METHODS

=cut

=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2012 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;
