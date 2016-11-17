# $Id$
# Perl module for Pangenome::Clusterset
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2012 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::Clusterset  - DESCRIPTION of Object

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


package Pangenome::Clusterset;

use vars qw(@ISA);
@ISA = qw();
@EXPORT_OK = qw();
use Carp;
use strict;


##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {   
	my $class = shift;
	my $self = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	return $self;
}  


# _init is where the heavy stuff will happen when new is called

sub _init {
	my($self,@args) = @_;
	$self->{method} = $args[0];
	$self->{gimap} = {};
	$self->{idmap} = {};
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

sub add_cluster {
	my ($self, $cluster) = @_;
	my $id = $cluster->get_id();
	
	if (exists $self->{idmap}->{$id}) {
		croak "Cluster id already exists\n";
	}else {
	   $self->{idmap}->{$id} = $cluster;
	}
	
	my @gis = $cluster->get_members_as_array();
	foreach my $g (@gis) {
		$self->{gimap}->{$g}->{$id} =1;
	}
}

sub get_cluster_id_by_gene {
	my ($self, $gene) = @_;
	if (exists $self->{gimap}->{$gene}) {
		my @clusters = keys %{$self->{gimap}->{$gene}};
		return @clusters;
	}else {
		return;
	}
}

sub get_cluster_by_gene {
	my ($self, $gene) = @_;
	my @cluster_ids = $self->get_cluster_id_by_gene($gene);
	my @return;
	foreach my $c (@cluster_ids) {
		push @return, $self->get_cluster_by_id($c);
		
	}
	return @return;
}

sub get_cluster_by_id {
	my ($self, $id) = @_;
	unless (exists $self->{idmap}->{$id}) {
		croak "Could not find cluster with id $id\n";
	}
	return $self->{idmap}->{$id};
}

sub resolve_conflict {
	my ($self, $blast) = @_;
	warn $self->{method}, "\n";
	foreach my $g (keys %{$self->{gimap}}) {
		my @clusters = $self->get_cluster_id_by_gene($g);
		next if (@clusters == 1);
		
		
		for my $i (@clusters) {
			for my $j (@clusters) {
				next if ($i eq $j);
				my $cluster1 = $self->get_cluster_by_id ($i);
				my $cluster2 = $self->get_cluster_by_id ($j);
				if ($cluster1->is_inside_of($cluster2)) {
					my $removal_id = $cluster1->get_id();
					my @removal_genes = $cluster1->get_members_as_array();
					foreach my $rg (@removal_genes) {
						delete $self->{gimap}->{$rg}->{$removal_id};
						
					}
					delete $self->{idmap}->{$removal_id};
					return 0;
				}
			}
		}
		
		my $best_cluster;
		my $best_score = 0;
		
		foreach my $c_id (@clusters) {
			my $cluster = $self->get_cluster_by_id ($c_id);
			my $score = $cluster->get_avg_score_by_gene($g, $blast) || 0;
			if ($score > $best_score) {
				$best_score = $score;
				$best_cluster = $c_id;
			}
		}
		warn "Best cluster for $g is $best_cluster\n";
		foreach my $c_id (@clusters) {
			next unless ($c_id eq $best_cluster);
			delete $self->{gimap}->{$g}->{$c_id};
			my $cluster = $self->get_cluster_by_id($c_id);
			$cluster->delete_member($g);
		}
		return 0;
	}
	return 1;
}

sub get_all_genes {
	my $self = shift;
	if (exists $self->{gimap}) {
		my @gis = keys %{$self->{gimap}};
		return \@gis;
	}else {
		return undef;
	}
}

sub get_all_clusters {
	my $self = shift;
	my @return;
	foreach my $id (keys %{$self->{idmap}}) {
		push @return, $self->{idmap}->{$id};
	}
	return @return;
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