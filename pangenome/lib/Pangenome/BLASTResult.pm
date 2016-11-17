# $Id$
# Perl module for Pangenome::BLASTResult
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2012 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::BLASTResult  - DESCRIPTION of Object

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


package Pangenome::BLASTResult;

use vars qw(@ISA);
@ISA = qw();
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
	my $self = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	$self->{_hits} = {};
	return $self;
}  


# _init is where the heavy stuff will happen when new is called

sub _init {
	my($self,@args) = @_;
	
	 
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

sub get_raw_score {
	my ($self, $q, $s) = @_;
	if (exists $self->{_hits}->{$q}->{$s}) {
		return $self->{_hits}->{$q}->{$s};
	}else {
		return;
	}
}

sub add_hit {
	my ($self, $q, $s, $score) = @_;
	if (exists $self->{_hits}->{$q}->{$s}) {
		if ($self->{_hits}->{$q}->{$s} < $score) {
			$self->{_hits}->{$q}->{$s} = $score;			
		}
	}else {
		$self->{_hits}->{$q}->{$s} = $score;
	}
}

sub get_avg_score {
	my ($self, $q, $s) = @_;
	my $score1 = 0;
	my $score2 = 0;
	if (exists $self->{_hits}->{$q}->{$s}) {
		$score1 =  $self->{_hits}->{$q}->{$s};
	}
	
	if (exists $self->{_hits}->{$s}->{$q}) {
		$score2 = $self->{_hits}->{$s}->{$q};
		
	}
	
	return ($score1 + $score2) /2;
}

sub get_normalized_raw_score {
	my ($self, $q, $s) = @_;
	my $self_score = $self->get_raw_score($q, $q) || 0;
	my $score = $self->get_raw_score($q, $s) || 0;
	return $self_score? $score/$self_score:undef;
	
	
}

sub get_normalized_avg_score {
	my ($self, $q, $s) = @_;
	my $self_score1 = $self->get_normalized_raw_score($q, $s) || 0;
	my $self_score2 = $self->get_normalized_raw_score ($s, $q) || 0;
	if ($self_score1 && $self_score2) {
		return ($self_score1 + $self_score2)/2;
	}else {
		return;
	}
}

sub get_avg_distance_score {
	my ($self, $g, $others) = @_;
	my $count = 0;
	my $tot_score = 0;
	foreach my $o (@{$others}) {
		next if ($o eq $g);
		my $score = $self->get_normalized_avg_score($g, $o);
		$count++;
		$tot_score += $score;
	}
	if ($count > 0) {
		return $tot_score/$count;
	}else {
		return 0;
	}
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