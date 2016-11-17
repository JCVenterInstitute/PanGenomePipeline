# $Id: Range.pm 382 2009-04-10 21:10:17Z malay $
# Perl module for SeqToolBox::Interval::Range
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2009 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

SeqToolBox::Interval::Range  - DESCRIPTION of Object

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

package SeqToolBox::Interval::Range;

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
	my $args  = shift;
	my $class = ref $args || $args;
	my $self  = {};
	bless $self, $class;
	$self->_init(@_);
	return $self;

}

# _init is where the heavy stuff will happen when new is called

sub _init {
	my $self = shift;
	my ( $name, $start, $end ) = @_;
	$self->{name}   = $name || undef;
	$self->{left}   = undef;
	$self->{right}  = undef;
	$self->{start}  = $start || undef;
	$self->{end}    = $end || undef;
	$self->{isroot} = 0;
}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

sub right {
	my $self = shift;
	return $self->{right};
}

sub left {
	my $self = shift;
	return $self->{left};
}

sub start {
	my $self = shift;
	return $self->{start};

}

sub set_start {
	my ( $self, $start ) = @_;
	$self->{start} = $start;
}

sub set_end {
	my ( $self, $end ) = @_;
	$self->{end} = $end;
}

sub end {
	my $self = shift;
	return $self->{end};
}

sub name {
	my $self = shift;
	return $self->{name};
}

sub intersects {
	my ( $self, $node ) = @_;

	#   print $self->start(),"\n";
	#    print $self->end(), "\n";
	#    print $node->start(), "\n";
	#     print $node->end()   , "\n";

	if ((  $self->start() >= $node->start() && $self->start() <= $node->end()
		)
		|| (    $node->start() >= $self->start()
			 && $node->start() <= $self->end() )

		)
	{
		return 1;
	}
	else {
		return undef;
	}
}

sub proximalto {
	my ( $self, $node ) = @_;

	if ( !$self->intersects($node) && $self->start() < $node->start() ) {
		return 1;
	}
	else {
		return undef;
	}
}

sub merge {
	my ( $self, $node ) = @_;
	my $start;
	my $end;

	if ( $self->start() >= $node->start() ) {
		$start = $node->start();
	}
	else {
		$start = $self->start();
	}

	if ( $self->end() >= $node->end() ) {
		$end = $self->end();
	}
	else {
		$end = $node->end();
	}
	$self->{start} = $start;
	$self->{end}   = $end;

}

sub replace_with {
	print STDERR "replace with called\n";
	my ( $self, $node, $tree ) = @_;
	$node->{right} = $self->{right};
	$node->{left}  = $self->{left};

	if ( defined $node->{right} ) {
		print STDERR "Node right undef\n";
		$node->{right}->{left} = $node;
	}

	if ( $self->{isroot} ) {
		$node->setasroot($tree);
	}
	else {
		print STDERR "Not root\n";
		$node->{left}->{right} = $node;
		print STDERR $node->{left}->{right}->to_string(), "\n";
		print STDERR "after setting left\n";
	}

}

sub setasroot {
	my ( $self, $tree ) = @_;
	$self->{isroot} = 1;
	$tree->{root}   = $self;
}

sub isroot {
	my $self = shift;
	return $self->{isroot};
}

sub to_string {
	my $self = shift;
	return $self->{name} . '[' . $self->start . ',' . $self->end . ']';
}

=head1 PRIVATE METHODS

=cut

=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2009 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;
