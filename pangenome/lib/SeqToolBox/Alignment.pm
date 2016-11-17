# $Id: Alignment.pm 162 2008-02-15 20:57:45Z malay $
# Perl module for SeqToolBox::Alignment
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2008 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

SeqToolBox::Alignment  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package SeqToolBox::Alignment;

use vars qw(@ISA);

#@ISA       = qw(Bio::SimpleAlign);
@ISA       = qw (SeqToolBox::SeqDB SeqToolBox::Root );
@EXPORT_OK = qw();
use strict;
use SeqToolBox::Root;
use SeqToolBox::SeqDB;
use SeqToolBox::Seq;
use Carp;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

#sub new {
#	my $class = shift;
#	my $self = {};
#	bless $self, ref($class) || $class;
#
#	#my @seq;
##	my ($file, $format) = $self->_rearrange(["FILE", "FORMAT"], @_);
##	unless ($file) {
##		croak "$file is a mandetory option\n";
##	}
##
##	$format = 'fasta' unless $format;
##
##
#
##	my $seqdb = SeqToolBox::SeqDB->new(-file => $file, -format =>$format);
##	while (my $s = $seqdb->next_seq) {
##		push @seq, $s;
##	}
##	$self->{_seqs} = \@seq;
#
#	#$self->_init(@_);
#	$self->SUPER::new(@_);
#	#return $alignobj;
#	return $self;
#}

# _init is where the heavy stuff will happen when new is called

#sub _init {
#	my ( $self, @args ) = @_;
#	my $make = $self->SUPER::_initialize;
#	return $make;
#}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

=head2 is_flush()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub is_flush {
	my ( $self, @args ) = @_;
	if ( exists $self->{_is_flush} ) {
		return $self->{_is_flush};
	}
	my $length;
	while ( my $s = $self->next_seq ) {
		unless ($length) {
			$length = $s->length;
			next;
		}
		unless ( $length == $s->length ) {
			$self->{_is_flush} = 0;
			return 0;
		}
	}
	if ($length) {
		$self->{_length}   = $length;
		$self->{_is_flush} = 1;
		return 1;
	}
	else {
		croak "No sequence found\n";
	}
}

=head2 length()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub length {
	my ( $self, @args ) = @_;
	if ( exists $self->{_length} ) {
		return $self->{_length};
	}
	else {
		if ( $self->is_flush ) {
			return $self->{_length};
		}
		else {
			croak "Alignment is not flushed\n";
		}
	}
}

=head2 get_conserved_positions()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_ungapped_positions {
	my ( $self, @args ) = @_;
	my $length = $self->length();
	my @con_pos;

	if ( $self->is_flush() ) {
		for ( my $i = 0 ; $i < $length ; $i++ ) {

			#print STDERR "iter $i\n";
			my $aln = $self->slice( $i + 1, $i + 1 );

			#my $seq = "";
			#my $count = 0;
			my $gap_found = 0;

			while ( my $s = $aln->next_seq ) {

				#print STDERR "SEQCOUNT " ,$count++, "\n\n";
				#print $s->get_seq(), "\n";
				if ( $s->get_seq() eq "-" ) {
					$gap_found = 1;
					last;
				}
			}

			#print STDERR $i, "\t", $seq,"\n";
			#return $seq;
			unless ($gap_found) {
				push @con_pos, $i + 1;
			}

		}
		return @con_pos;

		#	foreach my $seq ( $self->each_seq) {
		#		print STDERR $seq->seq;
		#	}
	}
	else {
		croak "Alignment is not flush\n";
	}
}

=head2 slice() 

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub slice {
	my ( $self, $start, $end ) = @_;
	my $new_align = SeqToolBox::Alignment->new();

	while ( my $seq = $self->next_seq ) {
		my $string = $seq->get_seq();
		my $s      = $start - 1;

		my $substring = substr( $string, $start - 1, ( $end - $start ) + 1 );
		my $new_seq = SeqToolBox::Seq->new(
											-id   => $seq->get_id,
											-desc => $seq->get_desc,
											-seq  => $substring
		);
		$new_align->add_seq($new_seq);
	}
	return $new_align;
}

=head2 get_column_as_string()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_column_as_string {
	my ( $self, $pos ) = @_;
	if ( $pos > $self->length || $pos < 1 ) {
		croak "Position out of range";
	}
	my $seq;

	if ( $self->is_flush() ) {
		my $aln = $self->slice( $pos, $pos );
		while ( my $s = $aln->next_seq() ) {
			$seq .= $s->get_seq();
		}

#		#print STDERR $i, "\t", $seq,"\n";
#		#return $seq;
#		unless ($gap_found) {
#			push @con_pos, $i + 1;
#		}

		#}
		return $seq;

		#	foreach my $seq ( $self->each_seq) {
		#		print STDERR $seq->seq;
		#	}
	}
	else {
		croak "Alignment is not flush\n";
	}

}

sub get_conserved_postions {
	my ( $self, @args ) = @_;
	my $length = $self->length();
	my @con_pos;

	if ( $self->is_flush() ) {
		for ( my $i = 0 ; $i < $length ; $i++ ) {

			#print STDERR "iter $i\n";
			my $aln = $self->slice( $i + 1, $i + 1 );

			#my $seq = "";
			#my $count = 0;
#			my $gap_found = 0;
			my $last_aa;
			my $conserved = 1;
			
			while ( my $s = $aln->next_seq ) {

				#print STDERR "SEQCOUNT " ,$count++, "\n\n";
			
				#print $s->get_seq(), "\n";
				if ($s->get_seq() eq '-' ) {
					$conserved = 0;
					last;
				}
				
				if (!$last_aa ) {
					$last_aa = $s->get_seq();
					next;
				}
				
				if ( $s->get_seq() ne $last_aa ) {
					$conserved = 0;
					last;
				}
			}

			#print STDERR $i, "\t", $seq,"\n";
			#return $seq;
			if ($conserved) {
				push @con_pos, $i + 1;
			}

		}
		return @con_pos;

		#	foreach my $seq ( $self->each_seq) {
		#		print STDERR $seq->seq;
		#	}
	}
	else {
		croak "Alignment is not flush\n";
	}
}
=head1 PRIVATE METHODS

=cut

=head1 SEE ALSO

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 COPYRIGHTS

Copyright (c) 2008 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;
