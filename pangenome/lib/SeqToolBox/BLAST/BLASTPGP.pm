# $Id: BLASTPGP.pm 14 2007-07-26 19:51:53Z malay $
# Perl module for SeqToolBox::BLAST::BLASTPGP
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::BLAST::BLASTPGP  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package SeqToolBox::BLAST::BLASTPGP;
use SeqToolBox::Root;
use SeqToolBox::BLASTParser;

use vars qw(@ISA);
@ISA = qw(SeqToolBox::Root SeqToolBox::BLASTParser);

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
	
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

=head2 get_gi_list()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_gi_list {
	my ($self,@args) = @_;
	if (exists $self->{gi_list}) {
		my @temp = @{$self->{gi_list}};
		return @temp;
	}else{
		return;
	}
}



=head1 PRIVATE METHODS

=cut

=head2 set_fh()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub set_fh {
	my ($self,@args) = @_;
	$self->{fh} = $args[0];
}

=head2 parse_file()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub _parse_file {
	my $self = shift;
	my $fh = $self->{fh};
	while (my $line = <$fh>) {
		last if ($line =~ /Sequences producing significant alignments\:/);
	}
	while (my $line = <$fh>){
		last if $line =~ /Sequences with E-value WORSE than threshold/;
		if ($line =~ /gi\|(\d+)\|/) {
			$self->_add_to_gi_list($1);
		}else {
			
		}
	}
	
	
}

=head2 _add_to_gi_list()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub _add_to_gi_list {
	my ($self,$gi) = @_;
	if (exists $self->{gi_list}) {
		my @temp = @{$self->{gi_list}};
		push @temp , $gi;
		$self->{gi_list} = \@temp;
	}else{
		my @temp;
		push @temp , $gi;
		$self->{gi_list} = \@temp;
	}
}

=head1 SEE ALSO

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 COPYRIGHTS

Copyright (c) 2007 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;