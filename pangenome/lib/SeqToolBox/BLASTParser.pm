# $Id: BLASTParser.pm 73 2007-10-24 19:29:03Z malay $
# Perl module for SeqToolBox::BLASTParser
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::BlastParser  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut



##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package SeqToolBox::BLASTParser;
use SeqToolBox::Root;
use IO::File;
use SeqToolBox::BLAST::BLASTPGP;
use Carp;
@ISA = qw(SeqToolBox::Root);
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
	my $obj = $self->_init(@_);
	return $obj;
}  


# _init is where the heavy stuff will happen when new is called

sub _init {
	my($self,@args) = @_;
	#print STDERR "@args\n";
	my ($file, $format, $handle) = $self->_rearrange(["FILE", "FORMAT","FH"], @args);
	#print STDERR "$file, $format\n";
	my $fh;
	
	if ($file) {
		$fh = IO::File->new ($file, "r");
	}
	
	if (defined $fh) {
		$self->{fh} = $fh;
	}elsif (defined $handle){
		$fh = $handle;
	}
	else{	croak "Could not open $file\n";
	}
	
	if ($format eq "PGP" || $format eq "pgp") {
		#print STDERR "*** PGP format\n";
		my $parser =  SeqToolBox::BLAST::BLASTPGP->new();
		$parser->set_fh ($fh);
		$parser->_parse_file();
		return $parser;
	}else {
		return $self;
	}
	
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
	#print STDERR "Gi list of blastparser\n";
}




=head1 PRIVATE METHODS

=cut

=head2 DESTROY()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub DESTROY {
	my ($self,@args) = @_;
	if (defined $self->{fh}){$self->{fh}->close();}
}

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2007 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;