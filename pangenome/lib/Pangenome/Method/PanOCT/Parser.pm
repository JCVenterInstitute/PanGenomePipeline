# $Id$
# Perl module for Pangenome::Method::PanOCT::Parser
# Author: Jason Inman <jinman@jcvi.org>
# Copyright (c) 2010 by JCVI. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::Method::PanOCT::Parser  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jason Inman <jinman@jcvi.org>


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package Pangenome::Method::PanOCT::Parser;

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

}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

sub parse_file {
	my ($self, $file) = @_;
	$self->{file} = $file;
}

sub dump {

	my $self = shift;
	my $dumper_script = File::Spec->catfile ($FindBin::Bin,File::Spec->updir(),"analysis", "parse_panoct.pl");
	my $file = $self->{file};
	system ("perl $dumper_script $file") == 0 || croak "Can't run dumper script on $file\n";

}


=head1 PRIVATE METHODS

=cut

=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2010 by JCVI. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;
