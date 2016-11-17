#$Id: SeqToolBox.pm 395 2009-06-12 19:45:06Z malay $
# Perl module for SeqToolBox
# Author: Malay < malay@bioinformatics.org >
# Copyright (c) 2006 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

=head1 NAME

SeqToolBox - Malay's sequence manipulation toolbox.

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here


=cut


package SeqToolBox;
#@ISA = qw(SeqToolBox::Root);
@EXPORT_OK = qw();

use SeqToolBox::Root;
use File::Spec;
use strict;

my $DB_DIR = File::Spec->tmpdir();
if ($ENV{SEQTOOLBOXDB}) {
		$DB_DIR = $ENV{SEQTOOLBOXDB};
}


=head1 CONSTRUCTOR

=cut

sub new {
   my $class = shift;

   my $self = {};
   bless $self, ref($class) || $class;
 #  $self->_init(@_);
	$self->{db_dir} = $DB_DIR;
   return $self;

}

sub _init {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
return $make; # success - we hope!
}


=head1 METHODS

=cut

=head2 get_dbdir()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_dbdir {
	my ($self,@args) = @_;
	return $self->{db_dir};
}


=head1 SEE ALSO


=head1 COPYRIGHTS

Copyright (c) 2006 by Malay <malay@bioinformatics.org>. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

1;
