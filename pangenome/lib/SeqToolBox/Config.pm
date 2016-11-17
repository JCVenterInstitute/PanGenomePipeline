# $Id: Config.pm 376 2009-03-24 22:33:22Z malay $
# Perl module for SeqToolBox::Config
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::Config  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package SeqToolBox::Config;

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
	#$self->_init(@_);
	return $self;
}  





##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

=head2 get_file_stamp ()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_file_stamp  {
	my ($self,@args) = @_;
	my ($package, $filename, $line) = caller();
	#print STDERR $filename;
	#my $get_time_stamp = $self->get_time_stamp();
	my $s = '# Time: '. localtime()."\n";
	$s .= '# Program: '. $filename. "\n";
#	for (my $i = 0; $i < @args; $i++) {
#		$s .= '# Option '. $i + 1 ." :". $args[$i] ."\n";
#	}
	$s .= '# Options: '. join (" ", @args)."\n";
	
	#print STDERR $s,"\n";
	return $s;
}


=head1 PRIVATE METHODS

=cut




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
