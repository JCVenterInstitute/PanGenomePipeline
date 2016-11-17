# $Id$
# Perl module for Pangenome::ParserFactory
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2010 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::ParserFactory  - DESCRIPTION of Object

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


package Pangenome::ParserFactory;

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
	return $self->_init(@_);
	
}  


# _init is where the heavy stuff will happen when new is called

sub _init {
	my($self,@args) = @_;
	my $method_name = $args[0];
	my $real_method_name;
	if ($method_name =~ /inparanoid/i) {
		$real_method_name = "Inparanoid";
	}elsif ($method_name =~ /orthomcl/i) {
		$real_method_name = "Orthomcl";
	}elsif ($method_name =~ /panoct/i) {
	    $real_method_name = "PanOCT";
	}elsif ($method_name =~ /sybil/i) {
	    $real_method_name = "Sybil";
	}else {
		croak "Method name $method_name is not a valid option\n";
	}
#	print STDERR "$real_method_name\n";
	my $module_name = "Pangenome::Method::".$real_method_name."::Parser";
#	print STDERR "$module_name\n";
	eval "require $module_name";
	if ($@) {
		croak "Could not load module for method $method_name\n";
	} 
	my $object = $module_name->new();
	unless ($object) {
		croak "Could not create object from method $method_name\n";
	}
	return $object;
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut





=head1 PRIVATE METHODS

=cut




=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2010 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;
