#$Id: Root.pm 3 2007-07-19 21:10:06Z malay $
# Perl module for SeqToolBox::Root
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2006 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

=head1 NAME

SeqToolBox::Root - Supplies some common method to all the classes.

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

package SeqToolBox::Root;    
use strict;

=head1 METHODS

=head2 _rearrange()

Describe your function here

 Usage   : This functions is adapted from BioPerl RootI.pm
 Args    : 
 Returns : 

=cut

sub _rearrange {
	my $self  = shift;
	my $order = shift;
	return @_ unless ( substr( $_[0] || '', 0, 1 ) eq '-' );
	push @_, undef unless $#_ % 2;
	my %param;
	while (@_) {
		( my $key = shift ) =~ tr/a-z\055/A-Z/d;
		$param{$key} = shift;
	}
	map { $_ = uc($_) } @{$order};
	return @param{@$order};
}

=head1 SEE ALSO


=head1 COPYRIGHTS

Copyright (c) 2006 by Malay <malay@bioinformatics.org>. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

1;
