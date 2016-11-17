# $Id: Root.pm 32957 2010-08-05 16:22:09Z mbasu $
# Perl module for Pangenome::Root
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2010 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

Pangenome::Root  - The Root class of Pangenome framework. All classes should inherit from this class.


=head1 SYNOPSIS

You are not supposed to use this directly. Always inherit from this class. 


=head1 DESCRIPTION

The sole purpose of the existance of this module is to provide some low
level functionality to the Pangenome package. This module is automatically called
by all the modules.


=head1 CONTACT

Malay <malay@jcvi.org>


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _.

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package Pangenome::Root;
use strict;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

# I liked this method from the bioperl modules. I just lifted it verbatim. Great works
# guys!!!!

=head2 _rearrange()

Stolen from Bioperl project. Originally written by Licoln Stein (?). This function provides the ability to supply the parameters with the style of Perl::Tk module like this: myfunction (-foo => $x, -bar => $y ). The purpose of this method is to provide a nice looking interface to method calling and return values. 

	Usage   : 
				myfunction (-foo=>$x, -bar =>$y);
				
				sub myfunction {
					
					my ($foo, $bar) = _rearrange (["FOO", "BAR"], @ARGV);
					# $foo is now $x and $bar is now $y.
					
				}
	
  	Args    : The function takes two arguments. The first one is the a array reference containing order in which all the named variables will be returned. The second one is the full arguments to the function. Don't forget to use shift if you are using the function as method in an object.
  	Returns : An array of the parsed values in the order given as the first parameter.
  	
=cut

sub _rearrange {

	my ( $self, $order, @param ) = @_;

	return unless @param;
	return @param unless ( defined( $param[0] ) && $param[0] =~ /^-/ );

	for ( my $i = 0; $i < @param; $i += 2 ) {
		$param[$i] =~ s/^\-//;
		$param[$i] =~ tr/a-z/A-Z/;
	}

	# Now we'll convert the @params variable into an associative array.
	local ($^W) = 0;    # prevent "odd number of elements" warning with -w.
	my (%param) = @param;

	my (@return_array);

	# What we intend to do is loop through the @{$order} variable,
	# and for each value, we use that as a key into our associative
	# array, pushing the value at that key onto our return array.
	my ($key);

	foreach $key ( @{$order} ) {
		my ($value) = $param{$key};
		delete $param{$key};
		push( @return_array, $value );
	}

#    print "\n_rearrange() after processing:\n";
#    my $i; for ($i=0;$i<@return_array;$i++) { printf "%20s => %s\n", ${$order}[$i], $return_array[$i]; } <STDIN>;

	return (@return_array);
}

=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2010 by Malay <mbasu@jcvi.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;

