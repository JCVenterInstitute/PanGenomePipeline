# $Id: Set.pm 585 2010-11-15 23:38:48Z malay $
# Perl module for SeqToolBox::Tools::Set
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::Tools::Set  - DESCRIPTION of Object

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


package SeqToolBox::Tools::Set;

use vars qw(@ISA);
@ISA = qw(SeqToolBox::Root);
@EXPORT_OK = qw();
use strict;
use SeqToolBox::Root;

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
#	$self->_init(@_);
	$self->{sets} = ();
	return $self;
}  



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

=head2 add_set()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub add_set {
	my ($self,@args) = @_;
	if ($self->{sets}){
		my @temp = @{$self->{sets}};
		push @temp, @args;
		$self->{sets} = \@temp; 
	}else {
		$self->{sets} = \@args;
	}
}

=head2 get_set()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_set {
	my ($self,$num) = @_;
	#print STDERR $self->{sets}->[0];
	if (defined $self->{sets}->[$num]){
		return @{$self->{sets}->[$num]};
	}else{
		return;
	}
}

=head2 get_all_union()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_all_union {
	my ($self) = @_;
	my %union;
	foreach my $i (@{$self->{sets}}){
		foreach my $j (@{$i}){
			$union{$j} = 1;
		}
	}
	my @temp = sort {$a cmp $b} keys %union;
	#print STDERR "@temp\n";
	return @temp;
}

##-------------------------------------------------------------------------##

=head2 get_all_intersection()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_all_intersection {
	my ($self,@args) = @_;
	my @intersection;
	my $length = scalar (@{$self->{sets}});
	
	foreach my $i (@{$self->{sets}}) {
		foreach my $j (@{$i}){
			my $found = 0;
			foreach my $k (@{$self->{sets}}){
				if ($self->_exists_in ($j, $k)){
					$found++;
				}
			}
			if ($found == $length){
				push @intersection, $j;
			}
		}
		last;	
	}
	
	return sort @intersection;	

}

=head2 get_all_intersection_as_hash()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_all_intersection_as_hash {
	my ($self,@args) = @_;
	my %intersection;
	my $length = scalar (@{$self->{sets}});
	
	foreach my $i (@{$self->{sets}}) {
		foreach my $j (@{$i}){
			my $found = 0;
			foreach my $k (@{$self->{sets}}){
				if ($self->_exists_in ($j, $k)){
					$found++;
				}
			}
			if ($found == $length){
				$intersection{$j}= 1;
			}
		}
		last;	
	}
	
	return \%intersection;	
	
	
}


##-------------------------------------------------------------------------##

=head2 get_intersection(){}

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_intersection {
	my ($self,@args) = @_;
	my @set1 = @{$self->{sets}->[$args[0]]};
	my @set2 = @{$self->{sets}->[$args[1]]};
	my @intersection;
	foreach my $i (@set1){
		if ($self->_exists_in($i, \@set2) ){
			push @intersection, $i;
		}
	}
	return sort @intersection;
}




=head1 PRIVATE METHODS

=cut

##-------------------------------------------------------------------------##

=head2 _exists_in()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub _exists_in {
	my ($self,@args) = @_;
	my $value = $args[0];
	my @array = @{$args[1]};
	#print STDERR "value= $value\t @array\n";
	
	foreach my $i (@array) {
		if ($i eq $value){
			return 1;
		}
	}
	return;
}

=head2 get_unique()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_uniques {
	my ($self,@args) = @_;
	my $intersection = $self->get_all_intersection_as_hash();
	my $result = ();
	for (my $i = 0; $i < @{$self->{sets}}; $i++) {
		for (my $j = 0; $j < @{$self->{sets}->[$i]}; $j++) {
			if (exists $intersection->{$self->{sets}->[$i]->[$j]}) {
				next;
			}else {
				push @{$result->[$i]}, $self->{sets}->[$i]->[$j];
			}
		}
	}
	
	return $result;	
}


=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2007 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;