# $Id$
# Perl module for Pangenome::Config
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2010 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::Config  - DESCRIPTION of Object

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


package Pangenome::Config;
use base ("Pangenome::Root");
use Config;
use Cwd;
use strict;
use FindBin;
use File::Spec;

use constant PATH => 'Inparanoid:mcl:MySQL:OrthoMCL:PanOCT:Sybil';



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


=head2 get_current_dir()

Returns the current directory.

	Usage   : my $current_dir = Pangenome::Config->get_current_dir();
  	Args    : None.
  	Returns : A scalar.
  	
=cut


sub get_current_dir {
	my ($self,@args) = @_;
	return cwd();
}

sub set_working_dir {
	my ($self, $wd) = @_;
	$self->{working_dir} = $wd;
}

sub get_working_dir {
	my $self = shift;
	if (exists $self->{working_dir}) {
		return $self->{working_dir};
	}else {
		return;
	}
}

sub set_starting_dir {
	my ($self, $sd) =@_;
	$self->{starting_dir} = $sd;
}

sub get_starting_dir {
	my $self = shift;
	if (exists $self->{starting_dir}) {
		return $self->{starting_dir};
	}else {
		return;
	}
}


sub find_bin {
    my ($self, $binary_name) = @_;
    my @f = split(/\-/,$Config{archname});
    my $arch = join('-', $f[0], $f[1]);
    my @paths = split ('\:',PATH);
   
    foreach my $p (@paths) {
#		my $p = '';
#		print STDERR $p, "\n";
	my @fullpaths; 
	push @fullpaths, File::Spec->catfile($FindBin::Bin,'..','ext',$p,$arch,$binary_name);
	push @fullpaths, File::Spec->catfile ($FindBin::Bin,'..','ext',$p,$binary_name);
	push @fullpaths, File::Spec->catfile ($FindBin::Bin,'..','ext',$p,'bin',$arch,$binary_name);
	push @fullpaths, File::Spec->catfile ($FindBin::Bin,'..','ext',$p,'bin',$binary_name);

	foreach my $fp (@fullpaths) {
	    if (-s $fp) {
		return $fp;
	    }
	}
    }
    return;
    
}

1;
