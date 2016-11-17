# $Id: RunI.pm 38688 2013-06-17 16:47:12Z ekelsey $
# Perl module for Pangenome::Method::RunI
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2010 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

Pangenome::Method::RunI  - Root interface of of all the Methods. 

=head1 SYNOPSIS

All the method modules should inherit from this.

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


package Pangenome::Method::RunI;
@ISA = qw(Pangenome::Root);
use Pangenome::Root;


use Carp;
use strict;


sub execute {
 	croak "This should never be called here\n";	
}


sub result_file {
	croak "This should never be called here\n";
}

sub set_fasta_files {
	my ($self, @files) = @_;
	$self->{files} = \@files;
}

sub get_fasta_files {
	my $self = shift;
	if (exists $self->{files}) {
		return @{$self->{files}}; 
	}else {
		return;
	}
}

sub get_base_fasta_files {
	my $self = shift;
	if (exists $self->{base_fasta_files}) {
		return @{$self->{base_fasta_files}};
	}else {
		my @basenames;
		foreach my $f ($self->get_fasta_files()) {
			my ($v, $d, $b) = File::Spec->splitpath($f);
			push @basenames, $b;
			
		}
		$self->{base_fasta_files} = \@basenames;
		return @basenames;
	}
}
sub set_blast_file {
    my ($self, $blast_file) = @_;
    $self->{blast_file} = $blast_file;	
}

sub get_blast_file {
    my $self = shift;
    if (exists $self->{blast_file}) {
	return $self->{blast_file};
    }else {
	return;
    }
}
sub set_percent_id{
    my ($self,$id) = @_;
    $self->{percent_id} = $id;
}

sub get_percent_id{
    my $self = shift;
    if (exists $self->{percent_id}) {
	return $self->{percent_id};
    }else {
	return;
    }
}
sub set_strict{
    my ($self,$level) = @_;
    
    if($level ne 'none'){
	$self->{strict} = ($level eq 'high') ? "M" : "Y";
    }elsif($level eq 'none'){
	$self->{strict} = "N";
    }
}

sub get_strict{
    my $self = shift;
    if (exists $self->{strict}) {
	return $self->{strict};
    }else {
	return;
    }
}
sub set_crib_file {
    my ($self, $crib_file) = @_;
    $self->{crib_file} = $crib_file;
}


sub get_crib_file {
	my $self = shift;
	if (exists $self->{crib_file}) {
		return $self->{crib_file};
	}else {
		return;
	}
}

sub set_gene_att_file {
    my ($self, $gene_att_file) = @_;
    $self->{gene_att_file} = $gene_att_file;   
}

sub get_gene_att_file {
	my $self = shift;
	if (exists $self->{gene_att_file}) {
		return $self->{gene_att_file};
	}else {
		return;
	}
}

sub set_config {
	my ($self, $config) = @_;
	unless (ref($config) eq "Pangenome::Config") {
		croak "Wrong object in set_config\n";
	}
	$self->{config} = $config;
}

sub get_config {
	my $self = shift;
	if (exists $self->{config}) {
		return $self->{config};
	}else {
		return;
	}
}

1;
