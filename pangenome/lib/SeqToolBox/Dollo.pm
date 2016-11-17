# $Id: Dollo.pm 376 2009-03-24 22:33:22Z malay $
# Perl module for SeqToolBox::Dollo
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2008 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::Dollo  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package SeqToolBox::Dollo;

use vars qw(@ISA);
@ISA = qw(SeqToolBox::Root);
use SeqToolBox::Root;
use Bio::TreeIO;
use Bio::Tree::Node;
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
	my ($file) = $self->_rearrange(['FILE'], @args);
	croak "Filename not given" unless $file;
	$self->_parse_file($file); 
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

sub get_tree {
	my $self = shift;
	if (exists $self->{tree}) {
		return $self->{tree};
	}else{
		return;
	}
}



=head1 PRIVATE METHODS

=cut

sub _parse_file {
	my ($self, $file) = @_;
#	print STDERR "Parse file called\n";
	open (FILE, $file) || croak "Couldn't open $file\n";
	my $inside = 0;
	my $tree;
	while (my $line = <FILE>) {
		chomp $line;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
#		print "Line $line\n";
		if (!defined($line) || $line eq ""){
			next;
		}
#		next unless defined($line) || $line eq "";
#		print "Line2 $line\n";
		
		if ($line =~ /^root/) {
			$inside = 1;
		}
		
		if ($inside){
#			chomp $line;
			#print STDERR $line, "\n";
			my ($first, $second, $step, @f) = split (/\s+/, $line);
#			print STDERR "$first \t $second @f\n";
			my @data = split(//,join("", @f));
			$self->{nodedata}->{$second} = \@data;
#			print STDERR "@data\n";
			if ($first eq 'root') {
				my $rootnode = Bio::Tree::Node->new(-id=> $second );
				#$rootnode->add_Descendent(Bio::Tree::Node->new(-id =>$second));
				$tree = Bio::Tree::Tree->new(-root => $rootnode);
			}else{
				my $node = Bio::Tree::Node->new(-id => $second);
#				my $target_id = $tree->findnode_by_id ($first);
				my @target_ids = $tree->find_node($first);
				croak "target node not found" unless $target_ids[0];
				my $target_id = $target_ids[0];
				$target_id->add_Descendent($node);
			}
		}
	}
	close (FILE);
	$self->{tree} = $tree;
}


sub get_node_data_by_node {
	my ($self, $node) = @_;
	my $id = $node->id();
	croak "Node id not found\n" unless $id;
	if (exists $self->{nodedata}->{$id}) {
		return @{$self->{nodedata}->{$id}};
	}else {
		croak "Node data not found for $id\n";
	}
}

sub get_node_data_by_id {
	my($self, $id) = @_;
	if (exists $self->{nodedata}->{$id}) {
		return @{$self->{nodedata}->{$id}};
	}else {
		croak "Node data not found for $id\n";
	}
}

sub get_node_gain_by_id {
	my ($self, $id) = @_;
	my @node_data = $self->get_node_data_by_id($id);
	my $gain = 0;
	foreach my $i (@node_data) {
		if ($i eq "1") {
			$gain++;
		}
	}
	return $gain;
}

sub get_node_loss_by_id {
	my ($self, $id) = @_;
	my @node_data = $self->get_node_data_by_id($id);
	my $loss = 0;
	foreach my $i (@node_data) {
		if ($i eq "0") {
			$loss++;
		}
	}
	return $loss;
}


=head1 SEE ALSO

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 COPYRIGHTS

Copyright (c) 2008 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;