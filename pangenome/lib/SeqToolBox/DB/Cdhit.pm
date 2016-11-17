# $Id: Cdhit.pm 382 2009-04-10 21:10:17Z malay $
# Perl module for SeqToolBox::DB::Cdhit
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2009 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

SeqToolBox::DB::Cdhit  - DESCRIPTION of Object

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

package SeqToolBox::DB::Cdhit;
use base "SeqToolBox::Root";

#use Rose::DB;
#use Rose::DB::Object;
use DBI;
use Carp;
use Data::Dumper;
use strict;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	return $self;
}

sub _init {
	my ( $self, @args ) = @_;
	my ($dbfile) = $self->_rearrange( ["DBFILE"], @args );
	croak "DBFILE is a required parameter in ", __PACKAGE__, "\n"
		unless $dbfile;
	croak "$dbfile not found" unless ( -s $dbfile );
	$self->{dbfile} = $dbfile;
	my $dbh
		= DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "",
						{ RaiseError => 1, AutoCommit => 1 } )
		or die $DBI::errstr;
	$self->{dbh} = $dbh;
	
	#	my $db = SeqToolBox::DB::CdhitDB->new( database => $dbfile );
	#	$self->{db} = $db;

#	$self->{cdhits} = SeqToolBox::DB::CdhitDB::CdhitObject::Manager->new(db=>$db);
}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

=head2 get_organisms()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_organisms {
	my $self = shift;
#	print STDERR "organisms called\n";
#	my @names = ();
#	my $cdhits = SeqToolBox::DB::CdhitDB::CdhitObject::Manager->get_cdhits(
#		db       => $self->{db},
#		distinct => 1,
#		select   => ['organism']
#
#	);
#	print STDERR $cdhits, "\n";
	my $sql = 'select distinct(organism) from cdhit';
	my $sth = $self->{dbh}->prepare_cached($sql);
	$sth->execute;
	my @array =  @{$sth->fetchall_arrayref};
#	print STDERR "@array\n";
#	print STDERR Dumper($array[0]);
	my @result = map { $_->[0]} @array;
#	print STDERR "@result\n";
##	print STDERR $result[0], "\n";
	$sth->finish();
	$sth = undef;
	return @result;	

}

=head1 PRIVATE METHODS

=cut

=head2 _get()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub _get {
	my ( $self, @args ) = @_;
}

sub DESTROY {
	my $self = shift;
	undef $self->{dbh};
}
=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2009 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;

#package SeqToolBox::DB::CdhitDB;
#use base qw(Rose::DB);
#__PACKAGE__->use_private_registry;
#__PACKAGE__->register_db( driver => 'sqlite' );
#1;
#
#package SeqToolBox::DB::CdhitDB::CdhitObject::Manager;
#use base ("Rose::DB::Object::Manager");
#sub object_class {'SeqToolBox::DB::CdhitDB::CdhitObject'}
#__PACKAGE__->make_manager_methods('cdhits');
#
#1;
#
#package SeqToolBox::DB::CdhitDB::CdhitObject;
#use base "Rose::DB::Object";
#__PACKAGE__->meta->setup(
#	table => "cdhit",
#
##	columns => [
##		qw (organism, query, cdd_id, domain, e_value, length, q_start, q_end, h_start, h_end)
##	]
#);
#__PACKAGE__->meta->auto_initialize;
#
#1;
