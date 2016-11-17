# $Id: Taxonomy.pm 566 2010-07-28 16:56:32Z malay $
# Perl module for SeqToolBox::Taxonomy
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

SeqToolBox::Taxonomy  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package SeqToolBox::Taxonomy;
use DBI;
use SeqToolBox;
use vars qw(@ISA);
@ISA       = qw();
@EXPORT_OK = qw();
use Carp;
use strict;
use File::Spec;

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

	#$self->_init(@_);
	#my $seqtoolbox_db = $ENV{SEQTOOLBOXDB} ? $ENV{SEQTOOLBOXDB} : ".";
	my $seqtoolbox_db = SeqToolBox->new()->get_dbdir();
	$self->{dbdir} = File::Spec->catdir( $seqtoolbox_db, "taxonomy" );
	my $gi_taxid_db
		= File::Spec->catfile( $self->{dbdir}, "gi_taxid_prot.db" );
	my $nodes_db = File::Spec->catfile( $self->{dbdir}, "nodes.db" );

	unless ( -s $gi_taxid_db || -s $nodes_db ) {
		croak
			'Taxonomy databases not found. Did you run "update_taxonomy.pl"?';
	}
	$self->{gi_taxid_db}     = $gi_taxid_db;
	$self->{nodes_db}        = $nodes_db;
	$self->{gi_taxid_handle} = DBI->connect( "dbi:SQLite:dbname=$gi_taxid_db",
							   "", "", { AutoCommit => 0, RaiseError => 1 } );
	$self->{gi_taxid_statement} = $self->{gi_taxid_handle}
		->prepare("select tax_id from gi_taxid_prot where gi = ?");
		
	$self->{nodes_db_handle} = 	my $dbh     = DBI->connect( "dbi:SQLite:dbname=$nodes_db", "", "",
							{ AutoCommit => 0, RaiseError => 1 } );
	$self->{nodes_db_statement} = $dbh->prepare(
					 "select parent_tax_id,rank from nodes where tax_id = ?");
	my %predefined = (
		"vertebrates" => 7742,
		"fungi"       => 4751,

		#				   "metazoa"       => 33208,
		"green_plants"  => 33090,
		"diplomonads"   => 5738,
		"eubacteria"    => 2,
		"cyanobacteria" => 1117,
		"archaea"       => 2157,
		"eumetazoa"     => 6072
	);
	$self->{common_ids} = \%predefined;
	$self->{defaults} = [ "vertebrates",  "fungi",
						  "green_plants", "cyanobacteria",
						  "eubacteria"
	];
	my %toplevels = ( 2157  => 1,
					  2     => 1,
					  2759  => 1,
					  12884 => 1,
					  10239 => 1,
					  28384 => 1,
					  12908 => 1
	);
	$self->{toplevels} = \%toplevels;

	return $self;
	my %c_cache;
	my %t_cache;
	$self->{c_cache} = \%c_cache;
	$self->{t_cache} = \%t_cache;
}

# _init is where the heavy stuff will happen when new is called

sub _init {
	my ( $self, @args ) = @_;
	my $make = $self->SUPER::_initialize;
	return $make;
}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

sub classify {
	my ( $self, $gi, @argv ) = @_;
	my %required;

	if ( !@argv ) {
		my @classes = @{ $self->{defaults} };

		foreach my $i (@classes) {
			$required{ $self->{common_ids}->{$i} } = $i;
		}
	}
	else {

		foreach my $i (@argv) {
			if ( exists $self->{common_ids}->{$i} ) {
				$required{ $self->{common_ids}->{$i} } = $i;
			}
			else {
				$required{$i} = $i;
			}
		}
	}

	#	my $gi_taxid = $self->{dbdir}.'/gi_taxid.db';
	my $gi_taxid = $self->{gi_taxid_db};
	my $dbh      = DBI->connect( "dbi:SQLite:dbname=$gi_taxid", "", "",
							{ RaiseError => 1 } );
	my $sth = $dbh->prepare("select tax_id from gi_taxid_prot where gi = ?");
	$sth->execute($gi);
	my $count = 0;
	my $taxid;

	while ( my @row = $sth->fetchrow_array ) {
		$count++;

		if ( $row[0] ) {
			$taxid = $row[0];
		}
	}
	$sth->finish;
	$sth = undef;
	$dbh->disconnect;
	croak "$gi returned more than one taxid\n" if $count > 1;

	#croak "$gi did not have a taxid\n" unless $taxid;
	if ( !$taxid ) {
		return;
	}

	#	my $nodesdb = $self->{dbdir}.'/nodes.db';
	my $nodesdb = $self->{nodes_db};
	$dbh = DBI->connect( "dbi:SQLite:dbname=$nodesdb", "", "",
						 { AutoCommit => 0, RaiseError => 1 } );
	$sth = $dbh->prepare("select parent_tax_id from nodes where tax_id = ?");
	my $class;

	while (1) {
		last unless $taxid;
		if ( exists $required{$taxid} ) { $class = $required{$taxid}; last; }
		if ( exists $self->{toplevels}->{$taxid} ) { last; }
		$sth->execute($taxid);
		my $count = 0;
		my $temp;

		while ( my @row = $sth->fetchrow_array ) {

			if ( $row[0] ) {
				$count++;
				$temp = $row[0];
			}
		}
		croak "more than one parent found for $taxid\n" if $count > 1;
		$taxid = $temp;
	}
	$sth->finish;
	$sth = undef;
	$dbh->disconnect;
	if   ($class) { return $class; }
	else          { return undef; }
}

sub classify_taxon {
	my ( $self, $gi, @argv ) = @_;
	my %required;
	my $taxid = $gi;

	if ( !@argv ) {
		my @classes = @{ $self->{defaults} };

		foreach my $i (@classes) {
			$required{ $self->{common_ids}->{$i} } = $i;
		}
	}
	else {

		foreach my $i (@argv) {
			if ( exists $self->{common_ids}->{$i} ) {
				$required{ $self->{common_ids}->{$i} } = $i;
			}
			else {
				$required{$i} = $i;
			}
		}
	}

	#	my $gi_taxid = $self->{dbdir}.'/gi_taxid.db';
	#my $gi_taxid = $self->{gi_taxid_db};
	#my $dbh      = DBI->connect( "dbi:SQLite:dbname=$gi_taxid", "", "",
#							{ RaiseError => 1 } );
#	my $sth = $dbh->prepare("select tax_id from gi_taxid_prot where gi = ?");
#	$sth->execute($gi);
#	my $count = 0;
#	my $taxid;
#
#	while ( my @row = $sth->fetchrow_array ) {
#		$count++;
#
#		if ( $row[0] ) {
#			$taxid = $row[0];
#		}
#	}
#	$sth->finish;
#	$sth = undef;
#	$dbh->disconnect;
#	croak "$gi returned more than one taxid\n" if $count > 1;
#
#	#croak "$gi did not have a taxid\n" unless $taxid;
#	if ( !$taxid ) {
#		return;
#	}

	#	my $nodesdb = $self->{dbdir}.'/nodes.db';
#	my $nodesdb = $self->{nodes_db};
#	$dbh = DBI->connect( "dbi:SQLite:dbname=$nodesdb", "", "",
#						 { AutoCommit => 0, RaiseError => 1 } );
#	$sth = $dbh->prepare("select parent_tax_id from nodes where tax_id = ?");
	my $class;

	while (1) {
		last unless $taxid;
		if ( exists $required{$taxid} ) { $class = $required{$taxid}; last; }
		if ( exists $self->{toplevels}->{$taxid} ) { last; }
		$self->{nodes_db_statement}->execute($taxid);
		my $count = 0;
		my $temp;

		while ( my @row = $self->{nodes_db_statement}->fetchrow_array ) {

			if ( $row[0] ) {
				$count++;
				$temp = $row[0];
			}
		}
		croak "more than one parent found for $taxid\n" if $count > 1;
		$taxid = $temp;
	}
#	$sth->finish;
#	$sth = undef;
#	$dbh->disconnect;
	if   ($class) { return $class; }
	else          { return undef; }
}

sub collapse_taxon {
	my ( $self, $id, $rank ) = @_;
	croak "ERROR: ID not defined" unless $id;
	croak "ERROR: RANK not defined" unless $rank;
#	my $nodesdb = $self->{nodes_db};
#	my $dbh     = DBI->connect( "dbi:SQLite:dbname=$nodesdb", "", "",
##							{ AutoCommit => 0, RaiseError => 1 } );
#	my $sth = $dbh->prepare(
#					 "select parent_tax_id,rank from nodes where tax_id = ?");
	
	if (exists $self->{c_cache}->{$id}->{$rank}) {
		return $self->{c_cache}->{$id}->{$rank};
	}
	my $class;
	my $taxid = $id;
	while (1) {
		last unless $taxid;

		#		if (exists $required{$taxid} ) {$class = $required{$taxid}; last;}

		if ( exists $self->{toplevels}->{$taxid} ) { last; }
		$self->{nodes_db_statement}->execute($taxid);
		my $count = 0;
		my $temp;
		my $temp_rank;

		while ( my @row = $self->{nodes_db_statement}->fetchrow_array ) {

			if ( $row[0] ) {
				$count++;
				$temp      = $row[0];
				$temp_rank = $row[1];

			}
		}
		croak "more than one parent found for $taxid\n" if $count > 1;
		
		
		if ( $rank =~ /^$temp_rank$/i ) {
			$class = $taxid;
			last;
		}
		$self->{c_cache}->{$id}->{$temp_rank} = $taxid;
		$taxid = $temp;
		
		#		last if ($ eq $rank);
	}
	
	$self->{c_cache}->{$id}->{$rank} = $class;
	
	if   ($class) { return $class; }
	else          { return undef; }

}

sub get_taxon {
	my ( $self, $gi ) = @_;
	croak "Gi required in get_taxon\n" unless $gi;
	my $id;

	if ( $gi =~ /gi\|(\d+)/i ) {
		$id = $1;
	}
	elsif ( $gi =~ /^(\d+)$/ ) {
		$id = $gi;
	}
	else {
		croak "Malformed GI: $gi\n";
	}

	if (exists $self->{t_cache}->{$id}) {
		return $self->{t_cache}->{$id};
	}
#	my $gi_taxid = $self->{gi_taxid_db};
#	my $dbh = DBI->connect_cached("dbi:SQLite:dbname=$gi_taxid", "","", {RaiseError=>1});
#	my $sth = $dbh->prepare_cached("select tax_id from gi_taxid_prot where gi = ?");
#	$sth->execute ($id);
	$self->{gi_taxid_statement}->execute($id);
	my $count = 0;
	my $taxid;

	while ( my @row = $self->{gi_taxid_statement}->fetchrow_array() ) {
		$count++;

		if ( $row[0] ) {
			$taxid = $row[0];
		}
	}

	#	$sth->finish;
	#	$sth = undef;
	#$dbh->disconnect;
	#print STDERR "$gi returned more than one taxid\n" if $count > 1;
	#croak "$gi did not have a taxid\n" unless $taxid;
	
	$self->{t_cache}->{$id} = $taxid;
	
	return $taxid;

}

=head1 PRIVATE METHODS

=cut

sub DESTROY {
	my $self = shift;
	if ($self->{gi_texid_statement}) {
	$self->{gi_taxid_statement}->finish();
	}
	$self->{gi_taxid_statement} = undef;
	eval {$self->{gi_taxid_handle}->disconnect();};
	eval{$self->{nodes_db_statement}->finish();};
	$self->{nodes_db_statement} = undef;
	eval {$self->{nodes_db_handle}->disconnect();};
}

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
