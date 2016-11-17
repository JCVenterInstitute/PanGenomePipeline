#!/usr/bin/env perl
# $Id: create_syntany_db.pl 32957 2010-08-05 16:22:09Z mbasu $

##---------------------------------------------------------------------------##
##  File: create_syntany_db.pl
##       
##  Author:
##        Malay <malay@bioinformatics.org>
##
##  Description:
##     
#******************************************************************************
#* Copyright (C) 2010 Malay K Basu <malay@bioinformatics.org> 
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

create_syntany_db.pl - One line description.

=head1 SYNOPSIS

create_syntany_db.pl [options] -o <option>


=head1 DESCRIPTION

Write a description of your prgram. 


=head1 ARGUMENTS 

=over 4

=item B<--option|-o>

First option.



=back

=head1 OPTIONS

Something here.


=head1 SEE ALSO

=head1 COPYRIGHT

Copyright (c) 2010 Malay K Basu <malay@bioinformatics.org>

=head1 AUTHORS

Malay K Basu <malay@bioinformatics.org>

=cut



##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use DBI;
use SeqToolBox::Interval::Tree;
use SeqToolBox::Interval::Range;

my $coord_genome_dir = shift;
opendir( DIR, $coord_genome_dir ) || die "Can't open $coord_genome_dir\n";
print STDERR "After opening\n";
my $coord;
while ( my $file = readdir(DIR) ) {
	print STDERR $file, "\n";
	next unless $file =~ /(\S+)\.coords/;
	my $basename = $1;
	$coord->{$basename}
		= read_coord( File::Spec->catfile( $coord_genome_dir, $file ) );
		print STDERR "KEYS: " , keys %{$coord}, "\n";

}
close DIR;
my $dbfile = "synteny.sqlite3";

my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "",
						 { AutoCommit => 0, RaiseError => 1 } );
$dbh->do(
		'create table synteny (
									organism text,
									node_name text,
									node_start numeric,
									node_end numeric,
									previous_node_name,
									next_node_name
									)'
	);
	
my $sth  = $dbh->prepare('insert into synteny values(?,?,?,?,?,?)');


foreach my $genome (keys %{$coord}) {
	my $root = $coord->{$genome}->getroot();
	while ($root) {
		my $name = $root->name();
		my $start = $root->start();
		my $end = $root->end();
		my $left;
		if ($root->left) {
			$left = $root->left()->name();
		}else {
			$left = 'NA';
		}
		my $right;
		if ($root->right){
			$right = $root->right()->name();
		}else {
			$right = 'NA';
		}
		$sth->execute($genome, $name, $start, $end, $left, $right);
		$root = $root->right();
	}
	
}
$sth->finish();
$dbh->do ('create index index1 on synteny (organism, node_name)');
$dbh->commit();
$sth = undef;
$dbh->disconnect();
exit 0;


sub read_coord {
	my $filename = shift;
	open( FILE, $filename ) || die "Can't open $filename\n";
	my $tree = SeqToolBox::Interval::Tree->new();

	while ( my $line = <FILE> ) {
		chomp $line;
		my ( $name, $start, $end, $dummy ) = split( /\t/, $line );

		if ( $start > $end ) {
			my $tmp = $start;
			$start = $end;
			$end   = $tmp;

		}
		my $range = SeqToolBox::Interval::Range->new( $name, $start, $end );
		$tree->insert($range);
	}

	return $tree;
}

##---------------------------------------------------------------------------##
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
##---------------------------------------------------------------------------##


