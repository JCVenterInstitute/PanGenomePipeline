#!/usr/local/bin/perl
use warnings;
use strict;
$|++;

=head1 NAME

create_contig_file.pl - create a multi-fasta of the non-pseudo contigs in the named sgd.

=head1 SYNOPSIS

    USAGE: create_contig_file.pl -D <dbname> || -d <db.list>
                                 -U <sybase user> 
                                 -P <sybase pass>

=head1 OPTIONS

B<--database, -D>   :   SYBASE sgd.

B<--db_list, -d>    :   list of sybase sgd names.

B<--user, -U>       :   SYBASE userid.

B<--pass, -P>       :   SYBASE password for the given userid.

=head1 DESCRIPTION

When supplied with the above credentials, will pull out all contig sequences for molecules whose
type is not like 'pseudo%' and dump them into a multi-fasta named 'dbname'.

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use DBI;
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use Pod::Usage;
use Text::Wrap qw(wrap $columns);
$columns = 80;

my %opts;
GetOptions( \%opts,
            'database|D=s',
            'db_list|d=s',
            'user|U=s',
            'pass|P=s',
            'help|h',
        ) || die "Problem getting options.\n";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $opts{help};

check_options( \%opts );

my $dbh = DBI->connect( "dbi:Sybase:server=SYBPROD", $opts{user}, $opts{pass}, { RaiseError => 1, PrintError => 1 } );

my $db_list = build_db_list( \%opts );

my $count = 0;
my $total = scalar( @{$db_list} );

for my $db ( @{$db_list} ) {

    $count++;
    print STDOUT "Working on $db. $count of $total\n";

    # get us a filehandle.
    open( my $ofh, '>', "./$db" ) || die "Can't open $db for writing: $!\n";
    select $ofh;

    $dbh->do( "USE $db" );

    # Set TEXTSIZE
    my $SELECT_MAX_DATALENGTH = 'SELECT MAX(a.sequence_datalength) FROM assembly a, asmbl_data d, stan s WHERE d.type not like "pseudo%" and d.id = s.asmbl_data_id and s.asmbl_id = a.asmbl_id';
    my ( $max_datalength ) = $dbh->selectrow_array( $SELECT_MAX_DATALENGTH );

    unless ( $max_datalength ) {
        print STDOUT "Couldn't get max datalength for $db!\n";
        next;
    }
    
    $dbh->do( "SET TEXTSIZE $max_datalength" );

    # Get sequences and push them onto the growing output file.
    my $SELECT_ASMBL_ID_AND_SEQ = "select a.asmbl_id, a.sequence from assembly a, asmbl_data d, stan s where d.type not like 'pseudo%' and d.id = s.asmbl_data_id and s.asmbl_id = a.asmbl_id";
    my $sth = $dbh->prepare( $SELECT_ASMBL_ID_AND_SEQ );
    $sth->execute();

    my ( $asmbl_id, $seq );
    $sth->bind_columns( \$asmbl_id, \$seq );

    while ( $sth->fetch ) {

        print wrap( '', '', (">$asmbl_id\n$seq\n") );

    }

}

exit(0);


sub build_db_list {

    my ( $opts ) = @_;
    my @db_list;

    if ( $opts->{database} ) {
        push @db_list, $opts->{database};
    } else {
        open( my $ifh, '<', $opts->{db_list} ) || die "Can't open $opts->{db_list}: $!\n";
        while (<$ifh>) {
            chomp;
            push @db_list, $_;
        }
    }

    return \@db_list

}


sub check_options {

    my ( $opts ) = @_;
    my $errors = '';

    $errors .= "Must supply a --database or --db_list\n"
        unless ( $opts->{database} || $opts->{db_list} );

    $errors .= "Must supply --database OR --db_list\n"
        if ( $opts->{database} && $opts->{db_list} );

    $errors .= "Must supply a --user\n"
        unless ( $opts->{user} );

    $errors .= "Must supply a --pass\n"
        unless ( $opts->{pass} );

    die $errors if $errors;

}
