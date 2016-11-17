#!/usr/local/bin/perl
=head1 NAME

make_gene_att_for_btm.pl - produce a gene_att style file for BTM given some db names

=head1 SYNOPSIS

USAGE: make_gene_att_for_btm.pl
            --db <sgd_database_name>
            --db_list <comma,seperated,list,of,sgd,database>
            --crib_file <path/to/file/containing/db.names>
            --user <sybase username>
            --password <sybase password>
            --att_file_prefix <prefix_for_att_file>


=head1 OPTIONS

Select one of -d, -D, and -c.

B<--db,-d>
    Sybase sgd database to connect to.

B<--db_list,-D>
    A comma-seperated list of sgd database names.

B<--crib_file,-c>
    A path for a file containing a single sgd database per line.

B<--att_file_prefix,-a>
    THe prefix for the generated att_file.  The att_file will be named <att_file_prefix>.att_file

B<--user,-U>
    User account with select, insert, and update privileges on the specified database.

B<--password,-P>
    Password for user account specified.

B<--server,-S>
    Optional.  Sybase server to connect to (default = SYBTIGR).

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script takes in a given set of sgd database names and provides an 'att_file' suitable for
use with the BTM clustering module of the pangenome pipeline.

=head1  CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use TIGR::FASTArecord;
use TIGR::FASTAwriter;
use DBI;
use Pod::Usage;
use SeqToolBox::SeqDB;
use Pangenome::DB::ProkSybase;
use Data::Dumper;

my %opts; 
my $DEBUG;
GetOptions(\%opts,  'db|d=s',
                    'db_list|D=s',
                    'crib_file|c=s',
                    'att_file_prefix|a=s',
                    'server|S=s',
                    'user|U=s',
                    'password|P=s',
                    'debug',
                    'help|h',
        ) || die "Error getting options: $!";

pod2usage( {-exitval => 0, -verbose => 2} ) if ($opts{help});

&check_opts;

my $db_list = &make_db_list;

my $dbh = DBI->connect("dbi:Sybase:server=$opts{server}",$opts{user},$opts{password},{PrintError=>1});
die("Can't connect to SYBTIGR as $opts{user}\n") unless defined $dbh;

my $pre_query = "SELECT  i.locus, f.end5, f.end3, i.com_name, f.asmbl_id, f.protein_datalength, f.feat_name
                FROM    ident i, asm_feature f, stan s 
                WHERE ";
                
my $post_query = " AND s.asmbl_id = f.asmbl_id 
                AND f.feat_name = i.feat_name 
                AND f.protein_datalength > 1 
                AND f.feat_type = \"ORF\"
                AND NOT EXISTS (select 1 from ORF_attribute o where o.feat_name=f.feat_name and o.att_type in (\"AFS\", \"APM\", \"DEG\"))
                AND NOT EXISTS (select 1 from role_link r where r.feat_name = f.feat_name and r.role_id = 270)
                ORDER BY f.asmbl_id, (f.end5 + f.end3)/2"; 
my ($afh,$rfh);
open ($afh,">$opts{att_file_prefix}.att_file") || die "Can't open $opts{att_file}.att_file: $!\n";
open ($rfh, ">role_id_lookup.txt") || die "Can't open role_id_lookup.txt: $!\n";

my $print_role_hsh;

foreach my $db (sort {$a cmp $b} keys %$db_list) {
   
    $dbh->do("use $db") || die("Could not connect to $db");
    my $a_id = $db_list->{$db}->{id};
    
    my $role_id_lookup = &create_role_id_lookup($db,$dbh,$a_id);
    my $hmm_lookup = &create_hmm_lookup($db,$dbh,$a_id);

    my $query = $pre_query;
    if($a_id eq 'ISCURRENT'){
	$query .= "s.iscurrent = 1";
    }else{
	$query .= "s.asmbl_id IN ($a_id)";
    }
    
    $query .= $post_query;

    my $sth = $dbh->prepare($query);

    $sth->execute();
    my ($locus, $end5, $end3, $com_name, $asmbl_id, $protein_datalength, $feat_name);
    $sth->bind_columns(\$locus, \$end5, \$end3, \$com_name, \$asmbl_id, \$protein_datalength, \$feat_name);

    while($sth->fetch()) {
	my $role_id;

	my $db_name;

	if(exists $db_list->{$db}->{name}){
	    $db_name = $db_list->{$db}->{name};
	}else{
	    $db_name = $db;
	}

	unless($db_name){
	    print Dumper($db_list);
	    print Dumper($db,$db_name);exit;
	}
	print $afh "$asmbl_id\t$locus\t$end5\t$end3\t$com_name\t$db_name\t$protein_datalength";

	#Print role_ids
	my $role_ids;
	if (exists $role_id_lookup->{$locus}){
	    my @role_ids = keys %{$role_id_lookup->{$locus}};
	    $role_ids = join("|",@role_ids);

	    foreach my $id(@role_ids){
		$print_role_hsh->{$id}->{main} = $role_id_lookup->{$locus}->{$id}->{mainrole};
		$print_role_hsh->{$id}->{subrole} = $role_id_lookup->{$locus}->{$id}->{sub1role};
	    }
	    

	}

	if($role_ids){
	    print $afh "\t$role_ids";
	}else{
	    print $afh "\t";
	}
	
	#Print hmms
	my $hmms;
	if (exists $hmm_lookup->{$locus}){
	    my @ids = keys %{$hmm_lookup->{$locus}};
	    $hmms = join("|",@ids);
	}

	if($hmms){
	    print $afh "\t$hmms\n";
	}else{
	    print $afh "\t\n";
	}
    }

}

#Add line to role_id_lookup that represents a none assigned role_id
#Using id=0
print $rfh "0\tNONE\tNONE\n";

foreach my $id (sort {$a <=> $b} keys %$print_role_hsh){
    print $rfh "$id\t$print_role_hsh->{$id}->{main}\t$print_role_hsh->{$id}->{subrole}\n";
}

exit(0);

######### SUBS ##########
sub create_role_id_lookup{
    my ($db,$dbh,$a_id) = @_;
    my $hsh;

    my $query = "select r.role_id,e.mainrole,e.sub1role, i.locus from $db..role_link r, egad..roles e, $db..asm_feature a, $db..stan s, $db..ident i "
        . "where s.asmbl_id = a.asmbl_id "
        . " and a.feat_name = i.feat_name "
        . "and i.feat_name = r.feat_name "
        . "and a.feat_type = \"ORF\" "
	. "and e.role_id = r.role_id ";

    if($a_id eq 'ISCURRENT'){
	$query .= "and s.iscurrent = 1";
    }else{
	$query .= "and s.iscurrent IN ($a_id)";
    }

    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $array_ref = $sth->fetchall_arrayref();

    foreach my $entry(@$array_ref){
	$hsh->{$entry->[3]}->{$entry->[0]}->{mainrole} = $entry->[1];
	$hsh->{$entry->[3]}->{$entry->[0]}->{sub1role} = $entry->[2];
    }

    return $hsh;
}
sub create_hmm_lookup{
    my ($db,$dbh,$a_id) = @_;
    my $hsh;

    my $query = "select e.accession, i.locus ".
	        "from ident i, stan s, asm_feature a, evidence e ".
		"where s.asmbl_id = a.asmbl_id ".
		"and a.feat_name = i.feat_name ".
		"and i.feat_name = e.feat_name ".
		"and e.ev_type=\"HMM3\" ";

    if($a_id eq 'ISCURRENT'){
	$query .= "and s.iscurrent = 1";
    }else{
	$query .= "and s.iscurrent IN ($a_id)";
    }

    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $array_ref = $sth->fetchall_arrayref();

    foreach my $entry(@$array_ref){
	$hsh->{$entry->[1]}->{$entry->[0]} = 1;
    }

    return $hsh;
}
sub make_db_list {

    my $db_list;

    if ($opts{db}) {
        $db_list->{$opts{db}}->{id} = "ISCURRENT";
	$db_list->{$opts{db}}->{name} = $opts{db};
    } elsif ($opts{crib_file}) {
        open (my $cfh, "<$opts{crib_file}") || die "Can't open crib_file: $opts{crib_file}: $!\n";
        while (<$cfh>) {
	    chomp;
            my ($name,$db,$id) = split(/\t/,$_);
	    my $sgd = ($db) ? $db : $name;

	    if($id){
		$db_list->{$sgd}->{id} = $id;
	    }else{
		$db_list->{$sgd}->{id} = "ISCURRENT";
	    }

	    $db_list->{$sgd}->{name} = $name;
	}
    } elsif($opts{db_list}) {
        my @db_list = split(',',$opts{db_list});
	foreach my $db (@db_list){
	    $db_list->{$db}->{id} = "ISCURRENT";
	    $db_list->{$db}->{name} = $db;
	}
    } else {
        die "Don't have any database to work with.\n";
    }

    return $db_list;

}

sub check_opts {

    my $errors .= '',

    $opts{server} = 'SYBPROD' unless $opts{server};
    $opts{user} = $opts{user};
    $opts{password} = $opts{password};
    $DEBUG = $opts{debug};

    unless ($opts{db} || $opts{db_list} || $opts{crib_file} ) {
        $errors .= "Need -d, -D, or -c\n";
    }

    $errors .= "Need a --att_file_prefix\n" unless $opts{att_file_prefix};

    die $errors if $errors;

}
