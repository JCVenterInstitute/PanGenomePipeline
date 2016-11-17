
=head1 NAME

ConsistencyChecks.pm - Module containing data integrity checks for the PanGenome pipeline

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

=head1 CONTACT

=cut

package Pangenome::ConsistencyChecks;
use strict;
use DBI;
use Data::Dumper;

sub new {
    my ( $class, $dbs, $server, $working_dir, $user, $password ) = @_;
    my $self = {};
    bless $self, ref($class) || $class;

    $self->{DBS}         = $dbs;
    $self->{SERVER}      = $server;
    $self->{WORKING_DIR} = $working_dir;
    $self->{USERNAME}    = $user;
    $self->{PASSWORD}    = $password;
    
    $self->{DBH} = DBI->connect( "dbi:Sybase:server=$server",$user, $password, { PrintError => 1 } );
    
    return $self;
}

sub check_loci_protein {
    my $self = shift;
    my $error;
    
    my $dbs    = $self->{DBS};
    my $db_missing_loci;
    my $db_missing_protein;
    
    foreach (sort @$dbs) {
	my($db,$asmbl_id) = split(/:/,$_);
	my @a_id = split(/,/,$asmbl_id); 
	my $a_id_string;
	map{$a_id_string .= "$_,"} @a_id;
	$a_id_string =~ s/,$//;
	
	my $dbh = $self->{DBH};
	$dbh->do("use $db");
	
	my $l_query = "SELECT i.locus,i.feat_name,a.protein_datalength " . 
	    "FROM ident i, asm_feature a, stan s " . 
	    "WHERE i.feat_name = a.feat_name " . 
	    "AND a.asmbl_id = s.asmbl_id " . 
	    "AND a.feat_type = \"ORF\" " . 
	    "AND i.locus is NULL ";
	
	my $p_query = "SELECT i.locus,i.feat_name,a.protein_datalength " . 
	    "FROM ident i, asm_feature a, stan s " . 
	    "WHERE i.feat_name = a.feat_name " . 
	    "AND a.asmbl_id = s.asmbl_id " . 
	    "AND a.feat_type = \"ORF\" " . 
	    "AND a.protein_datalength = 0 ".
	    "AND NOT EXISTS (select 1 from ORF_attribute o where o.feat_name=a.feat_name and o.att_type in (\"AFS\", \"APM\", \"DEG\")) ".
	    "AND NOT EXISTS (select 1 from role_link r where r.feat_name = a.feat_name and r.role_id = 270) ";
	
	if($asmbl_id eq "ISCURRENT"){
	    $l_query .= "AND s.iscurrent = 1";
	    $p_query .= "AND s.iscurrent = 1";
	}else{
	    $l_query .= "AND s.asmbl_id IN ($a_id_string)";
	    $p_query .= "AND s.asmbl_id IN ($a_id_string)";
	}
	
	my $sth = $dbh->prepare($l_query);
	$sth->execute();
	my ( $locus, $feat_name, $protein_length );
	$sth->bind_columns( \$locus, \$feat_name, \$protein_length);
	while ( $sth->fetch() ) {
	    $db_missing_loci->{$db} = 1;
	}
	$sth->finish;
	
	$sth = $dbh->prepare($p_query);
	$sth->execute();
	$sth->bind_columns( \$locus, \$feat_name, \$protein_length);
	while ( $sth->fetch() ) {
	    $db_missing_protein->{$db} = 1;
	}
    }
    
    if ( scalar keys %$db_missing_loci > 0 ) {
	$error .= "Features missing loci, please rerun locus loader for the following dbs:\n";
	$error .= join( "\n", keys %$db_missing_loci ) . "\n";
    }
    
    if ( scalar keys %$db_missing_protein > 0 ) {
	$error .= "Features missing protein, please rerun rewrite seqs for the following dbs:\n";
	$error .= join( "\n", keys %$db_missing_protein ) . "\n";
    }
    
    return $error if $error;
}

sub check_database_changes_since_data_pull {
    my $self = shift;
    my $error;
    
    #grab deflines from fasta file
    
    my $fasta_file = $self->{WORKING_DIR} . "/fastas/combined.fasta";
    die("Could not find combined.fasta") unless ( -s $fasta_file );
    
    my @deflines = `grep '>' $fasta_file`;
    my $fasta_loci;
    
    foreach my $defline (@deflines) {
	$defline =~ s/^>//;
	$defline =~ s/\s+$//;
	$fasta_loci->{$defline} = 1;
    }
    
    #pull loci from DB
    my $server = $self->{SERVER};
    my $db_loci;
    
    my $query = "SELECT  i.locus, f.end5, f.end3, i.com_name, f.asmbl_id, f.feat_name ".
	"FROM ident i, asm_feature f, stan s  ".
	"WHERE s.asmbl_id = f.asmbl_id  ".
	"AND f.feat_name = i.feat_name  ".
	"AND f.protein_datalength > 1 ";
    
    my $dbs = $self->{DBS};
    
    foreach (@$dbs) {
	my($db,$aid) = split(/:/,$_);
	my @a_id = split(/,/,$aid); 
	my $a_id_string;
	map{$a_id_string .= "$_,"} @a_id;
	$a_id_string =~ s/,$//;
	
	my $dbh = DBI->connect( "dbi:Sybase:server=$server", $self->{USERNAME}, $self->{PASSWORD}, { PrintError => 1 } );
	die("Can't connect to $server\n") unless defined $dbh;
	
	$dbh->do("use $db");
	
	if($aid eq "ISCURRENT"){
	    $query .= "AND s.iscurrent = 1 ";
	}else{
	    $query .= "AND s.asmbl_id IN ($a_id_string) ";
	}
	
	$query .= "ORDER BY f.asmbl_id, (f.end5 + f.end3)/2";
	
	my $sth = $dbh->prepare($query);
	$sth->execute();
	
	my ( $locus, $end5, $end3, $com_name, $asmbl_id, $feat_name );
	$sth->bind_columns( \$locus, \$end5, \$end3, \$com_name, \$asmbl_id, \$feat_name );
	while ( $sth->fetch() ) {
	    $db_loci->{$locus}->{DB} = $db if $locus;
	    $db_loci->{$feat_name}->{DB} = $db unless $locus;
	}
    }
    
    my ( @added, @deleted );
    
    foreach my $locus ( keys %$db_loci ) {
	unless ( exists $fasta_loci->{$locus} ) {
	    if ( exists $db_loci->{$locus} ) {
		push( @added, "$locus($db_loci->{$locus}->{DB})" );
	    }
	    else {
		push( @added, $locus );
	    }
	}
    }
    
    foreach my $locus ( keys %$fasta_loci ) {
	push( @deleted, $locus ) unless ( exists $db_loci->{$locus} );
    }
    
    $error .= "ERROR: Modifications have been made to the SGDs after initial pep files were pulled.";
    $error .= "If new features were added rerun the pep file generation and the blast searches.";
    $error .= "If features were deleted remove these features from the pep file and the blast results.\n";
    
    if ( scalar @deleted > 0 ) {
	$error .= "\nDELETED ORFs from SGD since blast was run\n";
	$error .= join( "\n", @deleted );
	$error .= "\n";
    }
    
    if ( scalar @added > 0 ) {
	$error .= "\nADDED ORFs to SGD since blast was run\n";
	$error .= join( "\n", @added );
	$error .= "\n";
    }
    
    return $error if $error;
}1;
