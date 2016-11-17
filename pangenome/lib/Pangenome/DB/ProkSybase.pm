package Pangenome::DB::ProkSybase;

BEGIN {
	eval { use DBD::Sybase };

	if ($@) {
		die "Could not load DBD::Sybase, make sure it is installed\n";

	}
}

use strict;
use warnings;
use DBI;
use Data::Dumper;
use Carp;

sub new {
	my ( $class, $params ) = @_;

	my $self = {};
	bless( $self, $class );

	$self->_init($params);

	return $self;
}

sub connect {
	my ( $self, $user, $pass, $database, $server ) = @_;
	my $retval = 1;

	my %vars = ( 'user'     => \$user,
				 'pass'     => \$pass,
				 'database' => \$database,
				 'server'   => \$server,
	);

	foreach my $var ( keys %vars ) {

		unless ( ${ $vars{$var} } ) {
			unless ( $self->{$var}  || $var eq 'database') {
				warn("$var should be set or passed in before connecting")
					unless ( $self->{$var} );
				$retval = 0;
				last;
			} else {
				${ $vars{$var} } = $self->{$var};
			}
		}
	}

	if ($retval) {
		$self->{'dbh'} = DBI->connect(
								   "dbi:Sybase:server=$server; packetSize=8092",
								   $user, $pass,
								   { 'RaiseError' => 1,
									 'AutoCommit' => 0
								   }
		) or die( "Database connection not made: " . DBI::errstr );
		$self->{'dbh'}->do("use $database") if $database;
		$retval = $self->{'dbh'};
	}
	$self->{'dbh'}->{syb_chained_txn} = 1;

	return $retval;
}

sub commit {
	my $self = shift;
	$self->{dbh}->commit();

}

sub rollback {
	my $self = shift;
	$self->{dbh}->rollback();

}

sub begin_work {
	my $self = shift;
	$self->{dbh}->begin_work();
}
######### HERE WE GET THINGS FROM THE DATABASE ##############################

sub get_organism_name {
	my ($self,$db) = @_;
	#if ( $self->{organism_name} ) { return $self->{organism_name}; }
	my $query    = "select organism_name from $db..new_project";
	my $arrayref = $self->_do_select_array($query);

	#	print STDERR Dumper($arrayref);
	if ( $arrayref->[0]->[0] ) {
		$self->{organism_name} = $arrayref->[0]->[0];
		return $arrayref->[0]->[0];
	} else {
		croak "Could not get organism name";
	}

	#	if (exists $hashref->{organism_name}) {
	#		print STDERR $hashref->{organism_name};
	#	}

}

sub get_sgd_meta_data{
	my ($self,$db) = @_;
	my $results;
	
	my $query = "select attribute_id, attribute_type, description from $db..metadata_attribute";
	my $arrayref = $self->_do_select_array($query);
	
	if(scalar @$arrayref > 0){
		foreach my $row (@$arrayref){   	
		   	$results->{$row->[0]}->{'type'} = $row->[1];
		   	$results->{$row->[0]}->{'description'} = $row->[2];
		}
		
		$query = "select attribute_id, attribute_value from $db..metadata";
		$arrayref = $self->_do_select_array($query);
	    
	    foreach my $row (@$arrayref){       
	        $results->{$row->[0]}->{'value'} = $row->[1];
	    }
	}
	
    return($results)
}

sub pan_get_meta_data{
	my ($self, $id) = @_;
	
	my $query = "SELECT attribute_id from metadata where genomes_data_id = $id";
	my $hashref = $self->_do_select_query( $query, 'attribute_id' );

    return $hashref;
}

sub pan_find_singleton_inserts{
	my $self = shift;
	my $method_id = $self->{method_id};
	
	my $query = "SELECT cluster_id,cluster_name ".
	            "FROM cluster ".
	            "where method_id=$method_id ".
	            "and cluster_name like \"S%\"";
	            
	my $hashref = $self->_do_select_query( $query, 'cluster_name' );
    return $hashref;
}
sub get_all_loci{
	my $self = shift;  
    my $results;
    
    my $query = "select a.feat_name, i.locus, a.asmbl_id ".
                "FROM asm_feature a, ident i, stan s ".
                "WHERE a.asmbl_id=s.asmbl_id ".
                "AND a.feat_name = i.feat_name ".
                "AND s.iscurrent = 1 ".
                "AND a.feat_type = \"ORF\"";
               
               
    my $hashref = $self->_do_select_query( $query, 'feat_name' );
    return $hashref;
}


sub pan_populate_meta_data{
    my ($self, $data, $genome_id, $bcp) = @_;	

    my $query = "SELECT attribute_type from metadata_attribute";
    my $arrayref = $self->_do_select_array($query);
    my $db_types = {};
    
    foreach my $row ( @$arrayref) {
        $db_types->{$row->[0]} = 1;
    } 

    my $id;
      
    foreach my $sgd_id (keys %$data){
    	my $type = $data->{$sgd_id}->{'type'};
        my $description = $data->{$sgd_id}->{'description'};
        
    	#Only load attributes that don't currently exists in metadata_attribute
    	unless(exists $db_types->{$type}){
		    my $i1 = $self->_do_insert_query("metadata_attribute", {attribute_type=>$type, 
	         	                                                   description=>$description});
    	}
    	
        my $query = "SELECT attribute_id from metadata_attribute where attribute_type = \"$type\"";      
        my $idref = $self->_do_select_array($query);
        
        $data->{$sgd_id}->{'value'} =~ s/(\"|\;)//g;
        
        my $insert_query = "INSERT into metadata ".
                           "(genomes_data_id, attribute_id, attribute_value) ".
                           "VALUES ($genome_id, $idref->[0][0], \"$data->{$sgd_id}->{'value'}\")";
        
        my $dbh = $self->{'dbh'};
        my $result;  
        eval {
	        my $sth = $dbh->prepare($insert_query);
	        $sth->execute;
	        $sth = $dbh->prepare('select @@identity');
	        $sth->execute;
	        my $row = $sth->fetchrow_arrayref();
            $result = $row->[0];
        };

	    if ($@) {
	        $self->{'dbh'}->rollback();
	        print STDERR "Could not insert data into method $@" . DBI->errstr;
	    }	                                          
        	                                                   
    }
}

sub pan_populate_genomes_data {
	my ( $self, $data ) = @_;
	my $d = $self->_do_insert_query( "genomes_data", $data );

	#	$self->{dbh}->commit()
	return $d ? $d : undef;

}

sub pan_get_genomes_data {
	my ( $self, $data ) = @_;
	
	if (exists $self->{genomes_data_id}) {
		#return $self->{genomes_data_id};
	}
	my $arrayref =
		$self->_do_conditional_select( 'genomes_data', 'genomes_data_id',
									   $data );

	if ( $arrayref->[0]->[0] ) {
		#$self->{genomes_data_id} = $arrayref->[0]->[0];	
		return $arrayref->[0]->[0];
	} else {
		return;
	}
}

#sub pan_get_members_id {
#	my ( $self, $db, $locus ) = @_;
#	my $dbh   = $self->{dbh};
#	my $query = "select member_id from members where db = ? and locus = ?";
#	my $sth   = $dbh->prepare_cached($query);
#	$sth->execute( $db, $locus );
#	my $array = $sth->fetchrow_arrayref();
#	$sth->finish();
#
#	if ( $array->[0] ) {
#		return $array->[0];
#	} else {
#		return;
#	}
#}

sub pan_populate_method {
	my ( $self, $data ) = @_;
	my $first_half  = "INSERT INTO method (";
	my $second_half = "values ( ";

	foreach my $key ( keys( %{$data} ) ) {
		$first_half .= " $key,";

		if ( $key eq 'type' || $key eq 'type_version' ) {
			$second_half .= "'" . $data->{$key} . "'" . ", ";
		} else {
			$second_half .= $data->{$key} . ", ";
		}

	}

	$first_half  =~ s/,$//;
	$second_half =~ s/,\s$//;

	my $query = $first_half . " ) " . $second_half . " )";

	my $dbh = $self->{'dbh'};

	#	print STDERR $query, "\n";
	my $result;
	eval {

		my $sth = $dbh->prepare($query);
		$sth->execute;
		$sth = $dbh->prepare('select @@identity');
		$sth->execute;
		my $row = $sth->fetchrow_arrayref();
#		print STDERR Dumper($row);
#		print STDERR $row->[0], "\n";
		$result = $row->[0];

		#	while (my $row = $sth->fetchrow_arrayref()) {
		#		print STDERR "IDENTITY value = $row->[0]\n";
		#	}
		#	$sth->finish();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		print STDERR "Could not insert data into method " . DBI->errstr;
	}
	
	$self->{method_id} = $result;
	return $result;

	#	my $query;
	#	my $d = $self->_do_insert_query ("method", $data);
	#	return $d ? $d : undef;
}
sub get_loci_feat_db_lookup{
	my ($self,$db) = @_;
	
	my $query = "SELECT a.feat_name, i.locus,a.asmbl_id ".
	            "from $db..asm_feature a, $db..ident i, $db..stan s ".
		    "where i.feat_name=a.feat_name ".
		    "and a.asmbl_id = s.asmbl_id ".
		    "and a.feat_type = \"ORF\" ".
		    "and s.iscurrent =1";

    my $hashref = $self->_do_select_query( $query, 'locus');
    return $hashref;
}
sub pan_populate_cluster_members {
	my ( $self, $data ) = @_;
	my $first_half  = "INSERT INTO cluster_members (";
	my $second_half = "values ( ";

	foreach my $key ( keys( %{$data} ) ) {
		$first_half .= " $key,";

		if ( $key eq 'cluster_name' ) {
			$second_half .= "'" . $data->{$key} . "'" . ", ";
		} else {
			$second_half .= $data->{$key} . ", ";
		}

	}

	$first_half  =~ s/,$//;
	$second_half =~ s/,\s$//;

	my $query = $first_half . " ) " . $second_half . " )";

	my $dbh = $self->{'dbh'};

	#	print STDERR $query, "\n";
#	my $result;
	eval {
		my $sth = $dbh->prepare($query);
		$sth->execute;
#		$sth = $dbh->prepare('select @@identity');
#		$sth->execute;
#		my $row = $sth->fetchrow_arrayref();
#		print STDERR Dumper($row);
#		print STDERR $row->[0], "\n";
#		$result = $row->[0];

		#	while (my $row = $sth->fetchrow_arrayref()) {
		#		print STDERR "IDENTITY value = $row->[0]\n";
		#	}
		#	$sth->finish();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		print STDERR "Could not insert data into method " . DBI->errstr;
		return;
	}
	return 1;
}

sub pan_populate_cluster_table {
	my ( $self, $data ) = @_;
	my $first_half  = "INSERT INTO cluster (";
	my $second_half = "values ( ";

	foreach my $key ( keys( %{$data} ) ) {
		$first_half .= " $key,";

		if ( $key eq 'cluster_name' ) {
			$second_half .= "'" . $data->{$key} . "'" . ", ";
		} else {
			$second_half .= $data->{$key} . ", ";
		}

	}

	$first_half  =~ s/,$//;
	$second_half =~ s/,\s$//;

	my $query = $first_half . " ) " . $second_half . " )";

	my $dbh = $self->{'dbh'};

	#	print STDERR $query, "\n";
	my $result;
	eval {
		my $sth = $dbh->prepare($query);
		$sth->execute;
		$sth = $dbh->prepare('select @@identity');
		$sth->execute;
		my $row = $sth->fetchrow_arrayref();
#		print STDERR Dumper($row);
#		print STDERR $row->[0], "\n";
		$result = $row->[0];

		#	while (my $row = $sth->fetchrow_arrayref()) {
		#		print STDERR "IDENTITY value = $row->[0]\n";
		#	}
		#	$sth->finish();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		print STDERR "Could not insert data into method " . DBI->errstr;
	}
	return $result;
	
}


sub pan_get_member_id {
	my ($self,$data) = @_;
	my $arrayref = $self->_do_conditional_select('members', 'member_id', $data);
	if ( $arrayref->[0]->[0] ) {
#		$self->{organism_name} = $arrayref->[0]->[0];
		return $arrayref->[0]->[0];
	} else {
		return;
	}
	
	
}
sub pan_populate_member_table {
	my ($self, $data) = @_;
	my $first_half  = "INSERT INTO members (";
	my $second_half = "values ( ";

	foreach my $key ( keys( %{$data} ) ) {
		$first_half .= " $key,";

		if ( $key eq 'feat_name' || $key eq 'locus' || $key eq 'assignby') {
			$second_half .= "'" . $data->{$key} . "'" . ", ";
		} else {
			$second_half .= $data->{$key} . ", ";
		}

	}

	$first_half  =~ s/,$//;
	$second_half =~ s/,\s$//;

	my $query = $first_half . " ) " . $second_half . " )";

	my $dbh = $self->{'dbh'};

	#print STDERR $query, "\n";
	
	my $result;
	eval {
		my $sth = $dbh->prepare($query);
		$sth->execute;
		$sth = $dbh->prepare('select @@identity');
		$sth->execute;
		my $row = $sth->fetchrow_arrayref();
#		print STDERR Dumper($row);
#		print STDERR $row->[0], "\n";
		$result = $row->[0];

		#	while (my $row = $sth->fetchrow_arrayref()) {
		#		print STDERR "IDENTITY value = $row->[0]\n";
		#	}
		#	$sth->finish();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		print STDERR "Could not insert data into method " . DBI->errstr;
	}
	return $result;
}
sub pan_get_strain_data{
    my $self = shift;
    
    my $query = "SELECT genomes_data_id, db FROM genomes_data";	
	my $hashref = $self->_do_select_query( $query, 'genomes_data_id' );
    return $hashref;
}

sub pan_get_member_data{
	my $self = shift;
	
	my $query = "SELECT m.member_id, m.genomes_data_id, m.locus, m.feat_name, c.cluster_id from members m, cluster_members c where c.member_id=m.member_id";
	my $hashref = $self->_do_select_query( $query, 'member_id' );
    return $hashref;
}
sub pan_get_member_count_by_genome_id{
    my ($self,$databases) = @_;
    
    my $dbh = $self->{dbh};
    
    my $query = "SELECT count(m.member_id) as members, m.genomes_data_id ".
		        "FROM members m ".
				"GROUP BY m.genomes_data_id";
				
	my $hashref = $self->_do_select_query( $query, 'genomes_data_id' );
    return $hashref;    
}

sub pan_get_method_id{
    my ($self,$method) = @_;
    my $query = "SELECT method_id,date from method where type=\"$method\""; 	
	my $arrayref = $self->_do_select_array( $query );
	
	$self->{method_id} = $arrayref->[0][0];
	$self->{date} = $arrayref->[0][1];
}
sub pan_find_clusters_with_member_count{
    my ($self,$count) = @_;
      
    my $method_id = $self->{method_id};
    my $dbh = $self->{dbh};
    my $comparison = ($count != 0) ? "= $count" : "> 1";
    
    my $query = "SELECT m.cluster_id,m.member_id ". 
				"FROM cluster_members m, cluster c ".
				"WHERE c.cluster_id=m.cluster_id ".
				"AND c.method_id=$method_id ".			
				"GROUP BY m.cluster_id HAVING count(m.member_id) $comparison";
	

	my $cluster_ids;
	my $sth = $dbh->prepare($query);
    $sth->execute;
    
	while (my ($id,$member_id) = $sth->fetchrow_array()){
        $cluster_ids->{$id} .= "$member_id,"
    }
    
    return $cluster_ids;    	
}

sub pan_get_members_per_cluster{
    my $self = shift;
    
    my $method_id = $self->{method_id};
    my $dbh = $self->{dbh};
    
    my $query = "SELECT m.cluster_id,s.locus, s.genomes_data_id,s.member_id ". 
                "FROM cluster_members m, cluster c, members s ".
                "WHERE c.cluster_id=m.cluster_id ".
                "AND c.method_id=$method_id ".
                "AND s.member_id=m.member_id";    

    my ($cluster_ids,$member_ids);
    
    my $sth = $dbh->prepare($query);
    $sth->execute;
    
    while (my ($id,$locus,$genome_id,$member_id) = $sth->fetchrow_array()){
    	$cluster_ids->{$id}->{$genome_id} .= "$locus,";
    	$member_ids->{$locus} = $member_id;
    }
   

    return ($cluster_ids,$member_ids);        
}
sub pan_get_cluster_ids_per_method_id{
    my ($self,$method_id) = @_;
    
    my $query = "SELECT c.cluster_id, c.cluster_name,s.stat_value,i.com_name,i.gene_sym,l.accession
					FROM cluster c, cluster_stat s, cluster_ident i, cluster_link l
					where c.method_id = $method_id
					and c.cluster_id = s.cluster_id
					and c.cluster_id = i.cluster_id
					and i.cluster_ident_id = l.cluster_ident_id
					and l.type=\"role_id\"";
	my $hashref = $self->_do_select_query( $query, 'cluster_id' );
    return $hashref;
}

sub pan_get_databases_for_cluster{
    my ($self,$members) = @_;
    
    my @results;
    
    my $query = "SELECT genomes_data_id,locus,member_id from members where member_id IN ($members)";
    my $hashref = $self->_do_select_query( $query, 'genomes_data_id' );
    return $hashref;  
}
sub get_frameshift_by_genome{
	my $self = shift;
	
	my $query = "SELECT i.locus, f.att_type ".
	            "FROM frameshift f, ident i, asm_feature a, stan s ".
	            "WHERE s.asmbl_id = a.asmbl_id ".
	            "AND s.iscurrent=1 ".
	            "AND a.feat_name=i.feat_name ".
	            "AND i.feat_name=f.feat_name";
	            
    my $hashref = $self->_do_select_query( $query, 'locus' );
    return $hashref;
}
sub get_feat_name_asm_id_by_locus {
	my ( $self, $locus,$db ) = @_;
	my $dbh = $self->{dbh};
	my $query
		= "select a.feat_name, a.asmbl_id from ident as i,asm_feature as a where a.feat_name = i.feat_name and i.locus = ?";
	my $sth = $dbh->prepare_cached($query);
	$sth->execute($locus);
	my $array = $sth->fetchrow_arrayref();
	$sth->finish();
	return $array->[0], $array->[1];

}

sub get_gene_sequences {
	my ( $self, $filter ) = @_;

	die( __PACKAGE__ . "::get_gene_sequences not yet implemented" );
}

sub get_user {
	my ($self) = @_;
	return $self->{'user'};
}

sub get_password {
	my ($self) = @_;
	return $self->{'pass'};
}

sub get_genes {
	my ( $self, $filter, $force ) = @_;

	if ( $self->{'genes'} ) {
		return $self->{'genes'} unless ($force);
	}

	my $query
		= "SELECT af.feat_id, af.feat_name, af.sequence, af.protein, af.end5, af.end3 "
		. "FROM asm_feature af, stan s, assembly a "
		. "WHERE af.feat_type = 'ORF'"
		. "AND af.asmbl_id = a.asmbl_id "
		. "AND a.asmbl_id = s.asmbl_id "
		. "AND s.iscurrent = 1"
		. "AND datalength(af.protein) <> 0";

	my $hashref = $self->_do_select_query( $query, 'feat_name' );
	$self->{'genes'} = $hashref;
	return $hashref;

}

sub get_gene_coords {
	my ( $self, $feat_name ) = @_;
	my $retval = [];

	unless ( defined( $self->{'genes'} ) ) {
		$self->get_genes();
	}

	unless ( exists( $self->{'genes'}->{$feat_name} ) ) {
		die("Could not parse out coordinates for feature: $feat_name");
	}

	$retval->[0] = $self->{'genes'}->{$feat_name}->{'end5'};
	$retval->[1] = $self->{'genes'}->{$feat_name}->{'end3'};
	return $retval;

}

sub get_genes_within_range {
	my ( $self, $left, $right ) = @_;

	my $query
		= "SELECT af.feat_name, af.end5, af.end3 "
		. "FROM asm_feature af, stan s "
		. "WHERE (af.end5 BETWEEN $left AND $right "
		. "OR af.end3 BETWEEN $left AND $right) "
		. "AND af.feat_type = 'ORF' "
		. "AND af.asmbl_id = s.asmbl_id "
		. "AND s.iscurrent = 1";

	my $hashref = $self->_do_select_query( $query, 'feat_name' );

	return $hashref;
}

sub get_protein {
    my ( $self, $db, $asmbl_id, $feat_name ) = @_;

    my $query = "SELECT protein FROM $db..asm_feature WHERE asmbl_id=$asmbl_id AND feat_name='$feat_name'";
    my $arrayref = $self->_do_select_array($query);
    return($arrayref->[0][0]);
}

sub get_iscurrent_annotation {
    my ( $self, $db ) = @_;

    my $query
        = "SELECT i.locus,a.feat_name, i.com_name,i.gene_sym,i.ec# "
        . "FROM $db..ident i, $db..asm_feature a, $db..stan s "
        . "WHERE s.iscurrent = 1 "
        . "AND a.asmbl_id = s.asmbl_id "
        . "AND a.feat_name = i.feat_name "
        . "AND a.feat_type = \"ORF\"";

    my $hashref = $self->_do_select_query( $query, 'locus');
    return $hashref;
}

sub get_role_ids_for_genome{
	my ($self,$db) = @_;
	
	my $query
        = "select r.role_id, i.locus from $db..role_link r, $db..asm_feature a, $db..stan s, $db..ident i "
        . "where s.asmbl_id = a.asmbl_id "
        . " and a.feat_name = i.feat_name "
        . "and i.feat_name = r.feat_name "
        . "and s.iscurrent = 1 "
        . "and a.feat_type = \"ORF\"";

    my $hashref = $self->_do_select_query( $query, 'locus');
    return $hashref;

}

sub get_protein_db {
	my ($self) = @_;

	if ( $self->{'protein_db'} ) {
		return $self->{'protein_db'};
	}

	unless ( $self->{'genes'} ) {
		$self->get_genes;
	}

	foreach my $feat_name ( keys %{ $self->{'genes'} } ) {
		$self->{'protein_db'}->{$feat_name}
			= $self->{'genes'}->{$feat_name}->{'protein'};
	}

	return $self->{'protein_db'};

}

sub get_gene_accessions {
	my ($self) = @_;

	unless ( defined( $self->{'genes'} ) ) {
		$self->get_genes;
	}

	return keys %{ $self->{'genes'} };
}

sub get_db_name {
	my ($self) = @_;
	return $self->{'database'};
}
######################## HERE WE PUT THINGS IN THE DATABASE ########################
sub store_annotation_rule_hit {
	my ( $self, $annot_rule, $orf_acc, $data ) = @_;
	my %already_loaded;

	my $rule_id = $annot_rule->get_id;
	foreach my $method ( $annot_rule->get_methods ) {
		next if ( exists( $already_loaded{$orf_acc}->{$rule_id} ) );
		$method->store_hit( $data, $self, $orf_acc, $annot_rule->{'id'} );
		$already_loaded{$orf_acc}->{$rule_id} = 1;
	}
}

sub insert_evidence {
	my ( $self, $values ) = @_;

	my $insert_id;

	my $query
		= "INSERT INTO evidence (feat_name, ev_type, accession, date, assignby, "
		. "end5, end3, rel_end5, rel_end3, m_lend, m_rend, curated, change_log, save_history, method ) "
		. "values( ?, ?, ?, getdate(), ?, 0, 0, ?, ?, ?, ?, 1, 1, 1, 'annotation_rule' )";

	my $dbh = $self->{'dbh'};

	eval {
		my $sth = $dbh->prepare($query);

		#This will check the values
		foreach my $column (
				qw/feat_name ev_type accession rel_end5 rel_end3 m_lend m_rend/)
		{
			die("Field $column is required to insert data into the evidence table"
				)
				unless ( exists( $values->{$column} )
						 && defined( $values->{$column} ) );
		}

		$values->{'assignby'} = $self->get_user;

		$sth->execute( $values->{'feat_name'}, $values->{'ev_type'},
					   $values->{'accession'}, $values->{'assignby'},
					   $values->{'rel_end5'},  $values->{'rel_end3'},
					   $values->{'m_lend'},    $values->{'m_rend'}
		);
		$dbh->commit();
	};

	if ($@) {
		warn "Could not insert evidence: " . DBI->errstr;
		$dbh->rollback;
		return;
	}

	#Now find the evidence we just inserted
	my $squery
		= "SELECT top 1 id from evidence where "
		. "feat_name = '"
		. $values->{'feat_name'}
		. "' and "
		. "ev_type = '"
		. $values->{'ev_type'}
		. "' and "
		. "accession = '"
		. $values->{'accession'}
		. "' and "
		. "assignby = '"
		. $values->{'assignby'}
		. "' and "
		. "rel_end5 = "
		. $values->{'rel_end5'} . " and "
		. "rel_end3 = "
		. $values->{'rel_end3'} . " and "
		. "m_lend = "
		. $values->{'m_lend'} . " and "
		. "m_rend = "
		. $values->{'m_rend'};

	my $id_hash = $self->_do_select_query( $squery, 'id' );

	($insert_id) = keys( %{$id_hash} );
	return $insert_id;
}

sub insert_feat_score {
	my ( $self, $values ) = @_;

	my $query = "INSERT INTO feat_score ( input_id, score_id, score )"
		. "values( ?, ?, ? )";

	my $dbh = $self->{'dbh'};

	eval {
		my $sth = $dbh->prepare($query);

		foreach my $column (qw/input_id score_id score/) {
			die("Field $column is required to insert data into the evidence table"
			) unless ( exists( $values->{$column} ) );
		}

		$sth->execute( $values->{'input_id'}, $values->{'score_id'},
					   $values->{'score'} );
		$dbh->commit();

	};

	if ($@) {
		warn "Could not insert into feat_score: " . DBI->errstr;
		$dbh->rollback;
	}

}
############################# REMOVE FROM DATABASE #############################

sub remove_all_rule_base {
	my ($self) = @_;

	my $select_ev_ids = "SELECT id from common..score_type where input_type = "
		. "'RULE_BASE'";

	my $id_ref = $self->_do_select_query( $select_ev_ids, 'id' );
	my ($score_id) = keys %{$id_ref};

	die("Could not get the score id from common..score_type")
		unless ( defined($score_id) );

	#Remove from evidence
	my $del_query = "DELETE from evidence where ev_type = 'RULE_BASE'";

	my $dbh = $self->{'dbh'};

	eval {

		my $dth = $dbh->prepare($del_query);
		$dth->execute
			or die(   "Could not remove past annotation rule data from "
					. "evidence [$del_query]"
					. DBI->errstr );

		#Remove from feat_score
		$del_query = "DELETE from feat_score where score_id = $score_id";

		$dth = $dbh->prepare($del_query);
		$dth->execute
			or die(   "Could not remove past annotation rule data from "
					. "evidence [$del_query]"
					. DBI->errstr );

		$dbh->commit();
	};

	if ($@) {
		warn "Could not delete all rule_base evidence: "
			. DBI->errstr
			. ". Rolling back";
		$dbh->rollback();
	}
}

sub remove_rule_base {
	my ( $self, $rule_ids ) = @_;

	if ( $rule_ids && $rule_ids eq 'all' ) {
		$self->remove_all_rule_base();
		return;
	}

	#Find evidence ids with these rule ids
	my @accs = keys %{$rule_ids};

	local $" = "', '";
	my $query
		= "SELECT id "
		. "FROM evidence "
		. "WHERE accession in ( '@accs' ) "
		. "AND ev_type = 'RULE_BASE'";

	my $hash_ref = $self->_do_select_query( $query, 'id' );
	my @ev_ids = keys(%$hash_ref);

  #Return and don't delete anything if we didn't find any evidence ids to delete
	return if ( @ev_ids == 0 );

	#Find the rule base score_type id
	$query = "SELECT id from common..score_type "
		. "WHERE input_type = 'RULE_BASE'";

	$hash_ref = $self->_do_select_query( $query, 'id' );
	my ($rule_base_id) = keys %$hash_ref;

	#delete everything
	eval {
		$" = ", ";
		$query
			= "DELETE from evidence "
			. "WHERE id in ( @ev_ids ) "
			. "AND ev_type = 'RULE_BASE'";

		my $dth = $self->{'dbh'}->prepare($query);
		$dth->execute;
		$dth->finish;

		$query
			= "DELETE from feat_score "
			. "WHERE input_id in ( @ev_ids ) "
			. "AND input_id = $rule_base_id";

		$dth = $self->{'dbh'}->prepare($query);
		$dth->execute;
		$dth->finish;

		$self->{'dbh'}->commit();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		die( "Could not delete annotation rule evidence from database [$@] "
			 . DBI->errstr );
	}
}

#################################### PRIVATE ####################################
sub _do_insert_query {
	my ( $self, $table, $values ) = @_;

	my $first_half  = "INSERT INTO $table (";
	my $second_half = "values ( ";

	foreach my $key ( keys( %{$values} ) ) {
		$first_half  .= " $key,";
		$second_half .= "'" . $values->{$key} . "'" . ", ";
	}

	$first_half  =~ s/,$//;
	$second_half =~ s/,\s$//;

	my $query = $first_half . " ) " . $second_half . " )";

	my $dbh = $self->{'dbh'};

	#	print STDERR $query, "\n";
	my $result;
	eval {
		my $sth = $dbh->prepare($query);
		$sth->execute or die("Could not insert row(s) into $table [$query]");
		$sth = $dbh->prepare('select @@identity');
		$sth->execute;
		my $row = $sth->fetchrow_arrayref();
#		print STDERR Dumper($row);
#		print STDERR $row->[0], "\n";
		$result = $row->[0];

		#	while (my $row = $sth->fetchrow_arrayref()) {
		#		print STDERR "IDENTITY value = $row->[0]\n";
		#	}
		#	$sth->finish();
	};

	if ($@) {
		$self->{'dbh'}->rollback();
		print STDERR "Could not insert data into $table " . DBI->errstr;
	}
	return $result;
}

sub _do_conditional_select {
	my ( $self, $table, $field, $values ) = @_;
	my $first_half  = "select $field from  $table where ";
	my $second_half = "";

	foreach my $key ( keys( %{$values} ) ) {
		$first_half .= $key . "=" . $values->{$key} . " and ";

		#		$second_half .= $values->{$key} . ", ";
	}

	$first_half =~ s/\s+and\s+$//;

	#	$first_half .= ';';
	#	$second_half =~ s/,\s$//;

	#	my $query = $first_half . " ) " . $second_half . " )";
#	print STDERR $first_half, "\n";
	return $self->_do_select_array($first_half);

	#	my $dbh = $self->{'dbh'};
	#	my $result;
	#	eval {
	#		my $sth = $dbh->prepare($first_half);
	#		print STDERR "After prepare\n";
	#		$sth->execute();
	#		print STDERR "After execute\n";
	#		$result =  $sth->fetchall_arrayref();
	#	};
	#	if ($@)  {
##		$self->{'dbh'}->rollback();
	#		croak "Could not insert data into $table ". DBI->errstr;
	#	}
	#	return $result;
}

sub _do_select_query {
	my ( $self, $query, $key ) = @_;

	unless ( defined $self->{'dbh'} ) {
		unless ( $self->connect ) {
			die("Could not connect to the database.  Check parameters and try again."
			);
		}
	}

	my $dbh = $self->{'dbh'};
	my $sth = $dbh->prepare($query);
	$sth->execute;

	return $sth->fetchall_hashref($key);

}

sub _do_select_array {
	my ( $self, $query ) = @_;
#	print STDERR $query;

	unless ( defined $self->{'dbh'} ) {
		unless ( $self->connect ) {
			die("Could not connect to the database.  Check parameters and try again."
			);
		}
	}

	my $dbh = $self->{'dbh'};
	my $sth = $dbh->prepare($query);
#	print STDERR "After prepare\n";
	$sth->execute;
	
    my $results = $sth->fetchall_arrayref();
    $sth->finish;
    
	return $results;
}

sub _init {
	my ( $self, $params ) = @_;

	foreach my $key ( keys %{$params} ) {
		$self->{$key} = $params->{$key};
	}

}

sub DESTROY {
	my $self = shift;

	# check for an overridden destructor...

	if ( defined( $self->{'dbh'} ) ) {
		$self->{'dbh'}->disconnect();
	}

	#if I decide to inherit, don't forget to call this.
	#$self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

1;
