=head1 NAME

Download.pm - Module allowing downloading of external resources

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

=head1 CONTACT

=cut

package Pangenome::Download;
use strict;
use DBI;
use File::Path qw(mkpath remove_tree);
use File::Glob qw(glob);
use File::Basename;
use File::Slurp;
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use Data::Dumper;
use TIGR::Foundation;
use TIGR::FASTAreader;
use TIGR::FASTArecord;
use TIGR::FASTAwriter;

sub new {
    my ( $class, $working_dir, $username,$password ) = @_;
    my $self = {};
    bless $self, ref($class) || $class;

    $self->{WORKING_DIR} = $working_dir;
    $self->{USER}   = $username;
    $self->{PWD}  = $password;
    $self->{PTT_FILE}  = $working_dir . "/fastas/" . "combined.ptt";   

    unlink $self->{PTT_FILE} if (-s $self->{PTT_FILE});

    return $self;
}
sub get_ptt_file{
    my $self = shift;

    if(exists $self->{PTT_FILE}){
	return $self->{PTT_FILE};
    }
}
sub set_db{
    my ($self,$db) = @_;
    $self->{DB} = $db;
}
sub set_asmbl_id{
    my ($self, $id) = @_;
    $self->{ASMBL_ID} = $id;
}
sub set_type{
    my ($self, $type) = @_;
    $self->{TYPE} = $type;
}
sub set_location{
    my ($self, $location) = @_;
    $self->{LOCATION} = $location;
}
sub get_db{
    my $self = shift;
    if (exists $self->{DB}) {
	return $self->{DB};
    }else {
	return;
    }
}
sub get_asmbl_id{
    my $self = shift;
    if (exists $self->{ASMBL_ID}) {
	return $self->{ASMBL_ID};
    }else {
	return;
    }
}
sub get_type{
    my $self = shift;
    if (exists $self->{type}) {
	return $self->{type};
    }else {
	return;
    }
}
sub get_working_dir{
    my $self = shift;
    if (exists $self->{WORKING_DIR}) {
	return $self->{WORKING_DIR};
    }else {
	return;
    }
}
sub get_location{
    my $self = shift;
    if (exists $self->{LOCATION}) {
	return $self->{LOCATION};
    }else {
	return;
    }
}
sub set_fasta_dir{
    my ($self,$dir) = @_;
    $self->{FASTA_DIR} = $dir;
}
sub set_accession_lookup{
    my ($self,$lookup) = @_;
    $self->{ACCESSION_NAME} = $lookup;
}
sub get_accession_name{
    my $self = shift;
    
    if(exists $self->{ACCESSION_NAME}){
	return $self->{ACCESSION_NAME};
    }else{
	return;
    }
}
sub get_sgd_list{
    my $self = shift;

    if(exists $self->{SGD_LIST}){
	return $self->{SGD_LIST};
    }else{
	return;
    }
}
sub fetch_sequences{
    my $self = shift;

    if($self->{TYPE} eq 'NCBI'){
	$self->download_from_ncbi;
    }elsif($self->{TYPE} eq 'JCVI'){
	$self->download_from_jcvi;
    }
}
sub download_from_jcvi{
    my $self = shift;

    my $bin_dir = $FindBin::Bin;
    my $pull_pep_exec   = "$bin_dir/pull_pep_fasta.pl";
    my $fasta_dir = $self->{WORKING_DIR} . "/fastas/";

    my $pull_pep_cmd = "$pull_pep_exec -P $fasta_dir -u $self->{USER} -p $self->{PWD}";

    my $db_param = " -D " . $self->{LOCATION};
    $db_param .= " -a \"$self->{ASMBL_ID}\""  unless($self->{ASMBL_ID} eq 'ISCURRENT');

    system( $pull_pep_cmd . $db_param ) == 0 || die ( "Couldn't pull pep: $? $!");
    
    $self->add_db_to_sgd_list;
}
sub add_db_to_sgd_list{
    my $self = shift;
    push(@{$self->{SGD_LIST}},"$self->{DB}\t$self->{LOCATION}\t$self->{ASMBL_ID}");
}
sub download_from_ncbi{
    my $self = shift;
    
    $self->make_fasta_dir;
    $self->download_ncbi_files;
    $self->store_downloaded_accession_names;

}
sub store_downloaded_accession_names{
    my $self = shift;

    my $dir = $self->{FASTA_DIR};
    my @faa_files = glob("$dir/*.faa");
    my $db_name_hsh;

    foreach my $file(@faa_files){
	my($name,$dir,$suffix) = fileparse($file);
	my($accession, $post) = split(/\.faa/,$name);
	$db_name_hsh->{$accession} = $self->{DB};
    }

    $self->set_accession_lookup($db_name_hsh);
}
sub make_fasta_dir{
    my $self = shift;

    my $dir = $self->{WORKING_DIR} . "/fastas/" . $self->{DB} . "/";
    mkpath($dir) unless (-d $dir);

    $self->set_fasta_dir($dir);
}
sub download_ncbi_files{
    my $self = shift;

    my $location = $self->{LOCATION};
    my $dir = $self->{FASTA_DIR};
    my $base_cmd = "wget -P $dir " . $location;

    my $cmd = $base_cmd  . "\\*faa";
    system($cmd) == 0 || die("Couldn't perform wget command: $!");
    
    $cmd = $base_cmd . "\\*ffn";
    system($cmd) == 0 || die("Couldn't perform wget command: $!");
    
    $cmd = $base_cmd . "\\*ptt";
    system($cmd) == 0 || die("Couldn't perform wget command: $!");

    sleep(1);
}
sub parse_files{
    my $self = shift;

    if($self->{TYPE} eq 'NCBI'){
	$self->parse_ncbi_files;
    }
    
}
sub parse_ncbi_files{
    my $self = shift;
	
    my $accessions = $self->{ACCESSION_NAME};
    my $dir = $self->{FASTA_DIR};
    my($ptt_id_hsh,$ptt_coord_hsh);

    foreach my $acc (keys %$accessions){

	if(-s $dir . "$acc.ptt"){

	    ($ptt_id_hsh,$ptt_coord_hsh) = $self->parse_ncbi_ptt_file($acc);
	    $self->convert_ptt_to_att_file($ptt_id_hsh,$acc,$accessions->{$acc});

	}else{

	    die("$acc.ptt does not exist in $dir");

	}

	if(-s $dir . "$acc.faa"){

	    $self->parse_ncbi_fasta_file($ptt_id_hsh,$ptt_coord_hsh,$acc,"faa");

	}else{

	    die("$acc.faa does not exist in $dir");

	}

	$self->parse_ncbi_fasta_file($ptt_id_hsh,$ptt_coord_hsh,$acc,"ffn") if(-s $dir . "$acc.ffn");
    }
}
sub parse_ncbi_fasta_file{
    my ($self,$ptt_id_hsh,$ptt_coord_hsh,$acc,$type) = @_;
   
    my $tf_object = new TIGR::Foundation;
    my @errors;

    my $in_file = $self->{FASTA_DIR} . $acc  . ".$type";
    my $fr = new TIGR::FASTAreader ($tf_object,\@errors, $in_file);

    my $fasta_file = new TIGR::FASTAwriter();
    my $post_fix = ($type eq 'faa') ? "pep": "seq" ;

    my $foh = $self->{WORKING_DIR} . "/fastas/" . "combined_ncbi_$type." . $post_fix;
    
    $fasta_file->open($foh,"a") || die "Can't open $foh: $!";

    open(my $pseudo,">", $self->{FASTA_DIR} . "$acc.pseudo");

    while($fr->hasNext()){
	my $seq_obj = $fr->next();
	my $header = $seq_obj->{header};
	
	if($type eq 'faa' && $header =~ /gi\|(\d+)/){

	    my $locus = $ptt_id_hsh->{$1}->{locus};
	    my $fasta_record = new TIGR::FASTArecord(">$locus",$seq_obj->{data_rec});
	    $fasta_file->write($fasta_record);

	}elsif($type eq 'ffn'){

	    if($header =~ /:([c]?\d+-\d+)/){

		my($end5,$end3) = split(/-/,$1);
		my $key;
		if($end5 =~ /^c/){
		    $end5 =~ s/^c//;
		    $key = "$end3..$end5";
		}else{
		    $key = "$end5..$end3";
		}

		if(exists $ptt_coord_hsh->{$key}){
		    my $locus = $ptt_coord_hsh->{$key};
		    my $fasta_record = new TIGR::FASTArecord(">$locus",$seq_obj->{data_rec});
		    $fasta_file->write($fasta_record);
		}else{
		    print $pseudo "$header\n";
		}
	    }
	}
    }
}
sub parse_ncbi_ptt_file{
    my ($self,$acc) = @_;

    my $file = $self->{FASTA_DIR} . "$acc.ptt";

    open(my $fh, "<", $file);
    my ($phsh,$chsh);

    while(my $line = <$fh>){
	$line =~ s/\s+$//;
	
	if($line =~ /^\d+\.\.\d+/){
	    my @values = split(/\t/,$line);
	    $phsh->{$values[3]}->{locus} = $values[5];
	    $phsh->{$values[3]}->{coords} = $values[0];
	    $phsh->{$values[3]}->{com_name} = $values[8];
	    $phsh->{$values[3]}->{strand} = $values[1];
	    $chsh->{$values[0]} = $values[5];
	}
    }
    close $fh;
    return ($phsh,$chsh);
}
sub convert_ptt_to_att_file{
    my ($self,$hsh,$acc,$name) = @_;
   
    my $file = $self->{PTT_FILE};
    open(my $fo, ">>",$file);

    foreach my $pid(sort  keys %$hsh){

	print $fo "$acc\t";
	print $fo "$hsh->{$pid}->{locus}\t";
	
	my $strand = $hsh->{$pid}->{strand};
	my $coords = $hsh->{$pid}->{coords};
	my($low,$high) = split(/\.\./,$coords);

	if($strand eq '-'){
	    print $fo "$high\t$low\t";
	}else{
	    print $fo "$low\t$high\t";
	}

	print $fo "$hsh->{$pid}->{com_name}\t";
	print $fo "$name";
	print $fo "\n";
    }

    close $fo;
}
1;
