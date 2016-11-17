=head1 NAME

NCBIDownload.pm - Module allowing downloading of NCBI files through eutils

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

=head1 CONTACT

=cut

package Pangenome::NCBIDownload;
use strict;
use DBI;
use File::Path qw(mkpath remove_tree);
use File::Glob qw(glob);
use File::Basename;
use File::Slurp;
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use Data::Dumper;
use LWP::Simple;
use Cwd;

my $MAX_RETURN = 500;
my $MAX_RETRIES = 3;

sub new {

    my ( $class,$output ) = @_;
    my $self = {};
    bless $self, ref($class) || $class;

    $self->{OUTPUT_DIR} = $output // cwd();
    $self->{EUTIL} = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

    return $self;
}

sub get_eutil{
    
    my $self = shift;

    return $self->{EUTIL} if(exists $self->{EUTIL});

}
sub get_output_dir{

    my $self = shift;

    return $self->{ACC} if(exists $self->{ACC});
}

sub set_accession{

    my ($self,$acc) = @_;

    $self->{ACC} = _trim($acc);
    $self->{ACC_TYPE} = $self->_determine_type;

}

sub get_accession{

    my $self = shift;

    return $self->{ACC} if(exists $self->{ACC});

}

sub get_accession_type{

    my $self = shift;

    return $self->{ACC_TYPE} if(exists $self->{ACC_TYPE});

}

sub set_accession_name{

    my ($self,$name) = @_;

    $name =~ s/\s+/_/g;
    $name =~ s/\//_/g;

    $self->{ACC_NAME} = _trim($name);

}

sub get_accession_name{

    my $self = shift;

    return $self->{ACC_NAME} if(exists $self->{ACC_NAME});

}

sub set_gb_file{

    my ($self,$file) = @_;

    $self->{GB_FILE} = $file;

}

sub get_gb_file{

    my $self = shift;

    return $self->{GB_FILE} if(exists $self->{GB_FILE});

}

sub set_fasta_file{

    my ($self,$file) = @_;

    $self->{FASTA_FILE} = $file;

}

sub get_fasta_file{

    my $self = shift;

    return $self->{FASTA_FILE} if(exists $self->{FASTA_FILE});

}

sub set_uids{

    my ($self,$uids) = @_;

    $self->{UIDS} = $uids;

}

sub get_uids{

    my $self = shift;

    return $self->{UIDS} if(exists $self->{UIDS});

}

sub find_uids{

    my $self = shift;
    my @ids;

    if ( $self->get_accession_type eq 'mixed' ) {

        # split accessions into two groups: default and other
        # run the esearch eutil twice, getting the appropriate ids for each and storing them all in @ids
        my ( @default_accs, @refseq_accs, @other_accs );
        for my $acc ( split( ',', $self->get_accession ) ) {

            if (    ( $acc =~ /^SAMN/ ) || 
                    ( $acc =~ /^[A-Z]{4}\d{2,8}/ ) ) {

                push @other_accs, $acc;

            } elsif ( $acc =~ /^[A-Z]{2}_[A-Z]{4}\d{2,8}/) {

                push @refseq_accs, $acc

            } else {

                push @default_accs, $acc;

            }

        }

        # build urls:
        my ( $defaults_search_url, $other_search_url, $refseq_search_url );
        my $search_url = $self->{EUTIL} . "esearch.fcgi?db=nucleotide&tool=gbdownload&email=jinman\@jcvi.org&retmax=5000&term=";

        if ( scalar( @default_accs ) ) {
            $defaults_search_url = $search_url . join(',', @default_accs);
        }
        if ( scalar( @other_accs ) ) {
            $other_search_url    = $search_url . join(',', @other_accs)  . '+AND+wgs_contig[prop]';
        }
        if ( scalar( @refseq_accs ) ) {
            $refseq_search_url   = $search_url . join(',', @refseq_accs) . '+AND+"meta%20Genome-Annotation-Data"[prop]';
        }

        foreach my $url ( $defaults_search_url, $refseq_search_url, $other_search_url ) {

            next unless $url;

            my $search_result = undef;
            my $retry_counter = 0;
            do {

                $search_result = get( $url );
                $retry_counter++;
                sleep 1;

            } until ( $search_result || $retry_counter == $MAX_RETRIES );

            if ( ! $search_result ) {

                warn "Esearch URL returned no data: $search_url\n";

            } else {

                my $result_count = 0;
                if ( $search_result =~ /<Count>(\d+)<\/Count>/s ) {
                    $result_count = $1;
                } else {
                    warn "Couldn't get uid count for ". $self->get_accession_name(). "!\nURL: $search_url\n";
                }

                if ( $result_count == 0 ) {
                    warn( "No uid returned for ". $self->get_accession_name() . "!\nURL: $search_url\n" );
                } elsif ( $result_count > $MAX_RETURN ) {
                    warn( "Number of uids found ($result_count) exceeds maxinum of $MAX_RETURN for ".$self->get_accession_name() . "! URL:\n$search_url\n");
                } else {

                    my @lines = split( /\n/, $search_result );

                    #Get IDs associated with accession
                    foreach my $line (@lines){

                        $line =~ s/\s+$//;

                        #Pull NCBI files for each id
                        if ( $line =~ /<Id>(\d+)<\/Id>/ ) {

                            push(@ids,$1);

                        }

                    }

                }
    
            }

        }
        

    } else { # ( if type is not 'mixed' )

        #Make search url and get results
        my $search_url = $self->{EUTIL} . "esearch.fcgi?db=nucleotide&tool=gbdownload&email=jinman\@jcvi.org&retmax=5000&term=" . $self->get_accession;
        #$search_url .= '+AND+"meta%20Genome-Annotation-Data"[prop]' if ( $self->get_accession_type eq 'refseq' );
	$search_url .= '+AND+srcdb_refseq[PROP]' if ($self->get_accession_type eq 'refseq');
        $search_url .= '+AND+wgs_contig[prop]' if ( ( $self->get_accession_type eq 'wgs' ) || 
                                                    ( $self->get_accession_type eq 'biosample' ) );
        my $search_result = get($search_url);

        if ( ! $search_result ) {

            warn "Esearch URL returned no data: $search_url\n";

        } else {

            my $result_count = 0;
            if ( $search_result =~ /<Count>(\d+)<\/Count>/s ) {
                $result_count = $1;
            } else {
                warn 'Could not get count of results for ' . $self->get_accession_name() . ". URL:\n$search_url\n";
            }

            if ( $result_count == 0 ) {
                warn( "No uids returned for ". $self->get_accession_name() . "!\nURL: $search_url\n" );
            } elsif ( $result_count > $MAX_RETURN ) {
                warn( "Number of uids found ($result_count) exceeds limit of $MAX_RETURN for " . $self->get_accession_name() . ". URL:\n$search_url\n");
            } else {

                my @lines = split(/\n/,$search_result);

                #Get IDs associated with accession
                foreach my $line (@lines){
                    
                    $line =~ s/\s+$//;
                    
                    #Pull NCBI files for each id
                    if($line =~ /<Id>(\d+)<\/Id>/){

                        push(@ids,$1);

                    }

                }

            }

        }

    }

    $self->set_uids( scalar( @ids ) ? join( ",", @ids ) : ''  );

}


sub download_files {

    my ($self,$download_type) = @_;

    $self->find_uids;

    my $lc_type = lc($download_type);

    $self->download_fasta() if ($download_type =~ /(fasta|all)/);
    $self->download_gb() if ($download_type =~ /(gb|all)/);

    #Must do sleep to ensure we dont exceed NCBI's limit an accessing
    sleep(1);
    
}

sub download_gb {
    
    my ($self,$history) = @_;

    my $outfile = $self->{OUTPUT_DIR} . "/" . $self->{ACC_NAME} . ".gb";

    if ( $self->get_uids ) {

        my $fetch_gb_url = $self->{EUTIL} . "efetch.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text&tool=gbdownload&email=jinman\@jcvi.org&id=" . $self->get_uids;

        getstore($fetch_gb_url,$outfile);

        $self->set_gb_file($outfile);

    }

}

sub download_fasta {

    my ($self,$history) = @_;

    my $outfile = $self->{OUTPUT_DIR} . "/" . $self->{ACC_NAME}  . ".fasta";

    if ( $self->get_uids ) {

        my $fetch_fasta_url = $self->{EUTIL} . "efetch.fcgi?db=nucleotide&rettype=fasta&tool=gbdownload&email=jinman\@jcvi.org&retmode=text&id="  . $self->get_uids;

        getstore($fetch_fasta_url,$outfile);
        
        #File clean up
        _clean_fasta($outfile);
        unlink($outfile . "_orig");
        unlink("cleanFasta.discard");
        $self->set_fasta_file($outfile);

    }

}

sub _clean_fasta{
 
    my $file = shift;
    my $exe = "/usr/local/common/cleanFasta";

    my $cmd = $exe . " $file";

    system($cmd) == 0 || warn("ERROR: Could not run $cmd");

}

sub _trim{
    
    my $value = shift;

    $value =~ s/\s+$//;
    $value =~ s/^\s+//;

    return $value;

}

sub _determine_type{

    my $self = shift;

    my $acc = $self->get_accession;

    my $acc_type = '';

    my $biosample_found = 0;
    my $wgs_found       = 0;
    my $refseq_found   = 0;
    my $default_found   = 0;

    #Determines what type of accession WGS or BioSample
    for my $acc ( split( ',', $self->get_accession ) ) {

        if ( $acc =~ /^SAMN/ ) {

            $biosample_found = 1;

        } elsif ( $acc =~ /^[A-Z]{4}\d{2,8}/ ) {

            $wgs_found = 1;

        } elsif ( $acc =~ /^[A-Z]{2}_[A-Z]{4}\d{2,8}/ ) {

            $refseq_found = 1;

        } else {

            $default_found = 1;

        }

    }

    if ( ( $biosample_found + $wgs_found + $refseq_found + $default_found ) > 1 ) {

        $acc_type = 'mixed';

    } elsif ( $biosample_found ) {

        $acc_type = 'biosample';

    } elsif ( $wgs_found ) {

        $acc_type = 'wgs';

    } elsif ( $refseq_found ) {

        $acc_type = 'refseq'

    } else {

        $acc_type = 'default';

    }

    return $acc_type;

}

1;
