#!/usr/local/bin/perl

###############################################################################
#                                                                             #
#       Copyright (C) 2016-2017 J. Craig Venter Institute (JCVI).             #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

use warnings;
use strict;
$|++;

=head1 NAME

parse_genbank_files.pl

=head1 SYNOPSIS

    USAGE: ./parse_genbank_files.pl -l <list of genbank files> [-o output/dir]

=head1 OPTIONS

B<--file_list, l>   :   File containing paths to all genbank files to parse

B<--output, o>      :   Output directory [Default: Current working dir]

B<--nuc, n>         :   Work in nucleotide space instead of peptide space.

B<--both>           :   Produce both nucleotide AND peptide files.

B<--length>         :   Features this long and shorter will be excluded from output files.

B<--no_check>       :   Disable the post-processing checks for duplicate loci, etc.

B<--no_dos2unix>    :   Do not run dos2unix on the input set.

B<--verbose, -v>    :   Print much, much more detailed output regarding skipped features, etc.

B<--help, -h>       :   Displays Help

=head1 DESCRIPTION

This script takes in a list of genbank files. It then parses each file and creates a protein fasta file
and a gene attribute file which can then be used as input to the Pangenome pipeline.  If --nuc is used, a
nucleotide fasta file is created instead, and the gene attribute file and fasta file will include more than
just genes; it will also contain various rna features and pseudogenes, etc., longer than --length (default: 30)

Once the set of files has been downloaded, runs B<'check_att_files.pl'> to identify and resolve any cases of:

 1. duplicate locus tag prefixes used between .gb files
 2. duplicate locus tags used within a give .gb file.

See the --help info for B<'check_att_files.pl'> for more information on the particulars of the resolution process.

=head1 OUTPUTS

 .nuc - A nucleotide fasta file with a sequence for each nuc feature.
 .pep - A protein fasta file with a sequence for each gene
 .natt - A nucleotide gene attribute file
 .patt - A protein gene attribute file
 .pseudo - A list of all found pseudo genes (only applies to protein-space)

genbank_printing_overview.txt - A file that lists the number of genes found for each genome and prints "NO ANNOTATION" if no annotation was found on that genome

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;
use File::Basename;
use FindBin;
#use lib ("/usr/local/devel/ANNOTATION/jinman/lib/perl5/");
use Bio::SeqIO;
use Bio::Seq;
use Cwd;
use Pod::Usage;

my $DEFAULT_LENGTH = 30;
my $DEFAULT_VERBOSE = 0;
my $BIN_DIR = $FindBin::Bin;
my $CHECK_ATT_FILES_EXEC    = "$BIN_DIR/check_att_files.pl";
my $CHECKER_INPUT_FILE      = "check_att_files.input";
my $DOS2UNIX_EXEC           = "/usr/bin/dos2unix";
my $NUC_ATT_SUFFIX          = 'natt';
my $PROT_ATT_SUFFIX         = 'patt';
my $cutoff_len;
#Grab passed in arguments/options
my %opts;
my $verbose;
GetOptions (\%opts, 
            'file_list|l=s',
            'output|o=s',
            'nuc|n',
            'both',
            'length=i',
            'no_check',
            'no_dos2unix',
            'verbose|v',
            'help|h',
            ) || die "Can't get options!\n$!\n";

pod2usage( { -verbose => 2, -exitval => 0 } ) if $opts{help};

my $OUTPUT = check_params();
my ($assembly,$features);

#Open List of GenBank file locations
open( my $lfh, "<", $opts{file_list}) || die "Can't open $opts{file_list}: $!\n";

# create genomes.list file if we will be running the locus checks:
my %seen_genomes;

my $sum_features;
# Go through each GenBank file that was passed in
while( my $file = <$lfh> ) {
    
    chomp( $file );

    #Remove existing files for this accssion
    my ( $filename, $path, $suffix ) = fileparse( $file );
    $filename =~ s/\..*//;
    remove_existing_files( $filename );

    unless ( $opts{ no_dos2unix } ) {
        my @cmd = ( $DOS2UNIX_EXEC, $file );
        system( @cmd ) && die "Problem running dos2unix on $file!\n";
    }

    if ( $opts{nuc} || $opts{both} ) { 
        parse_nuc_features( $file, \%seen_genomes ); 
    }
    if ( $opts{ both } || ! $opts{ nuc } ) {

        my $feature_counts; #keeps track of genes/pseudo genes for a set of files

        #Open file and parse file name, used as output name
        open( my $fh, "<", $file ) || die "Can't open $file: $!\n";
        #Go through each line of the current GenBank file and parse 
        #annotation and sequences
        while( my $line = <$fh>){
        
            #Parse out the keyword
            my $keyword = "";
            if ( $line =~ /^\s*([A-Z]+)/ ) {
                $keyword = $1;
            } else {
                next;
            }
	    
            #Store Keyword's value in assembly hash
            if ( $keyword eq "ACCESSION" ) {
              
                if ( $line =~ /^\s*ACCESSION\s*(\w+)/ ){
                    $assembly->{'accession'} = $1 unless($1 eq "unknown");
                }

            } elsif ( $keyword eq "VERSION" ) {

                if ( $line =~ /^\s*VERSION\s*(\S+)\s*(\S+)/ ) {
                    ( $assembly->{'version'}, $assembly->{'gi'} ) = ( $1,$2 )
                }

            } elsif ( $keyword eq "ORGANISM" ) {

                $assembly->{'organism'} = $1 if ( $line =~ /\s*ORGANISM\s*(.*)/ );
                $assembly->{'organism'} =~ s/\s+$//;

            } elsif ( $keyword eq "LOCUS" ) {

                if ( $line =~ /^LOCUS\s+(\w+)/ ) {
                    $assembly->{'accession'} = $1 
                }

            } elsif ( $keyword eq "FEATURES" ) {

                #We parse the features here. This is the hard part.
                $assembly->{'features'} =  &parse_features( $fh, $assembly );
		
		$seen_genomes{ $filename }++; # Add to the hash of seen genomes if we're doing checks.
	
                $feature_counts = &print_files($assembly->{'accession'},$assembly->{'features'},$filename,$feature_counts);


            }
        }

        close $fh;

        #Print overview file of gene/counts
        open( my $o_fh, ">>", $OUTPUT . "/genbank_printing_overview" . ".txt" );
        print $o_fh "$filename";

        if ( exists $feature_counts->{$filename} ) {

            print $o_fh "\tGenes:$feature_counts->{$filename}->{genes}" if(exists $feature_counts->{$filename}->{genes});
            print $o_fh "\tPseudo:$feature_counts->{$filename}->{pseudo}" if(exists $feature_counts->{$filename}->{pseudo});
            print $o_fh "\n";

        } else {

            print $o_fh "\tNO ANNOTATION\n";

        }

        close $o_fh;

    }

}


unless ( $opts{ no_check } ) {
    open( my $gfh, '>', "$OUTPUT/$CHECKER_INPUT_FILE" ) || die "Can't open $OUTPUT/$CHECKER_INPUT_FILE for writing: $!\n";
    print $gfh map{ "$_\n" } sort keys %seen_genomes;
}


unless ( $opts{ no_check } ) {

    my $cmd = "$CHECK_ATT_FILES_EXEC -g $OUTPUT/$CHECKER_INPUT_FILE -a $OUTPUT -f $OUTPUT --resolve";

    if ( $opts{ both } || ! $opts{ nuc } ) {
        my $check_att_files_log = "$OUTPUT/check_att_files.prot.log";
        my $prot_cmd = "$cmd > $check_att_files_log";
        system( $prot_cmd ) && die "There was a problem running check_att_files with the invocation:\n$cmd\n";
    }
    
    if ( $opts{ nuc } || $opts{ both } ) {
        my $check_att_files_log = "$OUTPUT/check_att_files.nuc.log";
        $cmd .= " --nuc > $check_att_files_log";
        system( $cmd ) && die "There was a problem running check_att_files with the invocation:\n$cmd\n";
    }

}

exit(0);

# Subs! 

sub parse_nuc_features {

    my ( $gb_file, $seen_genomes_ref ) = @_;

    #Open file and parse file name, used as output name
    open(my $fh, "<", $gb_file) || die "Can't open $gb_file: $!\n";
    my ($filename,$path,$suffix) = fileparse($gb_file);
    $filename =~ s/\..*//;
 
    # create a SeqIO object for the genbank file to make use of the parsing available.
    my $seqio_object = Bio::SeqIO->new( -file => $gb_file, -format => 'genbank' );
    my $is_first = 1;
    my $seq_object;
    
    my $features = {};
    my $feature_counts;
    my %genes;
    my %CDS;

    while ( $seq_object = $seqio_object->next_seq() ) {

        # only need to get these once, assume these are the same for all records in this file.
        if ( $is_first ) {

            # get sequence version
            ## Don't think this is needed anymore
            $assembly->{version} = $seq_object->seq_version;

            # Organism is a little more convoluted to get:
            $assembly->{organism} =  $seq_object->species->node_name;

            $is_first = 0;

        }

        # get accession
        my $accession = $seq_object->accession_number();

        for my $feat_object ( $seq_object->get_SeqFeatures ) {

            my ( $acc, $locus_tag, $product, $key, $end5, $end3, $sequence, $pseudo, $translation );

            # Skip these features in the feature table:
            if ( $feat_object->primary_tag =~ /(?:source|assembly_gap|STS|repeat_region|misc_binding|misc_feature)/ ) {
                warn "Skipping " . $feat_object->primary_tag . " feature at " . $feat_object->start .
                         '..' . $feat_object->end . "\n" if $verbose;
                $feature_counts->{ $filename }->{ skipped_type }++;
                next;

            }

            # get coords
            ( $end5, $end3 ) = ( $feat_object->start, $feat_object->end );
            # length check
            if ( $end3 - $end5 + 1 < $cutoff_len ) {
                unless ( $end3 == 1 ) {
                    warn "Skipping feature at $end5..$end3 because it is smaller than the cutoff of $cutoff_len.\n" if $verbose;
                    $feature_counts->{ $filename }->{ skipped_length }++;
                    next;
                }
            }

            $feature_counts->{ $filename }->{ $feat_object->primary_tag }++;

            $key = "$end5..$end3";
            ( $end5, $end3 ) = ( $end3, $end5 ) if ( $feat_object->strand == -1);

            # mark pseudo if present
            if ( $feat_object->has_tag( 'pseudo' ) ) {

                $pseudo++;
                $feature_counts->{ $filename }->{ pseudo }++;

            }

            # get feat locus
            if ( $feat_object->has_tag( 'locus_tag' ) ) {

                ( $locus_tag ) = $feat_object->get_tag_values( 'locus_tag' );

            } else {

                warn "Found a " . $feat_object->primary_tag . " feature with no locus_tag!\n" if $verbose;
                next;

            }

            # get product name
            if ( $feat_object->has_tag( 'product' ) ) {

                ( $product ) = $feat_object->get_tag_values( 'product' );

           # } elsif ( $pseudo && $feat_object->has_tag( 'notes' ) ) {
            } elsif ( $pseudo ) {

                if ( $feat_object->has_tag( 'note' ) ) {
                    my ( $note ) = $feat_object->get_tag_values( 'note' );
                    if ( $note =~ /([^;]+);/ ) { # Capture everything up to the first semi-colon in /notes for pseudo-genes.
                        $product = $1;
                    }
                } else {
                    warn 'Feature ' . $feat_object->primary_tag . " is /pseudo but has no /note field to parse\n" if $verbose;
                }

            } else {

                warn 'Skipping ' . $feat_object->primary_tag . " feature $locus_tag with no product!\n" if $verbose;
                next;

            }

            # get sequence
            ( $sequence ) = $feat_object->seq->seq();

            if ( $sequence =~ /^\s+$/ ) {

                warn "Skipping $locus_tag, can't get sequence from $end5..$end3\n" if $verbose;
                next;

            }

            # and the translation
            if ( $feat_object->has_tag( 'translation' ) ) {

                ( $translation ) = $feat_object->get_tag_values( 'translation' );

            }

            # Now put all the pieces in place
            $features->{ $key }->{acc}  = $accession;
            $features->{ $key }->{end5} = $end5;
            $features->{ $key }->{end3} = $end3;
            $features->{ $key }->{locus_tag} = $locus_tag;
            $features->{ $key }->{pseudo}++ if $pseudo;
            $features->{ $key }->{product} = $product;
            $features->{ $key }->{sequence}  = $sequence;
            $features->{ $key }->{translation} = $translation;

        }

    }

    # Now we have all the features.  Print them:
    print_nuc_files( $features, $filename );

    open( my $o_fh, '>>', $OUTPUT . '/genbank_printing_overview.txt' );
    print $o_fh "$filename";

    if ( exists $feature_counts->{ $filename } ) {

        for my $feat_type ( keys %{$feature_counts->{ $filename }} ) {

            print $o_fh "\t$feat_type:$feature_counts->{ $filename }->{ $feat_type }";

        }
        print $o_fh "\n";

        # Add it to the list to check
        $seen_genomes_ref->{ $filename }++;

    } else {

        print $o_fh "\tNO_ANNOTATION\n";

    }
    
}


sub print_nuc_files {

   my ( $features, $filename ) = @_;

    open( my $att_fh, '>>', $OUTPUT . "/$filename.$NUC_ATT_SUFFIX" ) || die "Can't open att_file: $!\n";
    my $nuc_obj = Bio::SeqIO->new(-file => '>>' . $OUTPUT . "/$filename" . '.nuc', -format => 'fasta' );

    foreach my $key ( sort { $features->{ $a }->{acc}.$a cmp $features->{ $b }->{acc}.$b } keys %$features ) {

        # print Fasta
        my ( $acc, $locus, $end5, $end3 ) = @{$features->{ $key }}{ qw(acc locus_tag end5 end3) };
        my $seq_obj = Bio::Seq->new( -seq => $features->{ $key }->{sequence},
                                     -display_id => $locus,
                                     -alphabet => 'dna' );
        $nuc_obj->write_seq( $seq_obj );

        # print att_file line
        my $product;
        if (!($features->{$key}->{product})) {

            if ($features->{$key}->{note}) {
                $product = $features->{$key}->{note};
            } elsif ($features->{$key}->{comment}) {
                $product = $features->{$key}->{comment};
            } elsif ($features->{$key}->{gene}) {
                $product = $features->{$key}->{gene};
            }

        } else {
            $product = $features->{$key}->{product};
        }

        unless ( $product ) {
            warn( "WARNING: Missing product name for $key from $acc\n") if $verbose;
            $product = "[NO PRODUCT NAME]";
        }
        unless ( $acc ) { die "NO ACC\n" }
        unless ( $locus ) { die "NO LOCUS, FOR $end5..$end3\n" }
        unless ( $end5 ) { die "NO END5\n" }
        unless ( $end3 ) { die "NO END3\n" }
        unless ( $filename ) { die "NO FILENAME\n" }

        print $att_fh "$acc\t$locus\t$end5\t$end3\t$product\t$filename\n";

    }

}


sub print_files{

    my ($acc,$features,$filename,$feature_counts) = @_;

    open(my $att_fh, ">>", $OUTPUT . "/$filename.$PROT_ATT_SUFFIX");
    open(my $pseudo_fh, ">>", $OUTPUT . "/$filename" . ".pseudo");

    my $pep_obj = Bio::SeqIO->new(-file => '>>' . $OUTPUT . "/$filename" . '.pep', -format => 'fasta' );

    foreach my $key (keys %$features) {

        unless (exists $features->{$key}->{not_gene}) {

            if (exists $features->{$key}->{translation}) {

                $feature_counts->{$filename}->{genes}++;

                $features->{$key}->{translation} =~ s/(\s+|\n+)//g;
                my $locus = $features->{$key}->{'locus_tag'};

                #Print Fasta For Locus
                my $seq_obj = Bio::Seq->new(-seq => $features->{$key}->{translation},                        
                                -display_id => $locus,                        
                                -alphabet => "dna" );
                
                $pep_obj->write_seq($seq_obj);

                #Print Gene Att
                my $product;

                if (!($features->{$key}->{product})) {

                    if ($features->{$key}->{note}) {
                        $product = $features->{$key}->{note};
                    } elsif ($features->{$key}->{comment}) {
                        $product = $features->{$key}->{comment};
                    } elsif ($features->{$key}->{gene}) {
                        $product = $features->{$key}->{gene};
                    }

                } else {
                    $product = $features->{$key}->{product};
                }

                if ($product) {
                    print $att_fh "$acc\t$locus\t$features->{$key}->{end5}\t$features->{$key}->{end3}\t$product\t$filename\n";
                } else {
                    die ("ERROR: Missing product name for $key from $acc\n");
                }

            } else {

                if (exists $features->{$key}->{pseudo}) {

                    if ($key) {

                        $feature_counts->{$filename}->{pseudo}++;
                        print $pseudo_fh "$key";
                        print $pseudo_fh "\t$features->{$key}->{'locus_tag'}" if exists ($features->{$key}->{'locus_tag'});
                        print $pseudo_fh "\n";

                    }

                }

            }

        }

    }


    return $feature_counts;
	
}


sub parse_features {

    my ($fh, $assembly) = @_;

    my $pre_feature = {};

    #Cycle through features section and store information.  Will filter it later
    my $key;  # Will hold the` coords string (i.e. 246..2043)
        
FEATURES: while( my $line = <$fh> ) {

        chomp( $line );

        #If we find a cds or gene feature
        if( $line =~ /^     (CDS|gene)\s+(.*)/ ) {
           
            #Only storing CDS information
            if ($1 eq 'gene') {
                $key = "";
            } else {

                my $coords = $2;
                my ($e5, $e3);
                
                #Check to see if we might have a multi-line value
                #if there is a trailing comma, the coords span lines
                #either due to a complement or join (or both)
                
                while( $coords =~ /,$/ ) {

                    my $tmp = <$fh>;
                    $tmp =~ s/^\s+(.*)\n/$1/;
                    $coords .= $tmp;

                }

                #Parse end5, end3 from coordinate line
                ($key, $e5, $e3) = &parse_coord_string( $coords );

		#Store the information
                $pre_feature->{$key}->{'end5'} = $e5;
                $pre_feature->{$key}->{'end3'} = $e3;
                
            }

        } elsif( $line =~ /(\d+\.\.\d+)/ ) {

            #This is not a gene or cds and we don't want it.
            #This is like an RNA or something
	    #Allo for genes that take up the whole assembly
	    unless($line =~ /(^     source|CONTIG)/ ){
		$pre_feature->{$1}->{'not_gene'} = $1;
		$key = "";
	    }

        } elsif($line =~ /\/pseudo/){
            $pre_feature->{$key}->{pseudo} = 1 ;
        }
 
        if ( $line =~ /^\s*ORIGIN/ ) {
            last FEATURES;
        }

	#From here and below, we only look into these if we have an active key
        #If our key is eq "", this means we don't want to parse the stuff
        next if( !defined( $key ) || $key eq "" );

        if ($line =~ /^\s+\/(\w+)\=(.*)/) {
            
            if ( $line =~ /^\s+\/(\w+)\=(.*)/ ) {

                my ($att, $val) = ($1,$2);
                &process_gene_line( $att, $val, $fh, $pre_feature->{$key} );

            } 

        } 
	
    }

    
    return $pre_feature;

}


sub process_gene_line {

    my ($att, $val, $fh, $pre_feature, $debug) = @_;

    #Check to see if we might have a multi-line value
    #if there are quotes and it's not closed, the we have a multiline value
    while( $val =~ /^\s*\"/ && $val !~ /\"\s*$/ ) {

        my $tmp = <$fh>;
        $tmp =~ s/(^\s+|\s+\Z)//g;
        $val .= " ".$tmp;

    }

    #We can remove the quotes now, as they have served their purpose
    $val =~ s/\"//g;
 
    $pre_feature->{$att} = $val;

} #end [sub process_gene_line {]


sub parse_coord_string {

    my ($coords) = @_;
    my ($key, $e5, $e3);

    #Get rid of partial gene indicators
    $coords =~ s/[\<\>]//g;

    #Get the key and end5 and end3
    if ( $coords =~ /join\((\d+)\.\.\d+,1\)/ ) {
        $e5=$1;
        $e3=1;
        $key="$e5..1";
    } elsif ( $coords =~ /join\((.*)\)/ ) {
    #if ( $coords =~ /join\((.*)\)/ ) {

        my $values = $1;
        $values =~ s/[\(\)]//;

        my @coords = split(/\,/,$values);
        my ($te5,$te3) = (0,0);

        foreach(@coords){

            my ($end5,$end3) = split(/\.\./,$_);
            
            $te5 = $end5 if($te5 == 0);
            $te5 = $end5 if($end5 < $te5);
            $te3 = $end3 if($end3 == 0);
            $te3 = $end3 if($end3 > $te3);

        }

        $e5 = $te5;
        $e3 = $te3;
        $key = "$e5.$e3";

    } elsif ( $coords =~ /(\d+).*\.\.(\d+)/ ) {

        $key = "$1..$2";
        ($e5,$e3) = ($1,$2);

    }

    unless (defined $key && defined $e5 && defined $e3 ) {
        die("Could not parse gene/CDS coords: $coords");
    }

    #Complement (if needed)
    if ( $coords =~ /complement/ ) {

        #Reverse the coords
        ($e5,$e3) = ($e3,$e5);    

    }

    return ($key, $e5, $e3);

}


sub check_params{

    my $output = $opts{output} // cwd;
    my $errors;
    
    mkdir($output) unless (-s $output);
    
    if ( $opts{file_list} ) {
        $errors .= "File does not exist or is size zero: $opts{file_list}\n" unless ( -s $opts{file_list});
    } else {
        $errors .= "Must provide parameter: --file_list\n" unless ( $opts{file_list} );
    }

    $cutoff_len = $opts{length} // $DEFAULT_LENGTH;
    $verbose    = $opts{verbose} // $DEFAULT_VERBOSE;

    unless ( $opts{ no_check } ) {

        unless ( -s "$BIN_DIR/check_att_files.pl" ) {
            $errors .= "Can't find ./check_att_file.pl, which is needed without --no_check\n";
        }

    }

    $errors .= "\nUsage: ./parse_genbank_files.pl -l genbank_file_list.txt [-o output/dir]\n" if $errors;

    die( $errors ) if $errors;

    return( $output );

}


sub remove_existing_files {

    my $filename = shift;

    unlink("$OUTPUT/$filename.nuc");
    unlink("$OUTPUT/$filename.natt");
    unlink("$OUTPUT/$filename.pep");
    unlink("$OUTPUT/$filename.patt");
    unlink("$OUTPUT/$filename.pseudo");

}
