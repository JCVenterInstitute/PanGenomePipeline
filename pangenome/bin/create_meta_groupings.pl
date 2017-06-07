#!/usr/local/bin/perl -w

use warnings;
use strict;
$|++;

=head1 NAME

create_meta_groupings.pl

=head1 SYNOPSIS

    USAGE: ./create_meta_groupings.pl -a gene_attribute.dat -c centroids.fasta -m panoct.result -f metadata_groupings.txt -d db.list [-r role_lookup.txt -t threshold_file -o output/dir --print_amb --help  'data_file|a=s']
	  

=head1 OPTIONS

B<--data_file, a>        :   The gene_attribute .dat file found in the pangenome project directory

B<--centroids, c>        :   centroids.fasta file found in pangenome results directory

B<--method_result, -m>   :   The panoct.result file found in the pangenome results directory

B<--metadata_file, f>    :   File specifiying the genome name and the associated label

B<--db_list, d>          :   db.list file used in the pangenome run

B<--threshold_file>      :   [Optional] File giving the thresholds to use for label cutoffs

B<--role_lookup, r>      :   [Optional] role_id_lookup.txt if it exists it should be in the pangenome results directory. It only generates if there are SGD genomes.

B<--print_amb>           :   [Optional] Prints the clusters that did not meet the true and false thresholds

B<--index>               :   [Optional] Include fGIs index file to produce grouping results as related to fGIs

B<--help>                
	    
=head1 DESCRIPTION

This scripts takes a metadata grouping file and based on each genome's group determines if a cluster's genome membership meets the true and false thresholds. Based 
on the thresholds the clusters are given a label. 

=head1 INPUTS

--metadata_file : This file is a tab delimited file to assing a genome a group in the format of <genome><tab><group><tab><label>. The genome name must match what's used in the db.list file in the pangenome run.

 %cat groupings.txt
 M0001	plaque_size	small
 CM0002	plaque_size	small
 CM0007	plaque_size	small
 CM0008	plaque_size	small
 CM0001var	plaque_size	large
 F	plaque_size	large
 McKrae	plaque_size	large
 KOS	plaque_size	large
 17	plaque_size	large

--threshold_file: A tab delimited file to specify what thresholds should be generated in the format of <true threshold><tab><false threshold>
For a cluster to be given a label it must contain at least the high threshold of genomes and no more than the false threhold of genomes. 

 %cat threshold.txt
 100  0
 75   25
 50   0

--fGIs index file: File produced by the pangenome pipeline. Tab delimited file containing the cluster number and it's associated fGIs. 

 1       87
 2       87
 3       87
 4       87
 5       87
 6       87
 7       87
 8       87
 9       87
 10      87

=head1 OUTPUTS

=over 1

<name>.clusters_<high threshold>_<low threshold>.txt - A tab delimited file for each threshold combination listing the clusters, their labels and the associated annotation

<name>.counts_<high threshold>_<low threshold>.txt - An overview file that gives an overview of the genomes and their associated label. As well as what number of genomes
must be present for a given label for that cluster to be said to have that label based on the thresholds.

role_id_counts_<high threshold>_<low threshold>.txt - Tab delimited file that displays how many clusters in a specifed label where assigned the specific role_ids. This only gets generated if role_ids were given as a paramter.

=back

=head1 CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Cwd;
use Data::Dumper;
use File::Slurp;
use MLDBM qw( DB_File Data::Dumper );
use Fcntl;
use File::Copy;
use File::Basename;
use File::Path qw(make_path remove_tree);

my %opts;

GetOptions(
    \%opts,              'data_file|a=s',    'centroids|c=s', 'method_result|m=s',
    'metadata_file|f=s', 'threshold_file=s', 'output|o=s',    'role_lookup|r=s',
    'db_list|d=s',       'print_amb',        'index=s',       'help|h'
) || die("Error getting options! $!");

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my $default_cutoffs = set_default_cutoffs();
my ( $THRESHOLDS, $OUTPUT, $threshold_cutoffs ) = check_params($default_cutoffs);

#opens att file tied hash
tie my %ANNOT_HSH, 'MLDBM', $opts{data_file} or die "Can't open MLDBM file: $!\n";

my ( $DATABASES, $db_array ) = process_db_file( $opts{db_list} );
my @dbs = @$db_array;

my ( $C_FGIS, $F_FGIS ) = parse_index_file( $opts{index} ) if $opts{index};

my $CENTROIDS = parse_centroids( $opts{centroids} );
my ( $META_DATA, $META_COUNTS, $META_LABELS, $META_OPTIONS ) = parse_metadata( $opts{metadata_file}, $threshold_cutoffs );
my ( $all_clusters, $meta, $cluster_role_id ) = parse_cluster_results( $opts{method_result} );
my ( $ROLE_ID_LOOKUP, $SUBROLE_LOOKUP ) = parse_role_id_lookup( $opts{role_lookup} ) if $opts{role_lookup};
print_meta_groups( $meta, $all_clusters, $threshold_cutoffs, $cluster_role_id, $C_FGIS, $F_FGIS );

exit(0);


sub parse_index_file {

    my $file = shift;
    my $cluster_hsh;
    my $fgi_hsh;

    open( my $fh, "<", $file );

    foreach (<$fh>) {

        my $line = $_;
        $line =~ s/\s+$//;

        my ( $cluster, $fig ) = split( /\t/, $line );
        $cluster_hsh->{$cluster} = $fig;
        $fgi_hsh->{$fig}->{$cluster} = 1;

    }

    return ( $cluster_hsh, $fgi_hsh );

}


sub parse_centroids {

    my $file = shift;

    my @def_lines = `grep \">\" $file`;

    my $hsh;
    foreach my $line (@def_lines) {

        $line =~ s/\s+$//;
        my @values     = split( /\s/, $line );
        my $identifier = shift(@values);
        my $locus      = shift(@values);

        my ( $prefix, $cluster_id ) = split( /_/, $identifier );
        $hsh->{$cluster_id}->{locus} = $locus;
        $hsh->{$cluster_id}->{com_name} = join( " ", @values );

    }

    return $hsh;

}


sub process_db_file {

    my $file  = shift;
    my @lines = read_file($file);

    my $databases;
    my @db_array;

    foreach my $line (@lines) {

        $line =~ s/\s+$//;
        my ( $name, $location, $asmbl_id ) = split( /\t/, $line );

        $databases->{$name} = 0;
        push( @db_array, $name );

    }

    return ( $databases, \@db_array );

}


sub print_meta_groups {

    my ( $group_hsh, $clusters, $thresholds, $cluster_role_ids, $cluster_fgis, $fgis_hsh ) = @_;
    my $role_id_count_hsh;

    #Use filename as the output file
    my ( $filename, $path, $suffix ) = fileparse( $opts{metadata_file} );
    $filename =~ s/\..*$//;

    #Print data for each set of thresholds to their own file
    foreach my $true ( keys %$thresholds ) {

        my $false = $thresholds->{$true};
        my $fgi_print_output;

        open( my $mfh, ">", "$OUTPUT/$filename" . "_clusters_" . $true . "_" . $false . ".txt" );
        open( my $cfh, ">", "$OUTPUT/$filename" . "_counts_" . $true . "_" . $false . ".txt" );
        open( my $ifh, ">", "$OUTPUT/$filename" . "_fGIs_" . $true . "_" . $false . ".txt" ) if ( $opts{index} );

        print $cfh "Threshold parameter:\n";
        print $cfh "\tTrue Cutoff: $true" . "%\n";
        print $cfh "\tFalse Cutoff: $false" . "%\n";

        print $cfh
"\nGroup Thresholds\n*Number of required genomes to be present in a cluster for that cluster to be assigned the label\n\n";

        # Print the groups and the number of genomes from each group that
        # satisifes the threshold
        foreach my $g ( keys %$META_COUNTS ) {

            foreach my $l ( keys %{ $META_COUNTS->{$g} } ) {

                print $cfh "$g - $l\n";
                print $cfh "\tTrue cutoff: $META_COUNTS->{$g}->{$l}->{$true}->{true} out of $META_COUNTS->{$g}->{$l}->{original_count} genomes \n";
                print $cfh "\tFalse cutoff: $META_COUNTS->{$g}->{$l}->{$true}->{false} out of $META_COUNTS->{$g}->{$l}->{original_count} genomes\n";

                foreach my $genome ( sort @{ $META_LABELS->{$g}->{$l} } ) {
                    print $cfh "$genome\n";
                }

                print $cfh "\n";
            }
        }

        #Print Headers
        print $cfh "\nLabel Counts\t Number of Clusters\t Number of fGIs\n\n";
        print $mfh "Label\t";
        print $mfh "fGI ID\t" if $opts{index};
        print $mfh "Cluster ID\tRole ID\tProtein Name\tCentroid genome\tCentroid Locus\t";

        my $db_list_print = join( "\t", @dbs );
        print $mfh $db_list_print . "\n";

        # Prints out the clusters that match the threshold
        # Keeps track of the counts for each label
        foreach my $group ( keys %$group_hsh ) {

            my $amb_count = 0;

            foreach my $label ( sort keys %{ $group_hsh->{$group} } ) {

                #Clusters can have multiple labels, split the string
                my @labels = split( ",", $label );
                my $fgi_label_count;

                if ( exists $group_hsh->{$group}->{$label}->{$true} ) {

                    my @cluster_ids = @{ $group_hsh->{$group}->{$label}->{$true}->{cluster_ids} };
                    my @counts      = @{ $group_hsh->{$group}->{$label}->{$true}->{label_count} };

                    #Will ignore ambigous clusters to print unless flag is set
                    unless ( $label =~ /amb/ && !( $opts{print_amb} ) ) {

                        print $cfh "$label\t" . scalar @cluster_ids;
                        print $cfh "\n" unless ( $opts{index} );

                        for ( my $i = 0 ; $i < scalar @cluster_ids ; $i++ ) {

                            my $id = $cluster_ids[$i];
                            my @count = split( ",", $counts[$i] );
                            my @roles;

                            @roles = @{ $cluster_role_id->{$id} } if $cluster_role_id;

                            my ( $main, $subrole );

                            if ( scalar @roles > 0 ) {

                                foreach my $role (@roles) {

                                    $main    = $ROLE_ID_LOOKUP->{$role}->{main};
                                    $subrole = $ROLE_ID_LOOKUP->{$role}->{subrole};

                                    #Set values for those clusters with no assigned role_id
                                    if ( $role == 0 ) {
                                        $main    = "NONE";
                                        $subrole = "NONE";
                                    }

                                    # Stores the counts for individaual role categories
                                    $role_id_count_hsh->{$main}->{$label}->{$true}->{total}->{count}++;
                                    $role_id_count_hsh->{$main}->{$label}->{$true}->{$subrole}->{count}++;
                                }

                            } else {

                                $main    = "NONE";
                                $subrole = "NONE";
                                $role_id_count_hsh->{$main}->{$label}->{$true}->{total}->{count}++;
                                $role_id_count_hsh->{$main}->{$label}->{$true}->{$subrole}->{count}++;

                            }

                            #Creates the label/count printing string
                            my ( $print_label, $print_counts );

                            for ( my $i = 0 ; $i < scalar @labels ; $i++ ) {

                                my $label_string = join( "|", @labels );

                                if ( $label_string =~ /NULL/ ) {
                                    $print_label .= "$labels[$i]";
                                } else {
                                    $print_label  .= "$labels[$i],";
                                    $print_counts .= "$count[$i],";
                                }

                            }

                            $print_label =~ s/,$//;
                            $print_counts =~ s/,$//;

                            #Get centroid information
                            my $centroid_locus = $CENTROIDS->{$id}->{locus};
                            my $com_name       = $ANNOT_HSH{$centroid_locus}->{com_name};
                            my $centroid_db    = $ANNOT_HSH{$centroid_locus}->{db};

                            my $print_role_st = "";
                            if ( scalar @roles > 0 ) {

                                foreach my $role (@roles) {

                                    $print_role_st .= "$role," if $role ne '0';

                                }

                            } else {

                                $print_role_st = "-";

                            }

                            $print_role_st =~ s/,$//;

                            print $mfh "$print_label($print_counts)\t";

                            #print FGIs if option was passed in
                            if ( $opts{index} ) {

                                if ( exists $cluster_fgis->{$id} ) {

                                    $fgi_label_count->{$label}->{ $cluster_fgis->{$id} } = 1;
                                    print $mfh "$cluster_fgis->{$id}\t";

                                }

                            }

                            print $mfh "$id\t$print_role_st\t" . "$com_name\t$centroid_db\t$centroid_locus\t";
                            my $db_string;

                            foreach my $db (@dbs) {

                                if ( exists $clusters->{$id}->{$db} ) {
                                    print $mfh "$clusters->{$id}->{$db}\t";
                                    $db_string .= "$clusters->{$id}->{$db}\t";
                                } else {
                                    print $mfh "\t";
                                    $db_string .= "\t";
                                }
                            }

                            #Store information specific for fgr
                            $fgi_print_output->{$print_label}->{ $cluster_fgis->{$id} }->{$id} = 1 if ( $opts{index} );

                            print $mfh "\n";
                        }    #End cluster id loop

                        print $cfh "\t" . scalar( keys %{ $fgi_label_count->{$label} } ) . "\n" if $opts{index};

                    } else {

                        $amb_count = $amb_count + scalar @cluster_ids;

                    }

                }

            }

            print $cfh "Ambiguous Clusters\t$amb_count" . "\n";

            print_fgis( $fgi_print_output, $ifh, $cfh ) if ( $opts{index} );

            if ( $opts{role_lookup} ) {

                #print role_id files
                my @labels = keys %{ $group_hsh->{$group} };
                open( my $rfh, ">", "$OUTPUT/role_id_counts_" . $true . "_" . $false . ".txt" );

                #print headers
                print $rfh "role_id\tmain role\tsub role\t";
                my @print_labels;

                foreach my $label ( sort keys %{ $group_hsh->{$group} } ) {

                    my @labels = split( ",", $label );

                    unless ( $label =~ /amb/ && !( $opts{print_amb} ) ) {

                        push( @print_labels, $label );

                    }

                }

                map { print $rfh "$_\t" } sort(@print_labels);
                print $rfh "\n";

                foreach my $mainrole ( sort { $a cmp $b } keys $SUBROLE_LOOKUP ) {

                    my @subroles = keys %{ $SUBROLE_LOOKUP->{$mainrole} };

                    foreach my $role ( sort @subroles ) {

                        my $id = $SUBROLE_LOOKUP->{$mainrole}->{$role};

                        print $rfh "$id\t$mainrole\t$role\t";

                        foreach my $label ( sort @print_labels ) {

                            if ( exists $role_id_count_hsh->{$mainrole}->{$label}->{$true}->{$role} ) {

                                print $rfh $role_id_count_hsh->{$mainrole}->{$label}->{$true}->{$role}->{count} . "\t";

                            } else {

                                print $rfh "0\t";

                            }

                        }

                        print $rfh "\n";
                    }

                    #Print Mainrole total count line
                    print $rfh "MRSUM\t$mainrole\t\t";
                    foreach my $label ( sort @print_labels ) {

                        if ( exists $role_id_count_hsh->{$mainrole}->{$label}->{$true}->{total}->{count} ) {
                            print $rfh "$role_id_count_hsh->{$mainrole}->{$label}->{$true}->{total}->{count}" . "\t";
                        } else {
                            print $rfh "0\t";
                        }
                    }

                    print $rfh "\n";

                }

            }

        }

    }

}


sub print_fgis {

    my ( $print_hsh, $fh, $cfh ) = @_;

    #print header
    print $fh "Label\tfGI ID\t% of fGI\tfGI size\n";

    #print fGIs summary in count file for quick overview
    print $cfh "\nfGIs Summary\n";

    foreach my $label ( sort { $a cmp $b } keys %{$print_hsh} ) {

        foreach my $fgi ( sort { $a <=> $b } keys %{ $print_hsh->{$label} } ) {

            my $fgi_coverage    = scalar keys %{ $print_hsh->{$label}->{$fgi} };
            my $fgi_total_count = scalar keys %{ $F_FGIS->{$fgi} };

            my $perc = sprintf( "%.2f", $fgi_coverage / $fgi_total_count * 100 );

            print $cfh "\t$label: $perc% of fGI $fgi (Length: $fgi_total_count clusters)\n";
            print $fh "$label\t$fgi\t$perc\t$fgi_total_count\n";

        }

        print $cfh "\n";

    }

}


sub parse_cluster_results {

    my ($file) = @_;
    my $all_clusters;
    my $meta;
    my $role_ids;

    open( my $fh, "<", $file );

    while (<$fh>) {

        my @values = split( /\t/, $_ );
        map { $_ =~ s/\s+$// } @values;
        my @id = shift(@values);
        my $db_meta_match_hsh;
        my $role_id_found_hsh;

        # Determines what genome each locus is associated with
        # and stores in a hash what label this genome belongs to
        foreach my $value (@values) {

            if ( length($value) == 0 ) {

                #ignore blank columns which indicate there was no gene for this genome
                next;

            }

            #Find stored genome name
            my $db = $ANNOT_HSH{$value}->{db};

            if ( !defined($db) ) {

                #All gene identifiers in the table should have been defined in the attribute file
                die("$value was not defined in the attribute file $opts{data_file}\n");

            }

            #Store cluster information in hash
            $all_clusters->{ $id[0] }->{$db} = $value;

            #Determines what label is associated with the genome
            foreach my $group ( sort keys %{ $META_DATA->{$db} } ) {

                if ( exists $META_DATA->{$db}->{$group} ) {

                    foreach my $type ( keys %{ $META_DATA->{$db}->{$group} } ) {
                        push( @{ $db_meta_match_hsh->{$group}->{$type} }, $db );
                    }

                } else {

                    push( @{ $db_meta_match_hsh->{$group}->{"none"} }, $db );

                }
            }

            #If role_ids were provided, find associated role ids
            #Find role ids associated with current member
            if ( $opts{role_lookup} ) {

                if ( defined $ANNOT_HSH{$value}->{role_id} ) {

                    my $role_values = $ANNOT_HSH{$value}->{role_id};
                    my @role_ids = split( /\|/, $role_values );

                    foreach my $id (@role_ids) {

                        $id =~ s/\s+$//;
                        $id =~ s/^\s+//;
                        push( @{ $role_id_found_hsh->{$id} }, $value );

                    }

                } else {

                    push( @{ $role_id_found_hsh->{0} }, $value );

                }

            }

        }

        #store role_ids
        my @role_ids;
        if ( $opts{role_lookup} ) {

            @role_ids = keys %$role_id_found_hsh;
            $role_ids->{ $id[0] } = \@role_ids;

        }

        #Passes the genomes found in the cluster and their labels
        #to determine the final combined label for the cluster
        my $groups = find_meta_group($db_meta_match_hsh);

        foreach my $group ( keys %$groups ) {

            foreach my $thresh ( keys %{ $groups->{$group} } ) {

                my @labels = @{ $groups->{$group}->{$thresh} };

                my ( $label, $count );
                foreach my $l ( sort @labels ) {

                    my ( $p_l, $p_c ) = split( ":", $l );
                    $label .= "$p_l,";
                    $count .= "$p_c,";

                }

                $label =~ s/,$//;
                $count =~ s/,$//;

                push( @{ $meta->{$group}->{$label}->{$thresh}->{cluster_ids} }, $id[0] );
                push( @{ $meta->{$group}->{$label}->{$thresh}->{label_count} }, $count );

            }

        }

    }

    close $fh;
    return ( $all_clusters, $meta, $role_ids );

}


sub find_meta_group {

    my ($cluster_group_db_matches) = @_;
    my $groups_match;

    foreach my $group ( keys %$cluster_group_db_matches ) {

        #Go through each available label, if this particular
        #cluster represents genomes from this label use the
        #cutoffs to determine if it's true, false, or ambiguous
        foreach my $type ( sort keys %{ $META_COUNTS->{$group} } ) {

            foreach my $true_threshold (@$THRESHOLDS) {

                if ( exists $cluster_group_db_matches->{$group}->{$type} ) {

                    my $max_count_group     = $META_COUNTS->{$group}->{$type}->{$true_threshold}->{original_count};
                    my $cluster_count_group = scalar @{ $cluster_group_db_matches->{$group}->{$type} };

                    my $true_cutoff  = $META_COUNTS->{$group}->{$type}->{$true_threshold}->{true};
                    my $false_cutoff = $META_COUNTS->{$group}->{$type}->{$true_threshold}->{false};

                    if ( $cluster_count_group >= $true_cutoff ) {

                        push( @{ $groups_match->{$group}->{$true_threshold} }, "$type:$cluster_count_group" );

                    } elsif ( $cluster_count_group <= $false_cutoff ) {

                        push( @{ $groups_match->{$group}->{$true_threshold} }, "not_" . "$type:$cluster_count_group" );

                    } else {

                        push( @{ $groups_match->{$group}->{$true_threshold} }, "amb_" . "$type:$cluster_count_group" );
                    }

                } else {

                    push( @{ $groups_match->{$group}->{$true_threshold} }, "not_" . "$type:0" );

                }

            }

        }

    }

    return $groups_match;

}


sub parse_metadata {

    my ( $file, $threshold_cutoffs ) = @_;

    #Check File Existance
    die "$file does not exist or is size zero" unless ( -s $file );

    my @lines = read_file($file);
    my ( $hsh, $counts, $labels );

    #Store values in hash
    foreach my $line (@lines) {

        $line =~ s/\s+$//;
        my @values = split( /\t/, $line );

        if ( scalar @values >= 3 ) {

            #If there is a short hand name to use for genome, set that for display in out files
            my ( $db, $group, $type, $short ) = @values;
            my $use_label = ($short) ? $short : $type;

            $hsh->{$db}->{$group}->{$use_label} = 1;
            $counts->{$group}->{$use_label}->{original_count}++;

            push( @{ $labels->{$group}->{$use_label} }, $db );

        } else {
            die("Problem parsing metadata file (line: $line). Must have genome<tab>group<tab>label.");
        }

    }

    foreach my $group ( keys %$counts ) {

        foreach my $label ( keys %{ $counts->{$group} } ) {

            my $total = $counts->{$group}->{$label}->{original_count};

            #set cutoffs depending on threshold
            foreach my $true_key ( keys %$threshold_cutoffs ) {
                my ( $true, $false ) = calculate_true_false_values( $total, $threshold_cutoffs->{$true_key}, $true_key );

                $counts->{$group}->{$label}->{$true_key}->{true}  = $true;
                $counts->{$group}->{$label}->{$true_key}->{false} = $false;

            }

        }

    }

    return ( $hsh, $counts, $labels );

}


sub calculate_true_false_values {

    my ( $total, $false_key, $true_key ) = @_;

    my $true  = ( $total * $true_key ) / 100;
    my $false = ( $total * $false_key ) / 100;

    #Set floor/ceiling values
    if ( ( $false > 0 ) && ( $false < 1 ) ) {
        $false = 1;
    }

    if ( ( $true < $total ) && ( $true > ( $total - 1 ) ) ) {

        $true = $total - 1;
        $true = 1 if $true == 0;

    }

    #If true and false are the same value, make false one less
    #than true
    if ( $true == $false ) {
        $false = $true - 1;
    }

    return ( $true, $false );

}


sub set_default_cutoffs {

    my $hsh;

    $hsh->{100} = 0;
    $hsh->{50}  = 0;
    $hsh->{90}  = 10;
    $hsh->{75}  = 25;

    return $hsh;

}


sub parse_role_id_lookup {

    my $file = shift;

    my $errors .= "File is size zero or does not exist: file\n" unless ( -s $file );
    my ( $hsh, $shsh );

    if ($errors) {
        die($errors);
    } else {

        my @role_lookup = read_file($file);

        foreach (@role_lookup) {
            my @values = split( /\t/, $_ );
            $values[1] =~ s/\s+$//;
            $values[2] =~ s/\s+$//;

            $hsh->{ $values[0] }->{main}    = $values[1];
            $hsh->{ $values[0] }->{subrole} = $values[2];

            $shsh->{ $values[1] }->{ $values[2] } = $values[0];

        }

    }

    return ( $hsh, $shsh );
}


sub check_params {

    my $default_cutoffs = shift;
    my ( $thresholds, $output, $errors, $threshold_cutoffs );

    #Make output directory if it doesn't exist
    $output = $opts{output} // cwd();
    make_path($output) unless ( -d $output );

    #Check thresholds and store cutoffs
    if ( $opts{metadata_file} ) {

        if ( -s $opts{metadata_file} ) {

            if ( $opts{threshold_file} ) {

                #Set thresholds based on file specifications
                if ( -s $opts{threshold_file} ) {

                    my @lines = read_file( $opts{threshold_file} );

                    foreach my $line (@lines) {

                        $line =~ s/\s+$//;
                        my ( $true, $false ) = split( /\t/, $line );
                        $threshold_cutoffs->{$true} = $false;
                        push( @$thresholds, $true );

                    }

                } else {
                    $errors .= "File does not exist or is size zero: $opts{threshold_file}\n";
                }

            } else {

                my @default = keys %$default_cutoffs;
                $thresholds = \@default;

                $threshold_cutoffs = $default_cutoffs;
            }

        } else {

            $errors .= "No metadata file provided or it's size zero/does not exist\n";

        }

    }

    return ( $thresholds, $output, $threshold_cutoffs );

}
