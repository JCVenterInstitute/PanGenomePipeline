#!/usr/local/bin/perl
#! /usr/local/bin/perl -w

=head1 NAME

method_comparison.pl

=head1 SYNOPSIS

    USAGE: cluster_comparison.pl --output_dir <output directory> 
                                  -U <sybase username>
                                  -P <sybase password>
                                  [--db_list <list.txt> || --fasta_files <file.txt>] 
                                  [--method_files <file.txt> || --pan_genome_db || method_result_file_1 method_result_file_2 ...] 
                                  [--method_comparison]
=head1 OPTIONS

B<--db_list_file,-d>    :   Path to a file containing one sgd name per line

B<--username,-U>        :   Sybase username

B<--password,-P>        :   Sybase password

B<--fasta_files, -f>    :   Path to a file containing full path to individual fasta files, one per line (Used for data not in a SGD)

B<--method_files,-m>    :   Path to file containing paths to cluster method output, one per line.

B<--pan_db, -g>         :   PanGenome db. (Require if no --method_files)

B<--method_comparison>  :   Will do a method by method comparison

B<--output_dir, -o>     :   Output directory for results

B<--help,-h>            :   Display help

=head1 ARGUMENTS
	
The individual comparison method results files. Output of pangenome_cluster.pl.
Not needed if you provide --method_list or --pan_db.
	
=head1  DESCRIPTION

This script compares the clusters from multiple clustering methods and generates standard statistics. 

=head1  OUTPUT

This script currently outputs the following files:

    comparison_differences.txt - File listing all the clustering differences between the methods
    comparison_similarities.txt - File listing the clusters that are identical across all methods
    comparison_perfect_clusters.txt - File listing clusters with one and only one representative for each genome

This script also outputs the files above for the two way comparison method clusters. For example:

    comparison_differences_btm_vs_sybil.txt
    comparison_similarities_btm_vs_sybil.txt
    comparison_perfect_clusters_btm_vs_sybil.txt

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use warnings;
use strict;
$|++;

use File::Path;
use File::Spec;
use FindBin;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use Getopt::Long;
use Cwd 'realpath';
use Cwd;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use DBI;
use File::Slurp;

my $program = realpath($0);
my $myLib   = dirname($program);

my ( $user, $password );
my @FILES;

my %opts;
my $results = GetOptions( \%opts, 'output_dir|o=s',
				  'user|U=s',        'password|P=s',
				  'db_list|d=s',     'method_files|m=s',
				  'fasta_files|f=s', 'pan_db|g=s',
				  'method_comparison',
				  'help|h',) || pod2usage();

my $DBTYPE             = "Sybase";
my $SERVER             = "SYBPROD";
my $DEFAULT_USER       = "access";
my $DEFAULT_PASS       = "access";
my $DEFAULT_OUTPUT_DIR = getcwd();
my $PAN_SERVER         = "SYBIL";

my $output_dir = '';
my @COMPARISON_METHODS;
my $frame_shift = 0;

&check_options();

my ($PG_DBPROC) = &ConnectToDb( $PAN_SERVER, $DBTYPE, $user, $password, $opts{pan_db} )
	if ( $opts{pan_db} );

# These store the cluster results information
my ( %METHOD_CLUSTERS_LOCI, %METHOD_CLUSTERS_ID );

# Total count information
## $TOTAL_MEMBERS stores the loci list per genome of locus found in the cluster files
my ( $TOTAL_CLUSTER_COUNTS, $TOTAL_MEMBERS );

#Finds the comparison methods based on parameters passed in
# $method_hsh_db is only returned if pulling results from the pangenom db
my $method_hsh_db = &find_comparison_methods;

## Create look up hashes that will later be used in comparison
## $DB_LOCI contains loci list from database for each db passed in
## $LOCI_DB_LOOKUP contains db, feat_name and hmm information for each loci in a genome's database
## $EV_LOOKUP contains each HMM & PRK info and what loci match them as well as the com_name
my ( $DB_LOCI, $LOCI_DB_LOOKUP, $EV_LOOKUP ) = &get_loci( $opts{'db_list'}, $opts{'fasta_files'} );

## Parse each cluster comparison file and store results
&parse_pangenome_files( \@FILES ) if @FILES;
&parse_pangenome_db( $opts{pan_db}, $method_hsh_db ) if $opts{pan_db};

&print_total_counts;

## Store the comparison methods names, parsed from the files.
my @genomes = keys %$DB_LOCI;

if ( scalar @FILES >= 2 || $opts{pan_db} ) {

	# These store the comparison results
	my ( $diff, $perfect, $same );

	## Compare all methods together
	open( SINGLE, ">$output_dir/singletons.txt" );
	foreach my $db (@genomes) {
		( $perfect, $same, $diff )
			= &get_comparision_for_one_genome_four_methods( $db, $diff, $perfect, $same );
	}

	close SINGLE;
	## Do prints for all method comparison
	&print_stats( $perfect, $same );
	&print_files( $same, $diff );

    if($opts{method_comparison}){
		## Compare each method against every other method
		for ( my $i = 0; $i < $#COMPARISON_METHODS; $i++ ) {
	
			for ( my $j = $i + 1; $j <= $#COMPARISON_METHODS; $j++ ) {
	
				print
					"\n\nDetailed comparison of $COMPARISON_METHODS[$i] and $COMPARISON_METHODS[$j]\n";
	
				my ( $same, $diff, $perfect );
				my $single_file
					= $output_dir
					. "/singletons_"
					. $COMPARISON_METHODS[$i] . "_"
					. $COMPARISON_METHODS[$j] . ".txt";
	
				open( SINGLE, ">$single_file" );
				foreach my $db (@genomes) {
	
					( $perfect, $same, $diff )
						= &get_comparision_for_one_genome_two_methods( $db, $COMPARISON_METHODS[$i],
													  $COMPARISON_METHODS[$j], $diff, $perfect, $same );
				}
				close SINGLE;
	
				&print_stats( $perfect, $same );
				&print_files( $same, $diff, $COMPARISON_METHODS[$i], $COMPARISON_METHODS[$j] );
			}
	
		}
    }
}

####################### SUBS #############################
sub find_comparison_methods {

	&find_comparison_methods_file if (@FILES);
	my $method_hsh_db = &find_comparison_methods_db if ( $opts{pan_db} );

	return ($method_hsh_db);
}

sub find_comparison_methods_file {
	foreach my $input (@FILES) {

		$input =~ s/\s+$//;

		my ( $name, $path, $suffix ) = fileparse($input);

		push( @COMPARISON_METHODS, $name );
	}
}

sub find_comparison_methods_db {
	my $method_hsh;

	my $query = "SELECT type, method_id " . "FROM method ";

	my @results = &do_sql( $PG_DBPROC, $query );

	foreach my $method (@results) {
		my ( $type, $id ) = split( /\t/, $method );

		push( @COMPARISON_METHODS, $type );

		$method_hsh->{$type} = $id;
	}

	return ($method_hsh);
}

sub find_cluster_ids {
	my ( $id, $method ) = @_;

	my $query = "SELECT cluster_id, cluster_name " . "FROM cluster " . "WHERE method_id = $id";

	my @results = &do_sql( $PG_DBPROC, $query );
	my $result_hsh;

	foreach my $result (@results) {
		my ( $c_id, $c_name ) = split( /\t/, $result );
		$result_hsh->{$method}->{$c_id} = $c_name;
	}

	return ($result_hsh);
}

sub find_cluster_members {
	my $id = shift;

	my $query
		= "SELECT m.locus,g.db "
		. "FROM cluster_members cm, members m, genomes_data g "
		. "WHERE cm.member_id = m.member_id "
		. "AND m.genomes_data_id = g.genomes_data_id "
		. "AND cm.cluster_id = $id";

	my @results = &do_sql( $PG_DBPROC, $query );

	my $result_hsh;

	foreach my $result (@results) {
		my ( $locus, $db ) = split( /\t/, $result );
		$result_hsh->{$locus} = $db;
	}

	return ($result_hsh);
}

sub get_comparision_for_one_genome_one_method {

	#Find stats for this one genome in all comparison methods
	my ( $db, $PERFECT_HSH_BY_METHOD, $print_hsh, $non_perfect_hsh ) = @_;

	my $db_count   = scalar keys %$DB_LOCI;
	my $loci_count = scalar( keys %{ $DB_LOCI->{$db} } ) - 1;

	# Goes through each locus of a genome/database
	foreach my $locus ( keys %{ $DB_LOCI->{$db} } ) {

		my $current_loci_hsh;

		if ( $locus ne 'count' ) {

			my $cluster_locations;

			foreach my $key (@COMPARISON_METHODS) {

				# Only add if a cluster for this loci exists for this method
				if ( exists $METHOD_CLUSTERS_LOCI{$key}->{$db}->{$locus}->{'orthoids'} ) {

					my $cluster = $METHOD_CLUSTERS_LOCI{$key}->{$db}->{$locus}->{'orthoids'};

					$cluster_locations->{$key} = $cluster;

				}

			}

			my $not_perfect_cluster
				= 0;    #Flag to see if cluster has no paralogs and a rep from each DB

			my $cluster_loci_lookup;

			# Merge loci for comparison
			# Go through each method that has a cluster containing the current locus
			foreach my $method ( keys %$cluster_locations ) {

				my $combined_loci;
				my $cluster_id     = $cluster_locations->{$method};
				my $method_perfect = 0;

				$cluster_id =~ s/\|//;

				# Merge all the loci from the cluster
				# Go through array to preserve DB order
				foreach my $db (@genomes) {

					my $key_lookup = $db . "_loci";

					# Check to see that loci exists in this cluster for this DB
					if ( exists $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup} ) {

						my $current_loci
							= $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup};
						$combined_loci .= $current_loci . "|";

			 # Set flag that it's not perfect because it ends in a pipe, meaning more than one locus
			 # per genome
						my @members = split( '\|', $current_loci );

						if ( scalar @members > 1 ) {
							$not_perfect_cluster = 1;
							$method_perfect      = 1;
						}
						else {
							$non_perfect_hsh->{$cluster_id}->{$db} = scalar @members;
						}
					}
					else {

						# Set flag indicating not perfect because loci are not present for each db
						$not_perfect_cluster = 1;
						$method_perfect      = 1;

						#$non_perfect_hsh->{$cluster_id}->{$db} = "no hits";
						#$non_perfect_hsh->{$cluster_id}->{'no_hit_count'}++;
					}
				}

				unless ($method_perfect) {
					$PERFECT_HSH_BY_METHOD->{$method}->{$cluster_id} = 1;
				}

				$cluster_loci_lookup->{$method}->{'loci'} = $combined_loci;

			}

		}
	}

	my $db_loci_count     = scalar keys %{ $TOTAL_MEMBERS->{$db} };
	my $no_cluster_number = $DB_LOCI->{$db}->{'count'} - $db_loci_count;

	$print_hsh->{$db}->{'loci'}       = $DB_LOCI->{$db}->{'count'};
	$print_hsh->{$db}->{'no_cluster'} = $no_cluster_number;

	return ( $PERFECT_HSH_BY_METHOD, $print_hsh, $non_perfect_hsh );
}

sub get_comparision_for_one_genome_four_methods {

	#Find stats for this one genome in all comparison methods
	my ( $db, $diff, $perfect, $same ) = @_;
	my $same_fs_hsh;

	my @headers = join( "\t", @COMPARISON_METHODS );

	my $db_count   = scalar keys %$DB_LOCI;
	my $loci_count = scalar( keys %{ $DB_LOCI->{$db} } ) - 1;

	my ( $same_count, $majority_method_count,, $no_loci_match_hsh );

	( $same_count, $majority_method_count, $perfect, $same, $diff, $no_loci_match_hsh )
		= &method_comparison( $db, $diff, $perfect, $same );

	print "\n\n$db\n";
	print "\tNumber of loci in the genome: $DB_LOCI->{$db}->{'count'}\n";

	my $db_loci_count = scalar keys %{ $TOTAL_MEMBERS->{$db} };

	my $no_cluster_number = $DB_LOCI->{$db}->{'count'} - $db_loci_count;

	print "\tNumber of loci without clusters: $no_cluster_number\n";
	print "\tNumber of loci with the same cluster results from all "
		. scalar @COMPARISON_METHODS
		. " methods: $same_count\n";

	if ( scalar @COMPARISON_METHODS > 2 ) {
		print "\tNumber of loci with the same cluster results for "
			. scalar @COMPARISON_METHODS - 1
			. " out of the "
			. scalar @COMPARISON_METHODS
			. " methods: $majority_method_count\n";
	}

	return ( $perfect, $same, $diff );
}

sub get_comparision_for_one_genome_two_methods {
	my ( $db, $method_1, $method_2, $diff, $perfect, $same ) = @_;

	my ( $prefix_1, $postfix_1 ) = split( /\./, $method_1 );
	my ( $prefix_2, $postfix_2 ) = split( /\./, $method_2 );

	my @comparison_methodS;
	push( @comparison_methodS, $method_1 );
	push( @comparison_methodS, $method_2 );

	my $majority_method_count;
	my $no_loci_match_hsh;

	my $same_count;

	( $same_count, $majority_method_count, $perfect, $same, $diff, $no_loci_match_hsh )
		= &method_comparison( $db, $diff, $perfect, $same );

	my $db_loci_count = scalar keys %{ $TOTAL_MEMBERS->{$db} };

	my $no_cluster_number;
	my $loci_number;

	print "\n\n$db\n";
	print "\tNumber of loci in the genome: $DB_LOCI->{$db}->{'count'}\n";

	print "\tNumber of loci with the same cluster results from all "
		. scalar @comparison_methodS
		. " methods: $same_count\n";

	my $method_1_count = 0;
	my $method_2_count = 0;

	foreach my $loci ( keys %$no_loci_match_hsh ) {
		if ( exists $no_loci_match_hsh->{$loci}->{$method_1} ) {
			unless ( exists $no_loci_match_hsh->{$loci}->{$method_2} ) {
				$method_1_count++;
			}
		}
		elsif ( exists $no_loci_match_hsh->{$loci}->{$method_2} ) {
			unless ( exists $no_loci_match_hsh->{$loci}->{$method_1} ) {
				$method_2_count++;
			}
		}
	}

	print
		"\n\tNumber of loci in a $method_1 clusters where no $method_2 cluster exists: $method_1_count\n";
	print
		"\tNumber of loci in a $method_2 clusters where no $method_1 cluster exists: $method_2_count\n";

	return ( $perfect, $same, $diff );
}

sub method_comparison {
	my ( $db, $diff, $perfect, $same ) = @_;

	my $same_genome;
	my $majority_method_count;
	my $no_loci_match_hsh;
    my $frame_shift;
    
	# Goes through each locus of a genome/database
	foreach my $locus ( keys %{ $DB_LOCI->{$db} } ) {

		my $current_loci_hsh;

		if ( $locus ne 'count' ) {

			my $cluster_locations;
			my @singleton;

			foreach my $key (@COMPARISON_METHODS) {

				# Only add if a cluster for this loci exists for this method
				if ( exists $METHOD_CLUSTERS_LOCI{$key}->{$db}->{$locus}->{'orthoids'} ) {

					my $cluster = $METHOD_CLUSTERS_LOCI{$key}->{$db}->{$locus}->{'orthoids'};

					my @ids = split( /\|/, $cluster );

					foreach my $id (@ids) {
						if ( $METHOD_CLUSTERS_ID{$key}->{$id}->{'count'} > 1 ) {

							$cluster_locations->{$key} .= $id . "|";
						}
					}

				}
				else {
					$no_loci_match_hsh->{$locus}->{$key} = 1;
					$no_loci_match_hsh->{$locus}->{'method_count'}++;

					push( @singleton, $key );
				}
			}

			#Print singletons (loci that didn't cluster in any method)
			print SINGLE "$locus\n" if ( scalar @singleton == scalar @COMPARISON_METHODS );

			my $not_matching_flag = 0;    #Flag set to 1 if the loci mismatch
			my $not_perfect_cluster
				= 0;    #Flag to see if cluster has no paralogs and a rep from each DB

			my $cluster_loci_lookup;

			# Merge loci for comparison
			# Go through each method that has a cluster containing the current locus
			foreach my $method ( keys %$cluster_locations ) {

				my $cluster_id     = $cluster_locations->{$method};
				my $method_perfect = 0;

				my @ids = split( /\|/, $cluster_id );

				# Merge all the loci from the cluster
				# Go through array to preserve DB order

				foreach my $id (@ids) {
					my $combined_loci;

					foreach my $db (@genomes) {

						my $key_lookup = $db . "_loci";

						# Check to see that loci exists in this cluster for this DB
						if ( exists $METHOD_CLUSTERS_ID{$method}->{$id}->{$key_lookup} ) {

							my $current_loci = $METHOD_CLUSTERS_ID{$method}->{$id}->{$key_lookup};
							$combined_loci .= $current_loci . "|";

							# Set flag that it's not perfect more than one locus per genome
							my @members = split( '\|', $current_loci );

							$not_perfect_cluster = 1 if scalar @members > 1;
							$method_perfect      = 1 if scalar @members > 1;

						}
						else {

						  # Set flag indicating not perfect because loci are not present for each db
							$not_perfect_cluster = 1;
							$method_perfect      = 1;

						}
					}

					unless ($method_perfect) {
						$perfect->{$method}->{$id} = 1;
					}

                   # if(exists $LOCI_DB_LOOKUP->{$locus}->{'fs'} ) {
                    #	$frame_shift->{$method}->{$id} = 'fs';
                    #}
                    
					$cluster_loci_lookup->{$method}->{'loci'}->{$id} = $combined_loci;

				}

			}

			my $method_match_count = 0;
			my $method_count;

			for ( my $i = 0; $i < $#COMPARISON_METHODS; $i++ ) {
				my $method = $COMPARISON_METHODS[$i];

				if ( exists $cluster_loci_lookup->{$method}->{'loci'} ) {

					foreach my $id ( keys %{ $cluster_loci_lookup->{$method}->{'loci'} } ) {

						my $loci = $cluster_loci_lookup->{$method}->{'loci'}->{$id};

						for ( my $j = $i + 1; $j <= $#COMPARISON_METHODS; $j++ ) {

							my $comp_method = $COMPARISON_METHODS[$j];

							if ( exists $cluster_loci_lookup->{$method}->{'loci'} ) {

								foreach my $id (
										  keys %{ $cluster_loci_lookup->{$comp_method}->{'loci'} } )
								{
									my $comp_loci = $cluster_loci_lookup->{$comp_method}->{'loci'}->{$id};

                                    
									if ( $loci eq $comp_loci ) {
										
										
										$method_count->{$loci}->{$comp_method} = 1;
										$method_count->{$loci}->{$method}      = 1;
										
										
									}
									else {
									
										$not_matching_flag = 1;

									}

								}

							}
						}
					}

				}
				else {

					next;

				}
			}

			if ( $method_count && scalar @COMPARISON_METHODS > 2 ) {
				foreach my $loci_list ( keys %$method_count ) {
					my $method_matched = scalar keys %{ $method_count->{$loci_list} };

					$majority_method_count++
						if ( $method_matched == ( scalar @COMPARISON_METHODS - 1 ) );
				}
			}

                    
			foreach my $method (@COMPARISON_METHODS) {

				if ($not_matching_flag) {
					if ( exists $cluster_locations->{$method} ) {

						foreach my $id ( keys %{ $cluster_loci_lookup->{$method}->{'loci'} } ) {
							my $loci = $cluster_loci_lookup->{$method}->{'loci'}->{$id};
							my @loci = split( /\|/, $loci );

							foreach my $locus (@loci) {
								$current_loci_hsh->{'loci'}->{$locus} = 1;
							}

							$current_loci_hsh->{'method'}->{$method}
								= $cluster_locations->{$method};
						}
					}
					else {
						$current_loci_hsh->{'method'}->{$method} = " ";
					}

				}
				else {
				    my $method_count = scalar( keys %$cluster_locations );       
					
					if($method_count > 1){
						foreach my $id ( keys %{ $cluster_loci_lookup->{$method}->{'loci'} } ) {
	
							my $same_loci = $cluster_loci_lookup->{$method}->{'loci'}->{$id};
	
							if ( exists $LOCI_DB_LOOKUP->{$locus}->{'fs'} ) {
								$same->{$same_loci}->{'fs'} = 1;
							}
	
							if ( $method_count == scalar @COMPARISON_METHODS ) {
	
								if ( exists $cluster_locations->{$method} ) {
									$same_genome->{$locus}->{$method} = $cluster_locations->{$method};
								}
	
								unless ( exists $same->{$same_loci} ) {
									$same->{'count'}++;
								}
	
								$same->{$same_loci}->{$method} = $cluster_locations->{$method}
									if ( exists $cluster_locations->{$method} );
	
								unless ($not_perfect_cluster) {
	    
	                                      
									if ( exists $cluster_locations->{$method} ) {
										$same->{$same_loci}->{'perfect'} = 1;
									}
	
								}
							}
						}
					}
				}

			}
			if ($current_loci_hsh) {
				my @cluster_loci = sort( keys %{ $current_loci_hsh->{'loci'} } );
				my $combined_loci = join( '|', @cluster_loci );

				my $combined_method;
				my @method = sort( keys %{ $current_loci_hsh->{'method'} } );

				foreach my $line (@method) {
					$combined_method .= "$current_loci_hsh->{'method'}->{$line}";
				}

				foreach my $key ( keys %{ $current_loci_hsh->{'method'} } ) {
					$diff->{$combined_method}->{$key} = $current_loci_hsh->{'method'}->{$key};
				}
			}
		}
	}

	my $same_count = scalar keys %$same_genome;

	return ( $same_count, $majority_method_count, $perfect, $same, $diff, $no_loci_match_hsh );

}

sub print_total_counts {
	print "\n";

	foreach my $method ( keys %$TOTAL_CLUSTER_COUNTS ) {
		print "Total $method clusters: $TOTAL_CLUSTER_COUNTS->{$method}\n";
	}
}

sub print_stats {
	my ( $perfect, $same ) = @_;

	my $cluster_id_count_hsh;
	my $total_count = 0;
    my $fs_count = 0;
    
	print "\nTotal Cluster Numbers\n";

	my $check_count = 0;

	foreach my $loci ( keys %$same ) {
		if ( $loci ne 'count' ) {
			$total_count++ if exists $same->{$loci}->{'perfect'};
            $fs_count++ if exists $same->{$loci}->{'fs'};
			
		}
	}

	print "Number of clusters where all methods agree: $same->{'count'}\n";
	print "Number of clusters where there are "
		. scalar @genomes
		. " members, one per genome and all methods agree: $total_count\n";
    print "Number of clusters where all methods agree and some loci are frameshifted: $fs_count\n\n";
    
	foreach my $method ( keys %$perfect ) {
		my $count = scalar keys %{ $perfect->{$method} };
		print "Number of clusters in $method only where there are "
			. scalar @genomes
			. " members and all methods agree: $count\n";
	}

}

sub pull_loci_from_db {
	my $db_list = shift;

	open( FILE, $db_list );

	my ( $db_loci, $loci_db_lookup, $ev_lookup );

	foreach my $genome (<FILE>) {

		$genome =~ s/\s+$//;

		my ($dbproc) = &ConnectToDb( $SERVER, $DBTYPE, $user, $password, $genome );

		my $query
			= "select i.locus, i.feat_name "
			. "FROM ident i, asm_feature a, stan s "
			. "WHERE i.feat_name = a.feat_name "
			. "AND a.asmbl_id = s.asmbl_id "
			. "AND s.iscurrent = 1 "
			. "AND a.feat_type = \"ORF\"";

		my @loci_results = &do_sql( $dbproc, $query );

		my $count = 0;
		$db_loci->{$genome}->{'count'} = scalar @loci_results;

		foreach my $loci (@loci_results) {
			my ( $locus, $feat_name ) = split( /\t/, $loci );
			$db_loci->{$genome}->{$locus} = 1;

			$loci_db_lookup->{$locus}->{'db'}        = $genome;
			$loci_db_lookup->{$locus}->{'feat_name'} = $feat_name;

			if ( scalar @COMPARISON_METHODS > 1 ) {
				my $query
					= "select e.accession, h.hmm_com_name "
					. " from $genome..evidence e, $genome..feat_score f, hmm..hmm3 h "
					. "where h.hmm_acc = e.accession "
					. " and e.feat_name = \"$feat_name\" "
					. " and h.iso_type in (\"equivalog\", \"PFAM_equivalog\")  "
					. " and f.input_id = e.id "
					. " and f.score_id = 51  "
					. " and convert(numeric, f.score) > h.trusted_cutoff";

				my @hmm_results = &do_sql( $dbproc, $query );

				my $frame_query
					= "select f.att_type "
					. "FROM $genome..ident i, $genome..frameshift f "
					. "WHERE i.feat_name = f.feat_name "
					. "AND i.locus = \"$locus\"";
				my @frame_results = &do_sql( $dbproc, $frame_query );

				my $prk_query
					= "SELECT e.accession, f.score "
					. "FROM evidence e, feat_score f, common..score_type s "
					. "WHERE s.id = f.score_id "
					. "AND e.ev_type = s.input_type "
					. "AND e.curated = 1 "
					. "AND e.ev_type = \"PRK\" "
					. "AND e.feat_name = \"$feat_name\" "
					. "AND s.id=98014 ";

				my @prk_results = &do_sql( $dbproc, $prk_query );
				if ( scalar @hmm_results > 0 ) {

					foreach my $results (@hmm_results) {
						my ( $acc, $com_name ) = split( /\t/, $results );
						$ev_lookup->{'hmm'}->{$acc}->{'feat_name'}->{$locus} = 1;
						$ev_lookup->{'hmm'}->{$acc}->{'com_name'}            = "$com_name";
						$loci_db_lookup->{$locus}->{'hmm'}->{$acc}           = 1;
					}
				}

				if ( scalar @frame_results > 0 ) {
					foreach my $results (@frame_results) {
						if ( $results =~ /(FS|PM|AMB|DEG)/ ) {
							$loci_db_lookup->{$locus}->{'fs'}->{$results} = 1;
						}
					}

				}

				if ( scalar @prk_results > 0 ) {
					foreach my $results (@prk_results) {
						my ( $acc, $com_name ) = split( /\t/, $results );
						$ev_lookup->{'prk'}->{$acc}->{'feat_name'}->{$locus} = 1;
						$ev_lookup->{'prk'}->{$acc}->{'com_name'}            = "$com_name";
						$loci_db_lookup->{$locus}->{'prk'}->{$acc}           = 1;
					}
				}
			}

		}

		$dbproc->disconnect;

	}

	return ( $db_loci, $loci_db_lookup, $ev_lookup );
}

sub pull_loci_from_fasta {
	my $fasta_file = shift;

	my @files = read_file($fasta_file);
	my ( $db_loci, $loci_db_lookup, $ev_lookup );

	foreach my $fasta (@files) {

		$fasta =~ s/\s+$//;

		my ( $genome, $path, $suffix ) = fileparse($fasta);

		my @deflines = `grep '>' $fasta`;

		if (@deflines) {
			my $count = scalar @deflines;
			$db_loci->{$genome}->{'count'} = $count;
		}
		else {
			die "Could not run grep -c '>' $fasta, make sure file exists";
		}

		foreach my $locus (@deflines) {
			my ( $accession, $header ) = split( /\s/, $locus );

			$accession =~ s/>//;
			$accession =~ s/\s+$//;

			$db_loci->{$genome}->{$accession} = 1;

			$loci_db_lookup->{$accession}->{'db'} = $genome;
		}

	}

	return ( $db_loci, $loci_db_lookup );

}

sub get_loci {
	my ( $db_list, $fasta_list ) = @_;

	my ( $db_loci, $loci_db_lookup, $ev_lookup );

	( $db_loci, $loci_db_lookup, $ev_lookup ) = &pull_loci_from_db($db_list) if $db_list;
	( $db_loci, $loci_db_lookup ) = &pull_loci_from_fasta($fasta_list) if $fasta_list;

	return ( $db_loci, $loci_db_lookup, $ev_lookup );
}

sub parse_pangenome_db {
	my ( $db, $method_hsh_db ) = @_;

	foreach my $method ( keys %$method_hsh_db ) {
		my $method_id = $method_hsh_db->{$method};

		my %db_method_id;

		my $id_results = &find_cluster_ids( $method_id, $method );

		foreach my $cluster_id ( keys %{ $id_results->{$method} } ) {

			my $cluster_results = &find_cluster_members($cluster_id);
			my @cluster_members = keys %$cluster_results;
			my $cluster_count   = scalar @cluster_members;

			if ( $cluster_count > 1 ) {    #accounts for BTM listing singletons as clusters

				$METHOD_CLUSTERS_ID{$method}->{ $id_results->{$method}->{$cluster_id} }->{'count'}
					= $cluster_count;

				foreach my $loci (@cluster_members) {

					$loci =~ s/\s+$//;

					my $db = $cluster_results->{$loci};

					$db_method_id{$db}->{'count'}++;

					$METHOD_CLUSTERS_ID{$method}->{ $id_results->{$method}->{$cluster_id} }
						->{ $db . "_loci" } .= "$loci|";
					$METHOD_CLUSTERS_LOCI{$method}->{$db}->{$loci}->{'orthoids'}
						.= "$id_results->{$method}->{$cluster_id}|";

					$TOTAL_MEMBERS->{$db}->{$loci} = 1;

				}

				foreach my $db ( keys %db_method_id ) {
					$METHOD_CLUSTERS_ID{$method}->{ $id_results->{$method}->{$cluster_id} }
						->{ $db . "_count" } = $db_method_id{$db}->{'count'};
				}

			}

		}
		foreach my $key ( sort_multihashkeys_by_value( $METHOD_CLUSTERS_ID{$method}, "count" ) ) {

			foreach my $value ( keys %{ $METHOD_CLUSTERS_ID{$method}->{$key} } ) {

				if ( $value =~ /\_loci$/ ) {
					$METHOD_CLUSTERS_ID{$method}->{$key}->{$value} =~ s/\|$//;
					&sort_pipes( $METHOD_CLUSTERS_ID{$method}, $key, $value );
				}

			}

		   #Only count clusters with more than one member, some methods count singletons as clusters
			$TOTAL_CLUSTER_COUNTS->{$method}++
				if $METHOD_CLUSTERS_ID{$method}->{$key}->{'count'} > 1;
		}

	}

}

sub parse_pangenome_files {
	my $files = shift;

	foreach my $file (@$files) {
		$file =~ s/\s+$//;

		my ( $name, $path, $suffix ) = fileparse($file);

		open( FILE, $file ) || die "can't open $file: $!";

		while (<FILE>) {
			my $cluster_line = $_;

			my %db_file;

			my @results         = split( /\t/, $cluster_line );
			my $cluster_id      = $results[0];
			my @cluster_members = splice( @results, 1, scalar @results );
			my $cluster_count   = scalar @cluster_members;

			if ( $cluster_count > 1 ) {    #accounts for BTM listing singletons as clusters

				$METHOD_CLUSTERS_ID{$name}->{$cluster_id}->{'count'} = $cluster_count;

				foreach my $loci (@cluster_members) {

					$loci =~ s/\s+$//;

					my $db = $LOCI_DB_LOOKUP->{$loci}->{'db'};

					$db_file{$db}->{'count'}++;

					$METHOD_CLUSTERS_ID{$name}->{$cluster_id}->{ $db . "_loci" } .= "$loci|";
					$METHOD_CLUSTERS_LOCI{$name}->{$db}->{$loci}->{'orthoids'}   .= "$cluster_id|";

					$TOTAL_MEMBERS->{$db}->{$loci} = 1;

				}

				foreach my $db ( keys %db_file ) {
					$METHOD_CLUSTERS_ID{$name}->{$cluster_id}->{ $db . "_count" }
						= $db_file{$db}->{'count'};
				}

			}
		}

		foreach my $key ( sort_multihashkeys_by_value( $METHOD_CLUSTERS_ID{$name}, "count" ) ) {

			foreach my $value ( keys %{ $METHOD_CLUSTERS_ID{$name}->{$key} } ) {

				if ( $value =~ /\_loci$/ ) {
					$METHOD_CLUSTERS_ID{$name}->{$key}->{$value} =~ s/\|$//;
					&sort_pipes( $METHOD_CLUSTERS_ID{$name}, $key, $value );
				}

			}

		   #Only count clusters with more than one member, some methods count singletons as clusters
			$TOTAL_CLUSTER_COUNTS->{$name}++ if $METHOD_CLUSTERS_ID{$name}->{$key}->{'count'} > 1;
		}

		close FILE;
	}
}

sub print_files {
	my ( $SAME_HSH, $DIFF_HSH, $method_1, $method_2 ) = @_;

	if ( $method_1 && $method_2 ) {

		$method_1 =~ s/\.\w+//;
		$method_2 =~ s/\.\w+//;

		&print_diff_hsh( $DIFF_HSH, "comparison_differences_" . $method_1 . "_vs_" . $method_2 );
		&print_same_cluster_hsh( $SAME_HSH,
								 "comparison_similarities_" . $method_1 . "_vs_" . $method_2 );
		&print_perfect_cluster_hsh( $SAME_HSH,
								  "comparison_perfect_clusters_" . $method_1 . "_vs_" . $method_2 );
	}
	else {

		&print_diff_hsh( $DIFF_HSH, "comparison_differences" );
		&print_same_cluster_hsh( $SAME_HSH, "comparison_similarities" );
		&print_perfect_cluster_hsh( $SAME_HSH, "comparison_perfect_clusters" );
	}
}

sub print_perfect_cluster_hsh {
	my ( $results_hsh, $name ) = @_;

	open( FILE, ">$output_dir/$name.txt" )
		|| die "can't open $output_dir/$name: $!";

	foreach my $locus ( keys %$results_hsh ) {
		if ( $locus ne 'count' ) {
			if ( exists $results_hsh->{$locus}->{'perfect'} ) {
				if ( $locus ne 'count' ) {

					print FILE "\nLoci: $locus\n";

					foreach my $method ( keys %{ $results_hsh->{$locus} } ) {
						if ( $method ne 'perfect' && $method ne 'fs' ) {
							my $cluster_id = $results_hsh->{$locus}->{$method};

							print FILE "\t$method: $cluster_id\n";
						}
					}
				}
			}
		}
	}

	close FILE;

}

sub print_same_cluster_hsh {
	my ( $results_hsh, $name ) = @_;

	open( FILE, ">$output_dir/$name.txt" )
		|| die "can't open $output_dir/$name: $!";

	foreach my $locus ( keys %$results_hsh ) {
		if ( $locus ne 'count' ) {

			print FILE "\nLoci: $locus\n";

			foreach my $method ( keys %{ $results_hsh->{$locus} } ) {
				if ( $method !~ /(perfect|fs)/ ) {
					my $cluster_id = $results_hsh->{$locus}->{$method};

					print FILE "\t$method: $cluster_id\n";
				}
			}

			print FILE "\n\tFrameshift evidence\n" if ( exists $results_hsh->{$locus}->{'fs'} );
		}
	}

	close FILE;
}

sub print_cluster_hsh {
	my ( $results_hsh, $name ) = @_;

	open( FILE, ">$output_dir/$name" ) || die "can't open $output_dir/$name: $!";

	foreach my $locus ( keys %$results_hsh ) {

		print FILE "\nLocus: $locus\n";

		foreach my $method ( keys %{ $results_hsh->{$locus} } ) {

			my $cluster_id = $results_hsh->{$locus}->{$method};

			if ($cluster_id) {
				print FILE "\t$method: $cluster_id\n" if ($cluster_id);

				foreach my $db (@genomes) {
					my $key_lookup = $db . "_loci";

					if ( exists $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup} ) {
						print FILE
							"\t\t$db: $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup}\n";
					}
				}
			}
		}
	}

	close FILE;
}

sub print_diff_hsh {
	my ( $results_hsh, $name ) = @_;

	open( FILE, ">$output_dir/$name.txt" ) || die "can't open $output_dir/$name: $!";

	my $count = 1;

	foreach my $ids ( sort { $b cmp $a } keys %$results_hsh ) {

		my ( $hmms_hsh, $info_hsh, $frame_hsh, $prk_hsh );

		print FILE "\nCluster Difference " . $count++ . "\n";

		foreach my $method ( sort keys %{ $results_hsh->{$ids} } ) {

			my $cluster_id = $results_hsh->{$ids}->{$method};
			$cluster_id =~ s/\|//;
			print FILE "\t$method: Cluster $cluster_id\n" if $cluster_id;

			foreach my $db ( sort @genomes ) {

				my $key_lookup = $db . "_loci";

				if ( exists $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup} ) {
					my @loci
						= split( '\|', $METHOD_CLUSTERS_ID{$method}->{$cluster_id}->{$key_lookup} );

					my @feat_name;

					foreach my $element (@loci) {
						my $feat_names;

						if ( exists $LOCI_DB_LOOKUP->{$element}->{'feat_name'} ) {
							$feat_names = "$element($LOCI_DB_LOOKUP->{$element}->{'feat_name'})";
						}
						else {
							$feat_names = $element;
						}

						push( @feat_name, $feat_names );

						if ( exists $LOCI_DB_LOOKUP->{$element}->{'hmm'} ) {
							my @hmms = keys %{ $LOCI_DB_LOOKUP->{$element}->{'hmm'} };

							foreach my $hmm (@hmms) {
								$hmms_hsh->{$hmm} = 1;
							}
						}

						if ( exists $LOCI_DB_LOOKUP->{$element}->{'fs'} ) {
							my @frame_shifts = keys %{ $LOCI_DB_LOOKUP->{$element}->{'fs'} };

							foreach my $frame (@frame_shifts) {
								$frame_hsh->{$element}->{$frame} = 1;
							}

						}

						if ( exists $LOCI_DB_LOOKUP->{$element}->{'prk'} ) {
							my @prk_shifts = keys %{ $LOCI_DB_LOOKUP->{$element}->{'prk'} };

							foreach my $prk (@prk_shifts) {
								$prk_hsh->{$prk} = 1;
							}

						}

					}

					my $feat_names = join( '|', sort @feat_name );
					print FILE "\t\t$db: $feat_names \n";
				}

			}

		}

		if ( scalar keys %$hmms_hsh > 0 ) {
			foreach my $hmm ( keys %$hmms_hsh ) {
				my @loci_match = keys %{ $EV_LOOKUP->{'hmm'}->{$hmm}->{'feat_name'} };
				my $loci_matches = join( '|', @loci_match );

				print FILE "\n\tEquivalog HMM: $hmm ($EV_LOOKUP->{'hmm'}->{$hmm}->{'com_name'}) \n";
				print FILE "\tLoci Matches: $loci_matches";

				print FILE "\n";
			}
		}

		if ( scalar keys %$prk_hsh > 0 ) {
			foreach my $prk ( keys %$prk_hsh ) {
				my @loci_match = keys %{ $EV_LOOKUP->{'prk'}->{$prk}->{'feat_name'} };
				my $loci_matches = join( '|', @loci_match );

				print FILE "\n\tPRK Clusters: $prk ($EV_LOOKUP->{'prk'}->{$prk}->{'com_name'}) \n";
				print FILE "\tLoci Matches: $loci_matches";

				print FILE "\n";
			}
		}

		if ( scalar keys %$frame_hsh > 0 ) {
			print FILE "\n\tFRAMESHIFT: \n";

			foreach my $loci ( sort keys %$frame_hsh ) {
				my @types = keys %{ $frame_hsh->{$loci} };
				foreach my $type (@types) {
					print FILE "\t\t$loci ($type) \n";
				}
			}

			print FILE "\n";

		}

	}

	close FILE;
}

sub sort_hashkeys_by_value {
	my ($hashref) = @_;
	return sort { $hashref->{$a} <=> $hashref->{$b} } keys %$hashref;
}

sub sort_multihashkeys_by_value {
	my ( $hashref, $key1 ) = @_;
	return sort { $hashref->{$a}->{$key1} <=> $hashref->{$b}->{$key1} } keys %$hashref;
}

sub sort_pipes {
	my ( $hashref, $key1, $key2 ) = @_;

	my @x = split( /\|/, $hashref->{$key1}->{$key2} );

	@x = sort(@x);

	my $new_loci_string;

	for ( my $i = 0; $i < @x; $i++ ) {

		$new_loci_string .= $x[$i] . "|";
	}

	if ( $new_loci_string =~ /\|\Z/ ) { chop($new_loci_string); }

	$hashref->{$key1}->{$key2} = $new_loci_string;

	return ($hashref);
}

sub sort_pipes_2 {
	my ($string) = @_;

	my @x = split( /\|/, $string );

	@x = sort(@x);
	@x = sort { $b cmp $a } (@x);

	my $new_loci_string;

	for ( my $i = 0; $i < @x; $i++ ) {

		$new_loci_string .= $x[$i] . "|";

	}

	if ( $new_loci_string =~ /\|\Z/ ) { chop($new_loci_string); }

	return ($new_loci_string);
}

sub check_options {
	my $error_message;

	if ( $opts{help} ) {
		&_pod(2);
		exit(0);
	}

	if ( $opts{method_files} ) {
		@FILES = read_file( $opts{method_files} );

		if ( $opts{pan_db} || @ARGV >= 1 ) {
			$error_message
				.= "Can only supply one of the following: --method_files || --pan_db || pass in the individual result files as arguments\n";
		}
	}
	elsif ( scalar @ARGV >= 1 ) {
		@FILES = @ARGV;

		if ( $opts{pan_db} ) {
			$error_message
				.= "Can only supply one of the following: --method_files || --pan_db || pass in the individual result files as arguments\n";
		}
	}
	elsif ( !$opts{method_files} && scalar @ARGV < 1 && !$opts{pan_db} ) {
		$error_message
			.= "No cluster methods found. Use either --method_list || --pan_db || argument list\n";
	}

	if ( !$opts{db_list} && !$opts{fasta_files} ) {
		$error_message .= "No sequence dat found. Use either --db_list || --fasta_file.\n";
	}

	if ( $opts{'db_list'} && $opts{'fasta_files'} ) {
		$error_message .= "--db_list and --fasta_file can not both be used.\n";
	}

	if ( $opts{user} && $opts{password} ) {
		$user     = $opts{'user'};
		$password = $opts{'password'};
	}
	else {
		$user     = $DEFAULT_USER;
		$password = $DEFAULT_PASS;
	}

	if ( $opts{output_dir} ) {
		$output_dir = $opts{output_dir};
	}
	else {
		$output_dir = $DEFAULT_OUTPUT_DIR;
	}

	die $error_message if $error_message;
}

sub _pod {
	my $level = shift;
	pod2usage( { -exitval => 0, -verbose => $level, -output => \*STDERR } );
}

sub do_sql {
	my ( $dbproc, $query, $delimeter ) = @_;
	my ( $statementHandle, @x,      @results );
	my ( $i,               $result, @row );

	$delimeter = "\t";

	$statementHandle = $dbproc->prepare($query);
	if ( !defined $statementHandle ) {
		die "Cannot prepare statement: $DBI::errstr\n";
	}
	$statementHandle->execute() || die "failed query: $query\n";
	if ( $statementHandle->{syb_more_results} ne "" ) {
		while ( @row = $statementHandle->fetchrow() ) {
			push( @results, join( $delimeter, @row ) );
		}
	}

	#release the statement handle resources
	$statementHandle->finish;
	return (@results);
}

sub RunMod {
	my ( $dbproc, $query ) = @_;
	my ( $statementHandle, $result );

	$statementHandle = $dbproc->prepare($query);
	if ( !defined $statementHandle ) {
		die "Cannot prepare statement: $DBI::errstr\n$query\n";
	}
	$statementHandle->execute() || print "failed query: $query\n";
	$statementHandle->finish;
	return ($result);
}

sub ConnectToDb {
	my ( $SERVER, $DBTYPE, $user, $passwd, $db ) = @_;
	my ($dbproc);

	$dbproc = ( DBI->connect( "dbi:$DBTYPE:server=$SERVER", $user, $passwd ) );

	if ( !defined $dbproc ) {
		die "Cannot connect to Sybase server: $DBI::errstr\n";
	}
	$dbproc->do("use $db");
	return ($dbproc);
}
