#!/usr/local/bin/perl -w

use warnings;
use strict;
$|++;

=head1 NAME

pangenome_statistics.pl - calculates general statistics and files

=head1 SYNOPSIS

    USAGE: pangenome_statistics.pl   -a <data hsh file>
                                     -m <method result file>
    	                             -l <genome list>
                                     -s <combined pep fasta file>
                                     -n <combined seq fasta file>
                                     [-f <frameshift file>
                                      --fusion <fusion file>
                                      -c <centroid>
                                      -r <genome>
                                      -i <role ids>
                                      -t <com_name to limit by>
                                      --role_lookup
                                      -o <output directory>
                                     ]
         
=head1 OPTIONS

B<--data_file,-a>       :   combined.att_file data hash file

B<--method,-m>          :   cluster method result file

B<--genomes_list,-l>         :   Genome list file

B<--combined_fasta, s>  :   Path to combined pep fasta file (combined.fasta)

B<--combined_seq, n>    :   Path to combined seq fasta file (combined.seq)

B<--centorid,-c>        :   Centroid fasta file  <optional>

B<--reference,-r>       :   Reference genome to limit results by <optional>

B<--role_ids,-i>        :   Role ids to limit results by <optional>

B<--terms, -t>          :   Com_name terms to limit results by <optional>

B<--frameshift,-f>      :   Path to a frameshift file <optional>

B<--fusion>             :   Path to fusion/fragment file <optional>

B<--role_lookup>        :   File of role_ids and their associated mainrole/subrole categories

B<--help,-h>            :   Display this help message.

=head1  DESCRIPTION

This program will give overview counts of the clusters generated from a specificed clustering method. This
program also generates output files listing the core, shared and unique clusters.

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Cwd;
use DBI;
use Cwd 'abs_path';
use feature qw(switch);
use File::Basename;
use File::Slurp;
use File::Temp;
use FindBin;
use MLDBM qw( DB_File Data::Dumper );
use Fcntl;
use Data::Dumper;
use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
use lib "$FindBin::Bin/../lib";
use File::Path;
use TIGR::FASTArecord;
use TIGR::FASTAreader;
use TIGR::FASTAwriter;
use TIGR::Foundation;
use File::Copy;

my $server = '';
my $method = '';
my %opts;

GetOptions( \%opts, 
	    'data_file|a=s',
	    'centroids|c=s',
	    'method_result|m=s',
	    'genomes_list|l=s',
	    'reference|r=s',
	    'role_ids|i=s',
	    'terms|t=s',
	    'frameshift|f=s',
	    'fusion=s',
	    'combined_fasta|s=s',
	    'combined_seq|n=s',
	    'output|o=s',
	    'hmm_file=s',
	    'role_lookup=s',
	    'help|h') || die "Error getting options! $!";
                    
pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my($ROLES_OF_INTEREST,@NAMES_OF_INTEREST,$REFERENCES,$OUTPUT,$THRESHOLDS); # set in check_params

&check_params;
#opens att file tied hash
tie my %ANNOT_HSH, 'MLDBM', $opts{data_file} or die "Can't open MLDBM file: $!\n";

#global variables
my ($TYPE_COUNTS,$COM_NAMES);
my ($HISTOGRAM,$FRAME_COUNT);
my ($INITIAL_COUNTS);
my ($IGNORED_FRAME_SING);

my ($DATABASES,$db_array) = &process_db_file($opts{genomes_list});
my @dbs = @$db_array;

#Parse HMM Info.
my @HMMs = read_file($opts{hmm_file}) if $opts{hmm_file};
my $HMM_MATCHES;

foreach my $hmm(@HMMs){
    $hmm =~ s/\s+$//;
    $HMM_MATCHES->{$hmm}->{core} = 1;
}

my($sing_c,$core_c,$shared_c,$all_c) = &parse_cluster_results($opts{method_result},scalar @dbs);
my ($SING_FRAME) = &parse_frameshift($opts{frameshift}) if $opts{frameshift};
my ($FUSION) = &parse_fusion($opts{fusion}) if $opts{fusion};

my $CENTROIDS = &parse_centroids($opts{centroids}) if $opts{centroids};

my ($ROLE_ID_LOOKUP,$SUBROLE_LOOKUP) = &parse_role_id_lookup($opts{role_lookup}) if $opts{role_lookup};
my ($CLUSTER_ATT,$ROLE_IDS) = &find_clusters_to_print($all_c);

&print_clusters_of_interest($core_c,"core");
&print_clusters_of_interest($sing_c,"singletons");
&print_clusters_of_interest($sing_c,"singletons_no_frame");
&print_clusters_of_interest($shared_c,"shared");
&print_clusters_of_interest($shared_c,"shared_no_ref") if $opts{reference};
&print_clusters_of_interest($shared_c,"core_no_ref") if $opts{reference};

&make_singleton_fasta_files($sing_c) if($opts{combined_fasta} || $opts{combined_seq});
&make_presence_absence_lists($all_c);

&print_overview;

untie %ANNOT_HSH,
exit(0);
############################################
sub print_role_id{
    my ($hsh,$groups) = @_;

    foreach my $thresh (@$THRESHOLDS){
	open(my $fh,">","$OUTPUT/role_id_counts_" . $thresh . ".txt");

	#print headers
	print $fh "role_id\tmain role\tsub role\t";
	map{print $fh "$_\t"} sort(@{$groups});
	print $fh "\n";
	
	foreach my $mainrole(sort{$a cmp $b} keys $SUBROLE_LOOKUP){
	    
	    my @subroles = keys %{$SUBROLE_LOOKUP->{$mainrole}};
	    
	    foreach my $role(sort @subroles){
		
		my $id = $SUBROLE_LOOKUP->{$mainrole}->{$role};
		
		print $fh "$id\t$mainrole\t$role\t";
		foreach my $label(sort @$groups){
		    
		    if(exists $hsh->{$mainrole}->{$label}->{$thresh}->{$role}){
			print $fh $hsh->{$mainrole}->{$label}->{$thresh}->{$role}->{count} . "\t";
		    }else{
			print $fh "0\t";
		    }
		}
		
		print $fh "\n";
	    }
	    
	    #Print Mainrole total count line
	    print $fh "MRSUM\t$mainrole\t\t";
	    foreach my $label(sort @$groups){
		
		if(exists $hsh->{$mainrole}->{$label}->{$thresh}->{total}->{count}){
		    print $fh "$hsh->{$mainrole}->{$label}->{$thresh}->{total}->{count}" . "\t";
		}else{
		    print $fh "0\t";
		}
	    }
	    
	    print $fh "\n";
	}
    }
}
sub find_clusters_to_print{
    #Based on the parameters passed to this program
    #This sub determines whether a cluster should be printed in the various
    #stat scripts. It also determines cluster attributes and stores those
    #in the same hash.

    my $all_clusters = shift;
    my ($print_hsh,$no_ref,$role_hsh);
    my $cluster_att;

    foreach my $cluster(sort {$a<=>$b} keys %$all_clusters){
	my @members;

	foreach my $db(keys %{$all_clusters->{$cluster}}){
	    push(@members,$all_clusters->{$cluster}->{$db});
	}

	my ($role_ids,$com_name,$includes_reference,$frameshift,$fusion_types,$hmms) = &find_cluster_attributes(\@members, $REFERENCES);

	#Set com_name to be centroid unless a reference was given and it was found in this cluster
	if($com_name && !($opts{reference})){

	    if($CENTROIDS){

		my $centroid = $CENTROIDS->{$cluster}->{locus};
		$com_name = $CENTROIDS->{$cluster}->{com_name} . "\t" . $ANNOT_HSH{$centroid}->{db};
	    }
	}

	#Set global com_name hash
	$COM_NAMES->{$cluster} = $com_name;

	my($match_role,$match_term);
	my @role_ids = @$role_ids;

     	if($opts{role_ids}){
	    
	    foreach (@role_ids){
		$match_role = 1 if($_ =~ /($ROLES_OF_INTEREST)/);
	    }
	}

      	if($opts{terms})
	{
	    foreach my $term(@NAMES_OF_INTEREST){
		$match_term =1 if $com_name =~ /$term/;
	    }
	    
	}

	if($opts{role_ids}){
	    
	    if($opts{terms}){
		
		if($opts{reference}){
		    
		    if ($match_role || $match_term){

			if($includes_reference){

			    $cluster_att->{$cluster}->{PRINT} = 1;

			}else{

			    $cluster_att->{$cluster}->{PRINT} = 0;
			}
		    }
		    
		}else{
		    
		    $cluster_att->{$cluster}->{PRINT} = 1 if ($match_role || $match_term);
		    		    
		}
		
	    }elsif($opts{reference}){
		
		if ($match_role){

		    if($includes_reference){

			 $cluster_att->{$cluster}->{PRINT} = 1;

		    }else{

			$cluster_att->{$cluster}->{PRINT} = 0;
		    }
		}
		
	    }else{
		
		$cluster_att->{$cluster}->{PRINT} = 1 if($match_role);
		
	    }
	    
	}elsif($opts{terms}){
	    
	    if($opts{reference}){
		
		if ($match_term){

		    if($includes_reference){

			 $cluster_att->{$cluster}->{PRINT} = 1;

		    }else{

			$cluster_att->{$cluster}->{PRINT} = 0;

		    }
		}
		
	    }else{
		
		$cluster_att->{$cluster}->{PRINT} = 1 if ($match_term);
		
	    }
	    
	}elsif($opts{reference}){

	    if($includes_reference){

		$cluster_att->{$cluster}->{PRINT} = 1;

	    }else{

		$cluster_att->{$cluster}->{PRINT} = 0;
	    }
	    
	}else{

	    $cluster_att->{$cluster}->{PRINT} = 1;

	}

	#Add role_id_info and frameshift
	if(scalar @role_ids > 0){
	    
	    $cluster_att->{$cluster}->{role_id} =  join(",",sort @role_ids);
	   
	    map{$role_hsh->{$_}++} @role_ids if $cluster_att->{$cluster}->{PRINT} == 1;
	    
	}else{
	    $role_hsh->{0}++;
	    $cluster_att->{$cluster}->{role_id} = 0;
	}

	$cluster_att->{$cluster}->{frameshift} = "FS" if $frameshift;
	$cluster_att->{$cluster}->{fusion} = $fusion_types if $fusion_types;
	$cluster_att->{$cluster}->{hmms} = $hmms;
	    
    }

    return ($cluster_att,$role_hsh);
}
sub make_singleton_fasta_files{
    my $clusters = shift;
    my $tf_object = new TIGR::Foundation;
    my @errors;

    my $fpr_obj = new TIGR::FASTAreader ($tf_object,\@errors, $opts{combined_fasta}) if ($opts{combined_fasta});
    $INITIAL_COUNTS->{members} = $fpr_obj->{num_records} if ($opts{combined_fasta});
    my $fsr_obj = new TIGR::FASTAreader ($tf_object,\@errors, $opts{combined_seq}) if ($opts{combined_seq});
    $INITIAL_COUNTS->{members} = $fsr_obj->{num_records} if ($opts{combined_seq});

    my $fpw = new TIGR::FASTAwriter($tf_object, \@errors) if ($opts{combined_fasta});
    my $fsw = new TIGR::FASTAwriter($tf_object, \@errors) if ($opts{combined_seq});
    my $fnfpw = new TIGR::FASTAwriter($tf_object, \@errors) if ($opts{combined_fasta});
    my $fnfsw = new TIGR::FASTAwriter($tf_object, \@errors) if ($opts{combined_seq});

    $fpw->open("$OUTPUT/singletons.pep") if ($opts{combined_fasta});
    $fsw->open("$OUTPUT/singletons.seq") if ($opts{combined_seq});
    $fnfpw->open("$OUTPUT/singletons_no_frame.pep") if ($opts{combined_fasta});
    $fnfsw->open("$OUTPUT/singletons_no_frame.seq") if ($opts{combined_seq});

    foreach my $cluster_id(sort {$clusters->{$a} cmp $clusters->{$b}} keys %$clusters){
	
	if(exists $CLUSTER_ATT->{$cluster_id}->{PRINT} == 1){
	    
	    my $presult = $fpr_obj->getRecordByIdentifier($clusters->{$cluster_id}) if ($opts{combined_fasta});
	    my $sresult = $fsr_obj->getRecordByIdentifier($clusters->{$cluster_id}) if ($opts{combined_seq});
	    
	    my $com_name = $ANNOT_HSH{$clusters->{$cluster_id}}->{com_name};
	    $com_name =~ s/\// /g;
	    
	    my $header = "$clusters->{$cluster_id}|$cluster_id|$ANNOT_HSH{$clusters->{$cluster_id}}->{db}|$com_name";
	    my $p_rec = new TIGR::FASTArecord($header,$presult->{data_rec}) if ($opts{combined_fasta});
	    my $s_rec = new TIGR::FASTArecord($header,$sresult->{data_rec}) if ($opts{combined_seq});
	    
	    $fpw->write($p_rec ) if ($opts{combined_fasta});
	    $fsw->write($s_rec ) if ($opts{combined_seq});
	    
	    unless(exists $IGNORED_FRAME_SING->{$cluster_id}){
		
		$fnfpw->write($p_rec) if ($opts{combined_fasta});
		$fnfsw->write($s_rec) if ($opts{combined_seq});
		
	    }
	    
	}
	
    }

}
sub print_overview{
    open(my $sfh, ">", "$OUTPUT/overview_stats.txt");
    select $sfh;

    print "Limited by:\n";
    print "$opts{reference}\n" if $opts{reference};
    print "$opts{role_ids}\n" if $opts{role_ids};
    print "$opts{terms}\n" if $opts{terms};
    print "\n";
    
    $TYPE_COUNTS->{singletons}->{cluster} = 0 unless(exists $TYPE_COUNTS->{singletons});
    $TYPE_COUNTS->{singletons}->{loci} = 0 unless(exists $TYPE_COUNTS->{singletons}->{loci});
    
    print "Before limits \n";
    print "Total Clusters: " . ($INITIAL_COUNTS->{cluster}) . "\n";
    print "Total Loci: " . ($INITIAL_COUNTS->{members}) . "\n" if defined $INITIAL_COUNTS->{members};
    print "Total Fragmented Genes: " . ($INITIAL_COUNTS->{fragments}) . "\n" if exists $INITIAL_COUNTS->{fragments};
    print "Total Genomes: " . scalar @dbs . "\n";
    print "\n";

    print "After limits\n";
    print "Total Clusters: " . ($TYPE_COUNTS->{shared}->{cluster} + $TYPE_COUNTS->{singletons}->{cluster}) . "\n";
    print "Total Loci (Ignores Fragments): " . ($TYPE_COUNTS->{shared}->{loci} + $TYPE_COUNTS->{singletons}->{loci}) . "\n";
    print "Total Genomes: " . scalar @dbs . "\n";
    print "\n";

    print "Type\tClusters\tLoci\n";
    my %print_rep = ("shared_no_ref" => "Shared clusters that do no include reference",
		     "core_no_ref" => "Clusters that include all genomes except reference",
		     "singletons" => "Singleton Clusters",
		     "core" => "Core Clusters",
		     "shared" => "Shared Clusters",
	             "singletons_no_frame" => "Singletons that are not frameshifted");
    

    foreach my $type (sort keys %$TYPE_COUNTS){
	my $display = $print_rep{$type};
	print  $display . ": " . $TYPE_COUNTS->{$type}->{cluster} . "\t" . $TYPE_COUNTS->{$type}->{loci} . "\n";
    }
    
    print "Clusters with Fragments: $FRAME_COUNT\n" if $FRAME_COUNT;
    
    #Print HMM info
    if($opts{hmm_file}){
	my($hmm_total,$hmm_matched) = (0,0);
	my $hmm_missed;

	$hmm_total = scalar keys %$HMM_MATCHES;

	foreach my $key(keys %$HMM_MATCHES){
	    my @values = keys %{$HMM_MATCHES->{$key}};
	    
	    if((scalar @values) == 2){
		$hmm_matched++;
	    }else{
		$hmm_missed .= "$key\n";
	    }
	}

	my $percentage = $hmm_matched/$hmm_total * 100;

	print "\nHMM Information\n";
	print "Total HMMs of interest: $hmm_total\n";
	print "HMMs of interest found in clusters: $hmm_matched (" . sprintf("%.2f", $percentage) . "%)\n";
	print "Missing HMMs:\n";
	print $hmm_missed;
	print "\n";
    }

    print "\nCluster Size Breakdown\n";
    
    my $combined_cluster;
    
    foreach my $key(keys %$HISTOGRAM){
	foreach my $num(keys %{$HISTOGRAM->{$key}}){
	    $combined_cluster->{$num} = $HISTOGRAM->{$key}->{$num};
	}
    }
    
    foreach my $count(sort {$a<=>$b} keys %$combined_cluster){
	print $count . "\t" . $combined_cluster->{$count} . "\n";
    }

    if($opts{role_lookup}){
	print "\nRole ID Breakdown\n";
	print "role_id\tMain Role\tSub Role\tCluster Count\n";
	
	foreach my $role_id(sort {$ROLE_IDS->{$b} <=> $ROLE_IDS->{$a}} keys %$ROLE_IDS){
	    if($role_id == 0){
		print $role_id . "\t" . "NONE ASSIGNED" . "\t" . "NONE ASSIGNED" . "\t" . $ROLE_IDS->{$role_id} . "\n";
	    }else{
		print $role_id . "\t" . $ROLE_ID_LOOKUP->{$role_id}->{main} . "\t" . $ROLE_ID_LOOKUP->{$role_id}->{subrole} . "\t" . $ROLE_IDS->{$role_id} . "\n";
	    }
	}
    }

    select STDOUT;
}
sub make_presence_absence_lists{
    my $complete_clusters = shift;

    open(my $ofh,">", "$OUTPUT/all_clusters_member_presence.txt");
    open(my $afh,">", "$OUTPUT/all_clusters_members.txt");
    open(my $hfh,">", "$OUTPUT/hmm_clusters.txt") if($opts{hmm_file});
    
    print $ofh "cluster_id\tprotein name\tcentroid genome\tcentroid locus\trole_ids\tframeshift\t";
    print $afh "cluster_id\tprotein name\tcentroid genome\tcentroid locus\trole_ids\tframeshift\t";
    print $hfh "cluster_id\tprotein name\tcentroid genome\tcentroid locus\trole_ids\tframeshift\tHMMs\t" if($opts{hmm_file});

    foreach my $db(@dbs){
	print $ofh "$db\t";
	print $afh "$db\t";
	print $hfh "$db\t" if($opts{hmm_file});
    }

    my ($ps_h,$ps_a,$ps_o);
    my $hmm_count; 

    foreach my $cluster(sort {$a<=>$b} keys %$complete_clusters){
	my $print_hmm = 1 if (defined $CLUSTER_ATT->{$cluster}->{hmms} && $opts{hmm_file});
	
	if($print_hmm){
	    foreach my $hmm(@{$CLUSTER_ATT->{$cluster}->{hmms}}){
		$HMM_MATCHES->{$hmm}->{count}++ if (exists $HMM_MATCHES->{$hmm});
	    }
	}

	$ps_o = "\n$cluster\t";
	$ps_a = "\n$cluster\t";
	$ps_h = "\n$cluster\t" if $print_hmm;

	#Print com names/protein names
	my($com_name,$com_name_db) = split(/\t/,$COM_NAMES->{$cluster});
	$ps_o .= $com_name . "\t";
	$ps_a .= $com_name . "\t";
	$ps_h .= $com_name . "\t" if $print_hmm;

	if($CENTROIDS){
	    #Finds centroid information to print
	    my $centroid_locus = $CENTROIDS->{$cluster}->{locus};
	    my $centroid_db = $ANNOT_HSH{$centroid_locus}->{db};
	    
	    $ps_o .= $centroid_db  . "\t" . $centroid_locus . "\t";
	    $ps_a .= $centroid_db  . "\t" . $centroid_locus . "\t";
	    $ps_h .= $centroid_db  . "\t" . $centroid_locus . "\t" if $print_hmm;
	}else{
	    $ps_o .= "\t\t";
	    $ps_a .= "\t\t";
	    $ps_h .= "\t\t" if $print_hmm;
	}	
		
	#Print role ids and HMM information
	my $role_ids = $CLUSTER_ATT->{$cluster}->{role_id};
	my $frameshift = $CLUSTER_ATT->{$cluster}->{frameshift} if exists $CLUSTER_ATT->{$cluster}->{frameshift};
	
	if($role_ids && $role_ids ne 'NONE'){

	    $ps_o .= $role_ids . "\t";
	    $ps_a .= $role_ids . "\t";
	    $ps_h .= $role_ids . "\t" if $print_hmm;

	}else{

	    $ps_o .= "\t";
	    $ps_a .= "\t";
	    $ps_h .= "\t" if $print_hmm;
	}
	
	if($frameshift){

	    $ps_o .= $frameshift . "\t";
	    $ps_a .= $frameshift . "\t";
	    $ps_h .= $frameshift . "\t" if $print_hmm;

	}else{

	   $ps_o .= "\t";
	   $ps_a .= "\t"; 
	   $ps_h .= "\t" if $print_hmm;
	}
	

	my $print_hmm_string = join("|",@{$CLUSTER_ATT->{$cluster}->{hmms}});
	$ps_h .= $print_hmm_string . "\t" if $print_hmm;
       
	#print the clusters presence/absence
	foreach my $db(@dbs){

	    if(exists $complete_clusters->{$cluster}->{$db}){

		$ps_o .= "1\t";
		$ps_a .= "$complete_clusters->{$cluster}->{$db}\t";
		$ps_h .= "$complete_clusters->{$cluster}->{$db}\t" if $print_hmm;

	    }else{

		$ps_o .= "0\t";
		$ps_a .= "\t";
		$ps_h .= "\t" if $print_hmm;

	    }	    
	}
	
	if(exists $CLUSTER_ATT->{$cluster}->{PRINT}){
	    print $ofh $ps_o if $CLUSTER_ATT->{$cluster}->{PRINT} == 1;
	    print $afh $ps_a if $CLUSTER_ATT->{$cluster}->{PRINT} == 1;
	}

	print $hfh $ps_h if $ps_h;
    }

    close $afh;
    close $ofh;
    close $hfh if $hfh;
}

sub parse_centroids{
    my $file = shift;

    my @def_lines = `grep \">\" $file`;

    my $hsh;
    foreach my $line(@def_lines){
	$line =~ s/\s+$//;
	my @values = split(/\s/,$line);
	my $identifier = shift(@values);
	my $locus = shift(@values);

	my($prefix,$cluster_id) = split(/_/,$identifier);
	$hsh->{$cluster_id}->{locus} = $locus;
	$hsh->{$cluster_id}->{com_name} = join(" ",@values);

    }

    return $hsh;
}
sub parse_frameshift{
    my $file = shift;
    my $hsh;

    open(FILE,"<",$file);

    while(<FILE>){
	my $line = $_;
	chomp $line;
	unless ($line =~ /^>/){
	    $INITIAL_COUNTS->{fragments}++;
	    my @values = split(/\t/,$line);
	    map{$hsh->{$_} = 1} @values;
	}
    }

    close FILE;

    return ($hsh);
}
sub parse_fusion{
    my $file = shift;
    my $hsh;

    open(my $fh,"<",$file);
    while(my $line = <$fh>){
	$line =~ s/\s+$//;
	my @values = split(/\t/,$line);

	if($values[1] =~ /Fragment in/){
	    $hsh->{$values[0]} = "FG-In";
	}elsif($values[1] =~ /Fragment out/){
	    $hsh->{$values[0]} = "FG-Out";
	}else{
	    $hsh->{$values[0]} = "FN";
	}
    }

    close $fh;
    return $hsh;
}
sub find_cluster_attributes{
    my ($members, $reference) = @_;
    my ($includes_reference,$includes_frameshift);
    my ($com_name,$db_member,$role_id_found_hsh);
    my ($fusion_hsh,$hmm_hsh);
    
    # Loop through each member of a cluster to determine
    # if that cluster has certain attributes, such 
    # as frameshifts or members of a reference
    foreach my $member(@$members){

	my $db = $ANNOT_HSH{$member}->{db};
	my $hmms = $ANNOT_HSH{$member}->{hmms} if (defined $ANNOT_HSH{$member}->{hmms});

	if($hmms){
	    my @hmm = split(/\|/,$hmms);
	    map{$hmm_hsh->{$_} = 1} @hmm;
	}

	$includes_frameshift = 1 if exists $SING_FRAME->{$member}; 

	if(exists $FUSION->{$member}){
	    my $type = $FUSION->{$member};
	    $fusion_hsh->{$type} = 1;
	}

	# If a reference was given and that db was found
	# set reference flag and set the com_name for cluster to be from that member
	if($reference && ($db =~ /($reference)/)){

	    $includes_reference = 1;
	    $com_name = $ANNOT_HSH{$member}->{com_name} . "\t" . "$db";

	}else{
	    
	    #Set com_name to current member if no reference has been found
	    $com_name = $ANNOT_HSH{$member}->{com_name} . "\t" . "$db" unless($includes_reference);

	}

	#Find role ids associated with current member
	if(defined $ANNOT_HSH{$member}->{role_id}){

	    my $role_values = $ANNOT_HSH{$member}->{role_id};
	    my @role_ids = split(/\|/,$role_values);

	    foreach my $id(@role_ids){
		$id =~ s/\s+$//;
		$id =~ s/^\s+//;
		push(@{$role_id_found_hsh->{$id}},$member);
	    }
	}
    }
    
    my $fusion_types = join(",", keys %$fusion_hsh);
    my @role_ids = keys %$role_id_found_hsh;
    my @hmm_ids = keys %$hmm_hsh;

    return (\@role_ids,$com_name,$includes_reference,$includes_frameshift,$fusion_types,\@hmm_ids);    
}
sub get_print_string{
    my($id,$cluster_info,$db_members,$type,$member_count) = @_;

    my $role_ids = $CLUSTER_ATT->{$id}->{role_id};
    my $com_name = $COM_NAMES->{$id};
    my $frameshift = $CLUSTER_ATT->{$id}->{frameshift} if exists $CLUSTER_ATT->{$id}->{frameshift};
    my $fusion = $CLUSTER_ATT->{$id}->{fusion} if exists $CLUSTER_ATT->{$id}->{fusion};

    my $string = $id . "\t";
    $string .= $member_count . "\t";

    #Parse com_name
    my($name,$com_name_db) = split(/\t/,$com_name);
    $string .= $name . "\t";

    #Find centroid info
    unless($type eq 'singletons'){
	my $centroid_locus = $CENTROIDS->{$id}->{locus};
	my $centroid_db = $ANNOT_HSH{$centroid_locus}->{db};
	$string .= $centroid_db . "\t" . $centroid_locus . "\t";
    }

    #Print role ids and frameshift info
    $string .= ($role_ids ne "NONE") ?  $role_ids . "\t" : "\t";
    
    if($frameshift){
	$string .= ($fusion) ? "$frameshift,$fusion\t" : $frameshift . "\t";
    }elsif($fusion){
	$string .= $fusion . "\t";
    }else{
	$string .= "\t";
    }
    
    if($type eq 'singletons'){

	foreach my $db (keys %$db_members){

	    my @members = keys %{$db_members->{$db}};
	    $string .= "$members[0]\t$db";

	}

    }else{

	foreach my $db(@dbs){
	    
	    if(exists $db_members->{$db}){

		my @values = keys %{$db_members->{$db}};

		my $values_string = join("|",sort @values);
		$string .= "$values_string\t";
	       
		#Want to keep track of number of loci per genome, but to avoid duplicate
		#counts focus on shared and singletons
		$DATABASES->{$db} += scalar @values if ($type =~ /\b(singletons|shared)\b/);

	    }else{

		$string .= "\t";

	    }
	}
    }

    return $string;
}
sub print_clusters_of_interest{
    my ($clusters,$type) = @_;
    my ($cluster_count,$loci_count);

    #open OUT file
    my $file = $type . "_clusters.txt";   
    open(my $ofh, ">","$OUTPUT/$file");
    select $ofh;

    #Print file headers
    if($type =~ /singletons/){

	print"cluster id\tnum of members\tprotein name\trole_ids\tattributes\tlocus\tdb";

    }else{

	print"cluster id\tnum of members\tprotein name\tcentroid genome\tcentroid locus\trole_ids\tattributes\t";
	#map{print "$_\t"} sort keys %$DATABASES;
	map{print "$_\t"} @dbs;

    }

    #Print cluster information
    foreach my $cluster (sort {$a <=> $b} keys %$clusters){

	my $members = $clusters->{$cluster};
	my @members = split(",",$members);
	my $db_member = &find_db_members(\@members);
	my $frameshift = $CLUSTER_ATT->{$cluster}->{frameshift} if exists $CLUSTER_ATT->{$cluster}->{frameshift};

      	my $string = &get_print_string($cluster,$clusters->{$cluster},$db_member, $type,scalar @members);

	if($type =~ /no_ref/){

	    # Prints clusters that do not contain the passed in reference, these are in the CLUSTER_ATT
	    # hash but have PRINT set to zero
	    if($type eq 'shared_no_ref'){

		print "\n" . $string if($CLUSTER_ATT->{$cluster}->{PRINT} == 0)
		
	    }elsif($type eq 'core_no_ref'){

		my $cluster_db_count = scalar keys %$db_member;

		if($CLUSTER_ATT->{$cluster}->{PRINT} == 0 && (($cluster_db_count == (scalar @dbs - 1))) ){

		    print "\n" . $string;

		}
	    }

	}else{

	    if($CLUSTER_ATT->{$cluster}->{PRINT} == 1){
		
		if($type =~ /singletons_no_frame/){
		    
		    #Do not print clusters that have a frameshift in them
		    if($frameshift){
			
			$IGNORED_FRAME_SING->{$cluster} = 1;
			
		    }else{
			
			my $num_of_dbs = scalar keys %$db_member;
			print "\n" . $string;
			$TYPE_COUNTS->{$type}->{cluster}++;
			$TYPE_COUNTS->{$type}->{loci} += scalar @members;
		    
		    }
		    
		}else{
		    
		    my $num_of_dbs = scalar keys %$db_member;
		    $TYPE_COUNTS->{$type}->{cluster}++;
		    $TYPE_COUNTS->{$type}->{loci} += scalar @members;
		    $HISTOGRAM->{$type}->{$num_of_dbs}++ if ($type =~ /\b(shared|singletons)\b/);
		    $FRAME_COUNT++ if $frameshift;
		    
		    print "\n" . $string;
		    
		}
		
	    }
	}
    }
    
    select STDOUT;
}
sub find_db_members{
    my $members = shift;
    my $hsh;

    foreach my $member(@$members){

	my $db = $ANNOT_HSH{$member}->{db};
	if (!defined($db)){
	    #All gene identifiers in the table should have been defined in the attribute file
	    die ("$member was not defined in the attribute file $opts{data_file}\n");
	}
	$hsh->{$db}->{$member} = 1; 

    }

    return $hsh;

}
sub process_db_file{
    my $file = shift;
    my @lines = read_file($file);
    
    my $databases;
    my @db_array; 

    foreach my $line (@lines){
	$line =~ s/\s+$//;#remove trailing white space
	#There should only be a single genome identifier on each line
	if (length($line) == 0){
	    #skip blank lines
	    next;
	}

	$databases->{$line} = 0;
	push(@db_array,$line);
    }

    return ($databases,\@db_array);
}
sub find_potential_clusters{
    my ($cluster_hsh,$num_dbs) = @_;

    my($sing,$core,$shared);

    my $max_num = 0;

    foreach my $id (keys %$cluster_hsh){
	my $count = $cluster_hsh->{$id}->{member_count};
        $max_num = $count if $count > $max_num;

	$sing->{$id} = $cluster_hsh->{$id}->{members} if ($count == 1);
	$core->{$id} = $cluster_hsh->{$id}->{members} if ($count == $num_dbs);
	$shared->{$id} = $cluster_hsh->{$id}->{members} if ($count > 1);
    }

    return($sing,$core,$shared);
}
sub parse_cluster_results{
    my ($results_file,$num_dbs) = @_;

    open(RESULT,$results_file);
    my($sing,$core,$shared,$complete);

    while(<RESULT>){
	my @initial_values = split(/\t/,$_);
	map{$_ =~ s/\s+$//} @initial_values;
	my @id = shift(@initial_values);
	my @values;
	
	#Remove blank columns from @values
	for (my $empty = 0, my $not_empty = 0; $empty <= $#initial_values; $empty++){
	    if (length($initial_values[$empty]) == 0){
		#ignore blank columns which indicate there was no gene for this genome
		next;
	    }
	    $values[$not_empty++] = $initial_values[$empty];
	    #Find stored genome name
	    my $db = $ANNOT_HSH{$initial_values[$empty]}->{db};

	    if (!defined($db)){
		#All gene identifiers in the table should have been defined in the attribute file
		die ("$initial_values[$empty] was not defined in the attribute file $opts{data_file}\n");
	    }
	    $complete->{$id[0]}->{$db} = $initial_values[$empty];
	}

	#CHECK IF ADDING MORE THRESHOLDS MATTER HERE
	my $count = scalar @values;
	my $members = join(",",sort @values);
	
	$sing->{$id[0]} = $members if ($count == 1);
	$core->{$id[0]} = $members if ($count == $num_dbs);
	$shared->{$id[0]} = $members if ($count > 1);
    
	$INITIAL_COUNTS->{cluster}++;
    }

    close RESULT;
	 
    return ($sing,$core,$shared,$complete);
}
sub parse_role_id_lookup{
    my $file = shift;
    
    my $errors .= "File is size zero or does not exist: file\n" unless(-s $file);
    my ($hsh,$shsh);

    if($errors){
	die($errors);
    }else{
	my @role_lookup = read_file($file);
	
	foreach (@role_lookup){
	    my @values = split(/\t/,$_);
	    $values[1] =~ s/\s+$//;
	    $values[2] =~ s/\s+$//;
	    
	    $hsh->{$values[0]}->{main} = $values[1];
	    $hsh->{$values[0]}->{subrole} = $values[2]; 

	    $shsh->{$values[1]}->{$values[2]} = $values[0];
	}
    }

    return ($hsh,$shsh);
}
sub check_params {

    my $errors = '';

    if($opts{data_file}){

	$errors .= "$opts{data_file} does not exist or is size zero\n" unless(-s $opts{data_file});

    }else{

	$errors .= "Must provide the attribute data file, --data_file\n";

    }

    if($opts{method_result}){

	$errors .= "$opts{method_result} does not exist or is size zero\n" unless(-s $opts{method_result});

    }else{

	$errors .= "Must provide --method_result file\n";

    }

    if($opts{genomes_list}){

	$errors .= "$opts{genomes_list} does not exist or is size zero\n" unless(-s $opts{genomes_list});

    }else{

	$errors .= "Must provide --genomes_list file\n";

    }

    $ROLES_OF_INTEREST = ($opts{role_ids}) ? join("|", split(",",$opts{role_ids})) : "";

    @NAMES_OF_INTEREST = ($opts{terms}) ? split(",",$opts{terms}) : "";
    map{$_ =~ s/^\s+//} @NAMES_OF_INTEREST;

    $REFERENCES = ($opts{reference}) ? join("|", split(",",$opts{reference})) : "";

    $OUTPUT = $opts{output} // cwd();
    mkpath($OUTPUT) unless (-d $OUTPUT);

    $errors .= "File is size zero or does not exist: $opts{data_file}\n" unless(-s $opts{data_file});

    if ( $opts{fusion} ) {
        $errors .= "File is size zero or does not exist: $opts{fusion}\n" unless(-s $opts{fusion});
    }
    if ( $opts{combined_fasta} ) {
        $errors .= "File is size zero or does not exist: $opts{combined_fasta}\n" unless(-s $opts{combined_fasta});
    }
    if ( $opts{combined_seq} ) {
        $errors .= "File is size zero or does not exist: $opts{combined_seq}\n" unless(-s $opts{combined_seq});
    }
    die $errors if $errors;
}
