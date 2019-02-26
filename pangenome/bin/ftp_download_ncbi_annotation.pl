#!/usr/bin/env perl

###############################################################################
#                                                                             #
#       Copyright (C) 2016-2017 J. Craig Venter Institute (JCVI).             #
#       All rights reserved.                                                  #
#       Written by John Doe, Ph.D.                                            #
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

ftp_download_ncbi_annotation.pl

=head1 SYNOPSIS

    USAGE: ./ftp_download_ncbi_annotation.pl 

=head1 OPTIONS

B<--organism_search_term, -s>   :   Term to be used to limit search based on the organism field in the assembly summary file.

B<--biosample_list>             :   Provide a list of biosample IDs to be used as input.

B<--accession_list>             :   Provide a list of assembly accession IDs to be used as input.

B<--name_list>                  :   Provide a list of 'infraspecific species' names to download.

B<--cg, --wgs>                  :   Limit output to only complete genomes (B<--cg>) or incomplete genomes (B<--wgs>) based on 'assembly level' in the assembly sumamry file.  Default is to not limit output.

B<--working_dir, -w>            :   Directory to store various intermediate files, and look for/place the mapping, download, and log files by default.  Default is cwd.

B<--output_dir, -o>             :   Path to location where downloads will go.

B<--log_file, -l>               :   Path to a file containing warnings, etc.

B<--loglevel>                   :   0 for less detailed logs, 1 for more detailed.  Default is 0.

B<--no_downlaod>                :   Don't attempt to download anything, just create the output files.

B<--download_only>              :   No parsing of the assembly summary file, just download the files. [Requires B<--download_file>]

B<--download_file, -d>          :   Give a path to the file that contains a list of ids and ftp urls to download.

B<--mapping_file, -m>           :   Path to a file mapping ids to biosample ids.

B<--output_prefix>              :   Use the supplied term as the prefix for all log files, including download and mapping files.

B<--kingdom>                    :   Choose a kingdom within refseq from which to work in.  [Default is 'bacteria']

B<--section>                    :   Choose a section of the ftp site to download from.  Can be: 'refseq' [Default] or 'genbank', or 'both'.

B<--preserve_assembly_summary>  :   Save a copy of the retrieved assembly_summary.txt file

B<--assembly_summary_file>      :   Path to a previously saved assembly_summary.txt file.

B<--id_length>                  :   Max length of newly generated ids. [Default is 10]

B<--id_type>                    :   Used like "B<--id_type> biosample" to use bioample IDs instead of automatically generated IDs.

B<--min_N50>                    :   Minimum N50 for download. [ Not used by default ]

B<--max_contigs>                :   Omit genomes with more than this level of contigs.  [ Not used by default ]

B<--fasta,--gb,--both>          :   Specify the type of files to download.  [Default is B<--gb>]

B<--cleanup_ids>                :   For use with B<--download_only>, will cleanup ids found in B<--download_file>.

B<--separate_downloads>         :   Place each downloaded file into a directory within B<--output_dir> named by the newly generated ID for that genome.

B<--exclusion_list>             :   File containing a list of assembly accessions that should NOT be included, even if they match other criteria.

B<--help, -h>                   :   Print this help.

=head1 DESCRIPTION

This script is used to interact with a specific portion of the NCBI FTP site.  Specifically: genomes/refseq/bacteria
Under normal operation, the script does the following:

=over 1 

=item  1. Retrieve the assembly_summary.txt file from that directory

=item  2. Parse the file for rows matching various user-defined specifications

=item  3. Retrieve and parse assembly_stats.txt for potential download candidates

=item  4. Create several files based on that parsing

=item  5. Retrieve gb, fasta, or both types of file from the urls associated with any identified genomes.

=back

In addition to the normal operational mode, the user may specify B<--no_download>, which will perform the first four
steps.  It will still download the assembly_summary.txt file, but will not retrieve the gb or fasta files in step 5.  
The user might instead specify B<--download_only>, in which case the first four steps are skipped, and the script will
only attempt to download the files from the urls in the file specified with B<--download_file>.  If B<--cleanup_ids> is 
used, any non-alpha-numeric characters will be stripped out of the ids in the download file.

By default, only sequences from RefSeq will be downloaded.  This can be changed to GenBank by using B<--section genbank>,
which will instead download from the GenBank section of the NCBI FTP site.  When using a pre-defined list of assemblies to
download, that is, B<--biosample_list> or B<--accession_list>, a third option is available, called 'both', which will direct
the script to look in RefSeq first, and Genbank second, avoiding duplicates (based on the master wgs accession).

It is possible to save a local copy of the assembly_summary.txt file downloaded by using the B<--preserve_assembly_summary> flag,
and to use a specific local copy of an assembly_summary.txt file using B<--assembly_summary_file>.

=head1 INPUTS

The B<--organism_search_term> provides a string that will be interpreted as a regular expression used against the
fields in the assembly_summary.txt file that hold organism name information.  If the term supplied is of the form:
 "word word" (for example, "Ruminococcus gnavus")
the term will be modifed to allow for the very few cases where the organism_name field of assembly_summary.txt includes
square brackets around the genus name.  (For example, "[Ruminococcus] gnavus". )  If more advanced regex are supplied,
or if a regex was used that required fancy shell-escaping to get passed in, it might be useful to run this script in
B<--no_download mode>, examine the resulting .mapping and/or .download files, and once deemed satisfactory, proceed with
operation in B<--download_only> mode (or simply rerun without B<--no_download>, as the script doesn't mind writing over 
previous runs' output files)

As an alternative to B<--organism_search_term>, the user may provide a list of BioSample IDs with B<--biosample_list>, a
list of assembly accession IDs with B<--accession_list>, or a list of 'infraspecific species names' with B<--name_list>.
For any ID in those files, the script will find matches in the 'biosample', 'assembly accession', or 'infraspecific species'
columns of assembly_summary.txt and use matching rows to create the download file.  Note that unlike B<--organism_search_term>,
items in these alternative lists must match EXACTLY to the respective fields in assembly_summary.txt to be considered for download.

If B<--min_N50> is used to supply a minimum value for N50, the script will compare it with each non-complete genome's
contig-N50 value from its _assembly_stats.txt file, and omit assemblies that do not equal or exceed the given value.

If B<--max_contigs> is used to supply a maximum number of allowed contigs, the script will (again, for each non-complete
genome) compare the supplied value with the value in _assembly_stats.txt for 'contig-count'.  It will reject genomes 
that have more than B<--max_contgs> contigs, and will allow genomes if they have a value of 'contig-count' equal to or 
less than B<--max_contigs>.

If B<--exclusion_list> is provided, this script will build a list of the included assembly accesssions and will omit
any assemblies from that list, even if those assemblies match other criteria for exclusion.

=head1 OUTPUTS

There are three files produced under normal operation and when run in B<--no_download> mode:

=over 1

=item 1. B<Mapping File> - by default, named "[SECTION]_download_[TIMESTAMP].mapping", but may be supplied as B<--mapping_file>; This file maps Ids generated from the strain name (or, if not present, the organism name) to:

=over 2 

=item 1. Assembly accession

=item 2. Bioproject ID

=item 3. Biosample ID (if present, else "NA" )

=item 4. wgs master record (if present, else "NA" ))

=item 5. organism name as found in assembly_summary.txt

=item 6. strain fields as found in assembly_summary.txt

=item 7. contig_n50 from the assembly's _assembly_stats.txt file

=item 8. contig_count  from the assembly's _assembly_stats.txt file

=item 9. total_length from the assembly's _assembly_stats.txt file

=item 10. type-strains will have 'type strain' here

=back

contig_n50 and contig_count will contain 'complete' if the assembly is listed as a 'Complete Genome' in the assembly_sumamry.txt file, as those assemblies do not have those columns in _assembly_stats.txt

The auto-generated ID length can be changed from the default of 10 using the B<--id_length> option.
Alternatively, the script can produce output using Biosample IDs instead of the auto-generated IDs by specifying "B<--id_type> biosample" on the command line.

=item 2. B<Download File> - by default, named the same as mapping file, but ending in .download instead of .mapping (or whataver .extention might be present on the end of the B<--mapping_file> value); the name may also be supplied as B<--download_file>; two columns: id and the ftp_url.  Used as the list provided by B<--download_file> when B<--download_only> is invoked.

=back

=over 1

=item 3. B<Log File> - named the same as mapping file, but ending in .log instead of .mapping (or whatever .extension might be present on the end of B<--mapping_file> value); the name may also be specified as B<--log_file>; a file containing various warnings related to file retrieval and other issues arising during operation of the script.

=back

Note that log file will also be written during B<--download_only> operation, and still follows the same naming conventions.

In addition to these files, the downloaded files will be placed in the path specified as B<--output_dir>.  If the user specifies B<--separate_downloads>, each downloaded file is placed in a directory specific to the genome to which it belongs, named by the ID created within this script, found in the B<--output_dir>.  Additionally, B<--separate_downloads> also creates a list file for each type of file downlaoded.  (output_dir/fasta.list and/or output_dir/gb.list).

=head1 CONTACT

    Jason Inman
    jinman@jcvi.org

=cut

use Cwd;
use File::Basename;
use File::Copy;
use File::Fetch;
use File::Slurp;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Pod::Usage;
use LWP::UserAgent;

use Data::Dumper;

my $DEFAULT_KINGDOM = 'bacteria';
my $DEFAULT_SECTION = 'refseq';
my $DEFAULT_OUTPUT_BASENAME = 'download';
my $TODAY = get_date();
my @DEFAULT_DOWNLOAD_LIST = qw( );
my $DEFAULT_LOG_LEVEL = 0;
my $DEFAULT_ID_LENGTH = 10;
my $DEFAULT_ID_TYPE   = 'normal';
my $working_dir = getcwd();

my $MAX_DOWNLOAD_ATTEMPTS = 3;
my $SLEEP_TIME = 1;

my %opts;
GetOptions( \%opts,
            'accession_list=s',
            'assembly_summary_file=s',
            'biosample_list=s',
            'both',
            'cg',
            'cleanup_ids',
            'download_file|d=s',
            'download_only',
            'exclusion_list=s',
            'fasta',
            'gb',
            'id_length=i',
            'id_type=s',
            'kingdom|k=s',
            'log_file|l=s',
            'loglevel=i',
            'mapping_file|m=s',
            'max_contigs=i',
            'min_N50=i',
            'name_list=s',
            'no_download',
            'organism_search_term|s=s',
            'output_dir|o=s',
            'output_prefix|p=s',
            'preserve_assembly_summary',
            'section=s',
            'separate_downloads',
            'wgs',
            'working_dir|w=s',
            'help|h',
            ) || die "Error getting options! $!";
pod2usage( {-exitval => 0, -verbose => 2} ) if $opts{ help };
check_options();

# take care of filnames
$opts{ output_dir } = $opts{ output_dir } // getcwd() . '/output_dir';
mkdir $opts{ output_dir } unless ( -d $opts{ output_dir } );

my %exclusion_list;
if ( $opts{ exclusion_list } ) {

    open( my $efh, '<', $opts{ exclusion_list } ) || die "Can't open exclusion_list: $!\n";
    while ( <$efh> ) {
        chomp;
        my $fields = scalar( split(/\s/,$_) );
        if ( ( $_ !~ /^GC[AF]_/ ) || $fields > 1 ) {
            die "This line in the exclusion_list doean't look like it contains only a single assembly_accession:\n$_\n";
        }
        $exclusion_list{ $_ }++;

    }

}

my $mapping_file = $opts{ mapping_file } // $opts{ output_prefix } . '.mapping';
my $download_file = $opts{ download_file } // (fileparse( $mapping_file, qr/\.[^.]*/ ))[0] . '.download';
my $log_file = $opts{ log_file }    // (fileparse( $mapping_file, qr/\.[^.]*/ ))[0] . '.log';
open ( my $lfh, '>', $log_file ) || die "Problem opening log file: $log_file: $!\n";

# need these (maybe)
my %biosamples;
my %dupe_ids; # store the ord for the char to be used as a suffix when dupes arise
my $ua = LWP::UserAgent->new;

my $data;
# Use a saved version of the summary file...
if ( $opts{ assembly_summary_file } ) {
    $data = read_file( $opts{ assembly_summary_file } );
} else {

    # ... or get the summary file via ftp.
    my @sections;
    push @sections, 'refseq'  if ( $opts{ section } eq 'both' || $opts{ section } eq 'refseq'  );
    push @sections, 'genbank' if ( $opts{ section } eq 'both' || $opts{ section } eq 'genbank' );

    for my $section ( @sections ) {

        my $summary_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/'. $section . '/'. $opts{ kingdom } . '/assembly_summary.txt';
        _log("Retrieving assembly_summary.txt for $section\n", 1);
        my $response = $ua->get( $summary_url );
        if ( $response->is_success ) {
            $data .= $response->decoded_content;
        } else {
            die "ERROR getting assembly_sumamry.txt for $section: " . $response->status_line;
        }

    }
    write_assembly_summary( $data ) if $opts{ preserve_assembly_summary };

}

parse_data( $data ) unless ( $opts{ download_only } );
download_files( $download_file ) unless ( $opts{ no_download } );

exit(0);


sub write_assembly_summary {
# If asked, wil create a copy of the assembly_summary file used for this run.

    my ( $data ) = @_;

    my $assembly_summary_file = "$opts{ output_prefix }.assembly_summary.$opts{section}.txt";

    open( my $afh, '>', $assembly_summary_file ) || die "Can't write to $assembly_summary_file: $!\n";

    binmode $afh, ":encoding(UTF-8)";

    print $afh $data;

}


sub parse_data {

    my ( $data ) = @_;

    my %download;
    my %mapping;

    # load ids if necessary
    my ( $biosample_ids, $assembly_accs, $names );
    my %still_unseen; # Needed for --biosample_list checking.
    if ( $opts{ biosample_list } ) {

        $biosample_ids = load_input_ids( $opts{ biosample_list } ); 

        # Let me explain this next line... we want a hash with biosample ids as the keys
        # to both check for the existance of a given biosample and also make it easy to
        # remove from the set of keys to check.  By the end of the process of looking
        # for lines in assembly_summary.txt for matching lines, we'll have a hash with
        # keys that represent biosample ids that are missing.
        # Now, notice here we're assigning 'undef' to a hash slice (since we don't need values)
        # Also notice we're using an array reference to tell the hash slice which keys will be
        # set to undef.
        # So this is just assigning undef to a hash slice keyed by an arrary reference is all. :)
        @still_unseen{ @$biosample_ids } = undef;  

    }
    if ( $opts{ accession_list } ) {
        $assembly_accs = load_input_ids( $opts{ accession_list } );
    }
    if ( $opts{ name_list } ) {
        $names = load_input_ids( $opts{ name_list } );
    }

    my %matched_assemblies;
    my %candidates;
    my @field_list = qw/bioproject biosample wgs_master infraspecific_name isolate version_status assembly_level ftp_link type_strain/;
    for my $line ( split( "\n", $data ) ) {

        next if ( $line =~ /^#/ );

        my @line_array = split( "\t", $line );

        # Need columns 1-4, 8, 9, 11, 12, 20, 22 (index 0-3, 7, 8, 10, 11, 19, 21)
        my $assembly_accession  = $line_array[0];
        my $bioproject          = $line_array[1];
        my $biosample           = $line_array[2];
        my $wgs_master          = $line_array[3];
        my $infraspecific_name  = $line_array[7];
        my $isolate             = $line_array[8];
        my $version_status      = $line_array[10];
        my $assembly_level      = $line_array[11];
        my $ftp_link            = $line_array[19];
        my $type_strain         = $line_array[21];

        # check if it matches organism_search_term or our other search criteria
        if ( $opts{ organism_search_term } ) {

            # This hack is needed for the few unfortunate souls attempting to retrieve
            # genomes that have somehow made it into assembly_summary.txt with BRACKETS IN THE NAME.  Sheesh.
            if ( $opts{ organism_search_term } =~ /^([[:alpha:]]+)(\s+)([[:alpha:]]+)$/ ) { # Only do this when the search term is like "word word" to avoid gunking up other cases.
                $opts{ organism_search_term } = "([\\\[]?)$1([\\\]]?) $3";              # Surround the first word in the term with brackets that may or may not be there.
            }

            if ( $infraspecific_name =~ qr($opts{ organism_search_term }) ) {
                next if exists $matched_assemblies{ $assembly_accession };
                $matched_assemblies{ $assembly_accession }++;
            } else {
                next;
            }
        } elsif ( $opts{ biosample_list } ) {

            unless ( grep { /^$biosample$/ } keys %still_unseen ) {
                next;
            }
            # Don't need to check this anymore:
            delete $still_unseen{ $biosample };

        } elsif ( $opts{ accession_list } ) {
            unless ( grep { /^$assembly_accession$/ } @$assembly_accs ) {
                next;
            }
        } elsif ( $opts{ name_list } ) {
            unless ( grep { /^$infraspecific_name$/ } @$names ) {
                next;
            }
        }

        # only interested in latest genomes
        unless ( $version_status eq 'latest' ) {
            _log( "Skipping because non-latest:\t$line\n", 1 );
            next;
        }

        # might only be interested in complete genomes.
        if ( $opts{ cg } && $assembly_level ne "Complete Genome" ) {
            _log( "Skipping because --cg used and this isn't a Complete Genome:\t$line\n", 1 );
            next;
        }

        # might only be interested in incomplete genomes.
        if ( $opts{ wgs } && $assembly_level eq 'Complete Genome' ) {
            _log( "Skipping because --wgs used and this is a Complete Genome:\t$line\n", 1 );
            next;
        }

        # make sure we have a defined value for type_strain
        $type_strain = $type_strain // 'NA';
        $line_array[21] = $type_strain;

        # Also wgs_master
        $wgs_master = $wgs_master // 'NA';
        $line_array[3] = $wgs_master;

        # At this point, this genome has passed the first round of 'cuts'.  Move it to the list of genomes for which to download assembly_stats.txt
        # That is, unless it's on the --exclusion_list
        if ( exists $exclusion_list{ $assembly_accession } ) {
            _log( "Skipping because $assembly_accession is on the exclusion list:\t$line\n", 1 );
            next;
        }
        @{$candidates{ $assembly_accession }}{ @field_list } = @line_array[1,2,3,7,8,10,11,19,21];
        $candidates{ $assembly_accession }->{ line } = $line;
        
    }

    if ( $opts{ biosample_list } && scalar( keys %still_unseen ) )  {

        _warn( "The following biosamples were not found in assembly_summary.txt.  These links may help determine why:\n", 0 );
        for my $biosample ( keys %still_unseen ) {
            _warn( "$biosample\thttps://www.ncbi.nlm.nih.gov/assembly/?term=$biosample\n", 0 );
        }

    } 

    # Error if no candidates at this point.
    unless ( scalar %candidates ) {
        _warn( "No candidates for download found... please check inputs and try again.\n");
        exit(1);
    }

    filter_on_assembly_data( \%candidates, \@field_list, \%mapping, \%download );
    print_files( \%mapping, \%download );

}


sub filter_on_assembly_data {

    my ( $candidates, $field_list,  $mapping, $download ) = @_;

    my $astats_dir = "$working_dir/assembly_stats";
    mkdir $astats_dir unless ( -d $astats_dir );
    my $curr_dir = getcwd();
    chdir $astats_dir;

    # print curl config file:
    my $curl_config_file = "./$opts{ output_prefix }.assembly_stats_curl.config";
    open( my $cfh, '>', $curl_config_file ) || die "Can't open $curl_config_file: $!\n";

    my %seen;
    my $attempts = 0;

    while ( $attempts < $MAX_DOWNLOAD_ATTEMPTS && (scalar( keys %seen) != scalar( keys %$candidates ) ) ) {

        $attempts++;

        for my $assembly_accession ( keys %$candidates ) {

            if ( not exists $seen{ $assembly_accession } ) {

                # Neat how we can use File::Basename to parse urls, too :)
                my $ftp_link = $candidates->{ $assembly_accession }->{ ftp_link };
                my $genome_id = basename( $ftp_link );
                my $stats_file_url = "$ftp_link/$genome_id" . "_assembly_stats.txt";
                $candidates->{ $assembly_accession }->{ stats_file } = basename( $stats_file_url );

                print $cfh qq(-O\nurl = "$stats_file_url"\n);

            } 

        }

        # grab curl data:
        system( "curl -K $curl_config_file >& ./$opts{ output_prefix }.curl_log" );

        # check for retries:
        for my $assembly_accession ( keys %$candidates ) {
            my $expected_file = $candidates->{ $assembly_accession }->{ stats_file };
            $seen{ $assembly_accession }++ if ( -s $expected_file );
        }

    }    
        
    if ( ! scalar( keys %seen ) ) {
        die "Couldn't get any assembly_stats.txt files.  There is nothing more to do.";
    } elsif ( scalar( keys %seen ) != scalar( keys %$candidates ) ) {
        _warn("Looks like there were some assembly_stats.txt files that did not get pulled. These will be skipped:\n",0);
        for my $assembly_accession ( keys %$candidates ) {
            my $expected_file = $candidates->{ $assembly_accession }->{ stats_file };
            unless ( -s $expected_file ) {
                _warn( "$expected_file\n",0);
                delete $candidates->{ $assembly_accession };
            }
        }
    }

    chdir $curr_dir;

    my $dupe_index = 1;

    for my $assembly_accession ( keys %$candidates ) {

        my $assembly_stats_file = "$astats_dir/" . $candidates->{ $assembly_accession }->{ stats_file } ;
        my @asmbl_data = read_file( $assembly_stats_file );

        # Just a little assignment of values into a list using an array reference as a hash reference slice:
        my ( $bioproject, $biosample, $wgs_master, $infraspecific_name, $isolate, $version_status, $assembly_level, $ftp_link, $type_strain ) = 
                @{$candidates->{ $assembly_accession }}{ @$field_list };
        my $line = $candidates->{ $assembly_accession }->{ line };

        my ( $N50, $contig_count, $total_length, $is_complete );
        $is_complete = ( $assembly_level eq 'Complete Genome' ) ? 1 : 0;

       ( $N50, $contig_count, $total_length ) = get_assembly_stats( \@asmbl_data, $is_complete );

        if ( $opts{ min_N50 } && ! $is_complete ) {

            if ( $N50 < $opts{ min_N50 } ) {
                # failed.  log it.
                _log( "$infraspecific_name rejected because it had N50 of $N50, lower than requirement of $opts{ min_N50 }\n", 0 );
                delete $candidates->{ $assembly_accession };
                next;
            }

        }

        if ( $opts{ max_contigs } && ! $is_complete ) {

            if ( $contig_count > $opts{ max_contigs } ) {
                # failed.  log it.
                _log( "$infraspecific_name rejected because it had $contig_count contigs, more than the maximum of $opts{ max_contigs }\n", 0 );
                delete $candidates->{ $assembly_accession };
                next;
            }

        }

        my $id = check_row( $biosample, $infraspecific_name, $isolate, $line, $download, $mapping, \$dupe_index );

        if ( $id ) {
            add_rows( $download, $mapping, $id, $assembly_accession, $bioproject, $wgs_master, $biosample, $ftp_link, $infraspecific_name, $isolate, $N50, $contig_count, $total_length, $type_strain );
        }

        unlink $assembly_stats_file;

    }
}


sub check_row {

    my ( $biosample, $infraspecific_name, $isolate, $line, $download, $mapping, $dupe_index ) = @_;
    my $new_id = undef;

    $new_id = create_new_id( $isolate, $infraspecific_name, $line, $download, $mapping, $dupe_index );

    if ( $biosample =~ /^$/ ) {
        # log missing biosample
        _log( "Missing biosample from line for $new_id, $infraspecific_name.\n", 0 );
    } elsif ( $opts{ id_type } eq 'biosample' ) {
        $new_id = $biosample;
    }

    return $new_id;

}


sub create_new_id {

    my ( $isolate, $infraspecific_name, $line, $download, $mapping, $dupe_index ) = @_;
    my $new_id;

    # Attempt to give short, somewhat meaningful ids:
    if ( $isolate =~ /strain=(.*)/ ) {
        # get id from $1
        $new_id = $1;
    } elsif ( $isolate =~ /^$/ ) {
        # get id from name
        $new_id = $infraspecific_name;
    } else {
        # log unparsable strains
        _log( "Unable to parse strain from line: $line.\n", 0 );
    }

    # strip out undesirables. then take that last x whatever is left, as per id_length.
    $new_id =~ s/[^[:alnum:]]//g;
    $new_id = substr( $new_id, 0 - $opts{ id_length }, $opts{ id_length } );

    # make sure we have an id that isn't already used.
    $new_id = resolve_dupe( $new_id, $download, $mapping, $dupe_index, $isolate );

    return $new_id;

}


sub resolve_dupe {

    my ( $id, $download, $mapping, $dupe_index, $isolate ) = @_;

    if ( exists $download->{ $id } || exists $dupe_ids{ $id } ) {

        # if we just discovered this, we'll need to correct the initial instance of id in
        # the download_file AND the mapping file
        unless ( exists $dupe_ids{ $id } ) {

            my $new_id;
            $dupe_ids{ $id } = 'a';

            # add a lowercase letter suffix that increments, and shift everything over to preserve length.
            do {

                $new_id = $id . $dupe_ids{ $id };
                $new_id = substr( $new_id, 0 - $opts{ id_length }, $opts{ id_length } );
                $dupe_ids{ $id }++;

            } until ( not exists $download->{ $new_id } ); 

            $download->{ $new_id } = $download->{ $id };
            delete $download->{ $id };
            $mapping->{ $new_id } = $mapping->{ $id };
            delete $mapping->{ $id };

            # Make note of the new id for this id:
            _warn( "WARNING! Attempt to use $id as id more than once! Renaming to $new_id\n", 0);

        }

        # rename newly found instance of this id and warn the user
        my $new_id;
        do {

            $new_id = $id . $dupe_ids{ $id };
            $new_id = substr( $new_id, 0 - $opts{ id_length }, $opts{ id_length } );
            $dupe_ids{ $id }++;

        } until ( not exists $download->{ $new_id } );

        _warn( "WARNING! Attempt to use $id as id more than once! Renaming to $new_id\n", 0);

        $id = $new_id;

    }

    return $id;


}


sub get_assembly_data {

    my ( $ftp_link ) = @_;

    # Neat how we can use File::Basename to parse urls, too :)
    my $genome_id = basename( $ftp_link );
    my $stats_file_url = "$ftp_link/$genome_id" . "_assembly_stats.txt";

    my $assembly_stats = retrieve_assembly_stats_file( $stats_file_url );
    return undef unless ( defined $assembly_stats );

    my @asmbl_data = split( "\n", $assembly_stats );

    return \@asmbl_data;

}


sub retrieve_assembly_stats_file {

    my ( $stats_file_url ) = @_;

    my $data;
    my $tries = 0;

    while ( $tries <= $MAX_DOWNLOAD_ATTEMPTS ) {

        $tries++;

        my $data;
        my $response = $ua->get( $stats_file_url );
        if ( $response->is_success ) {
            $data = $response->decoded_content;
        } else {
            _warn( "Didn't get $stats_file_url: ".$response->status_line."\n", 0 );
        }
        sleep $SLEEP_TIME; 

        if ( defined $data ) {

            last;

        } else {

            _log( "No data for $stats_file_url after $tries attempt(s)\n", 0 );
            if ( $tries > $MAX_DOWNLOAD_ATTEMPTS ) {
                _warn( "Exceeded max number of download attempts ( $MAX_DOWNLOAD_ATTEMPTS ) for $stats_file_url\n", 0 );
            }

        }

    }

    return $data;

}


sub get_assembly_stats {

    my ( $asmbl_data, $is_complete ) = @_;
    my ( $N50, $contig_count, $total_length );

    if ( $is_complete ) {
        $N50 = 'complete';
        $contig_count = 'complete';
    }

    for my $line ( @$asmbl_data ) {

        next if ( $line =~ /^#/ );
        if ( $line =~ /all\tall\tall\tall\tcontig-count\t(\d+)/ ) {
            $contig_count = $1;
        }
        if ( $line =~ /all\tall\tall\tall\tcontig-N50\t(\d+)/ ) {
			$N50 =  $1;	
		}
        if ( $line =~ /all\tall\tall\tall\ttotal-length\t(\d+)/ ) {
            $total_length = $1;
        }

        last if ( defined $N50 && defined $contig_count && defined $total_length );

    }

    return ( $N50, $contig_count, $total_length );

}


sub add_rows {

    my ( $download, $mapping, $id, $assembly_accession, $bioproject, $wgs_master, $biosample, $ftp_link,
         $infraspecific_name, $isolate, $N50, $contig_count, $total_length, $type_strain ) = @_;

    # store line for the .downlaod file
    $download->{ $id } = $ftp_link;


    # store line for the .mapping file
    if ( exists $biosamples{ $biosample } ) {

        unless ( $biosample =~ /^$/ ) {

            # whooooopsie.  log duplicate biosamples.
            _warn( "WARNING! Duplicate biosample: $biosample\n", 0 );

        }

    }

    $isolate =~ s/strain=(.*)/$1/;

    $mapping->{ $id } = "$assembly_accession\t$bioproject\t$biosample\t$wgs_master\t$infraspecific_name\t$isolate\t$N50\t$contig_count\t$total_length\t$type_strain";
    $biosamples{ $biosample }++;

}


sub print_files {

    my ( $mapping, $download ) = @_;

    open( my $mfh, '>', $mapping_file ) || die "Can't open $mapping_file: $!\n";
    print $mfh map { "$_\t$mapping->{ $_ }\n" } keys %$mapping;

    open( my $dfh, '>', $download_file ) || die "Can't open $download_file: $!\n";
    print $dfh map { "$_\t$download->{ $_ }\n" } keys %$download;

}


sub download_files {

    my ( $download_file ) = @_;

    open( my $ofh, '<', $download_file ) || die "Can't open $download_file: $!\n";

    my $curl_config_file = "./$opts{ output_prefix }.output_file_curl.config";
    open( my $cofh, '>', $curl_config_file ) || die "Can't open $curl_config_file: $!\n";

    my %downloads;
    while (<$ofh>) {

        next if ( /^Identifier/ );

        chomp;
        my ( $id, $ftp_url ) = split( /\t/, $_ );
        cleanup_id( \$id ) if ( $opts{ cleanup_ids } );

        my $new_gb_file = $opts{ output_dir };
        if ( $opts{ separate_downloads } ) {
            $new_gb_file .= "/$id";
            mkdir "$opts{ output_dir }/$id" unless ( -d "$opts{ output_dir }/$id" );
        }
        $new_gb_file .= "/$id.gb";
        my $new_fasta_file = $opts{ output_dir };
        if ( $opts{ separate_downloads } ) {
            $new_fasta_file .= "/$id";
            mkdir "$opts{ output_dir }/$id" unless ( -d "$opts{ output_dir }/$id" );
        }
        $new_fasta_file .= "/$id.fasta";

        my $ftp_gb_url;
        my $ftp_fasta_url;
        if ( $ftp_url =~ /.*\/([^\/]+)$/ ) {
            $ftp_gb_url = $ftp_url . '/' . $1 . '_genomic.gbff.gz';
            $ftp_fasta_url = $ftp_url . '/' . $1 . '_genomic.fna.gz';
        } else {
            # log unable to parse ftp url
            _log( "Unable to parse ftp url from $ftp_url for $id", 0 );
        }

        if ( $opts{ gb } || $opts{ both } ) {

            print $cofh qq(-O\nurl= "$ftp_gb_url"\n);
            $downloads{ $ftp_gb_url }->{ target_file } = $new_gb_file;
            $downloads{ $ftp_gb_url }->{ id } = $id;

        }

        if ( $opts{ fasta } || $opts{ both } ) {

            print $cofh qq(-O\nurl= "$ftp_fasta_url"\n);
            $downloads{ $ftp_fasta_url }->{ target_file } = $new_fasta_file;
            $downloads{ $ftp_fasta_url }->{ id } = $id;

        }

    }

    retrieve_gb_and_fasta_files( \%downloads );

}


sub retrieve_gb_and_fasta_files {

    my ( $downloads ) = @_;

    # These are only needed when --separate_downloads is used
    my ( $flfh, $glfh, %fasta_list, %gb_list );
    
    if ( $opts{ separate_downloads } ) {
        open( $flfh, '>', "$opts{ output_dir }/fasta.list" ) || die "Can't open fasta.list for writing!";
        open( $glfh, '>', "$opts{ output_dir }/gb.list"    ) || die "Can't open gb.list for writing.";
    }

    my %seen;
    my $attempts = 0;

    while ( $attempts < $MAX_DOWNLOAD_ATTEMPTS && (scalar( keys %seen ) != scalar( keys %$downloads )) ) {

        $attempts++;

        my $curl_config_file = "./$opts{ output_prefix }.output_file_curl.config";
        my $cofh;
        if ( $attempts != 1 ) {
            open( $cofh, '>', $curl_config_file ) || die "Can't open $curl_config_file: $!\n";
            for my $url ( keys %$downloads ) {
                if ( not exists $seen{ $url } ) {
                    print $cofh qq(-O\nurl = "$url"\n);
                }
            }
        }

        # grab files with curl
        system( "curl -K $curl_config_file >& ./$opts{ output_prefix }.curl_output_files_log" );

        # Check for retries:
        for my $url ( keys %$downloads ) {
            my $expected_file = basename( $url );
            $seen{ $url }++ if ( -s $expected_file );
        }

    }

    # follow up about missing ones here:
    if ( ! scalar( keys %seen ) ) {
        die "Couldn't get any gb/fasta files.  There is nothing more to do.";
    } elsif ( scalar( keys %seen ) != scalar( keys %$downloads ) ) {
        _warn( "Looks like there were some files that did not get downloaded.  These will be skipped:\n", 0 );
        for my $url ( keys %$downloads ) {
            my $expected_file = basename( $url );
            unless ( -s $expected_file ) {
                _warn( "$expected_file\n", 0 );
                delete $downloads->{ $url };
            }
        }
    }

    # move them to where they should be now.
    for my $url ( keys %$downloads ) {

        my $zipped_file = basename( $url );

        my $target_file = $downloads->{ $url }->{ target_file };

        # unzip and rename in one move, then delete the zipped file
        gunzip $zipped_file => $target_file or die "gunzip failed: $GunzipError\n";
        unlink $zipped_file;

        if ( $target_file =~ /\.gb$/ ) {
            # CDS check
            my $expected_cds = get_expected_cds( $target_file );
            my $found_CDS = count_cds( $target_file );
            if ( $found_CDS ) {

                # If we don't have an expected cds, print what we have and move on.
                if ( $expected_cds == 0 ) {
		    _log( "Don't know how many CDS to look for in $target_file, but found $found_CDS\n", 0 );
                } elsif ( $found_CDS != $expected_cds ) {
		    _log( "Found $found_CDS CDS but expected $expected_cds CDS from $target_file!\n", 0 );
                } elsif ($found_CDS == $expected_cds){
		    _log ("Found $found_CDS CDS, same as expected $expected_cds from $target_file\n", 0 );
		}

            } else {
                _log( "Found ZERO CDS features in $target_file!\n", 0 );
            }

            # create the list file of all gb files if asked.
            if ( $opts{ separate_downloads } ) {
                my $id = $downloads->{ $url }->{ id };
                my $abs_gb_file = Cwd::abs_path( $target_file );
                print $glfh "$id\t$abs_gb_file\n";
            }

        } elsif ( $target_file =~ /\.fasta$/ ) {

            # create the listfile of all fasta files if asked.
            if ( $opts{ separate_downloads } ) {
                my $id = $downloads->{ $url }->{ id };
                my $abs_fasta_file = Cwd::abs_path( $target_file );
                print $flfh "$id\t$abs_fasta_file\n";
            }

        }

    }
}


sub fetch_file {

    my ( $ftp_url ) = @_;

    my $retval;
    my $attempt = 1;

    until ( $retval || $attempt > $MAX_DOWNLOAD_ATTEMPTS ) {

        my $ff = File::Fetch->new( uri => $ftp_url );
        $retval = $ff->fetch( to => $opts{ output_dir } );
        $attempt++ unless $retval;

        sleep $SLEEP_TIME;

    }
if ( $attempt > 1 && $retval ) {
    print "Took $attempt tries, but we got it :)\n";
}

    return $retval;

}



sub get_expected_cds {

    my ( $file ) = @_;
    my $expected_cds_count = 0;

    open( my $cds_fh, '<', $file ) || die "Error opening $file!\n$!\n";
    my $total_flag = 0;

    while ( <$cds_fh> ) {

	#Match if CDS includes total
	if ( / {12}CDS\s+\(total\)\s+[:]{0,2}\s+([\d]+\,?[\d]+)/ ) {
            ( $expected_cds_count = $1 ) =~ s/,//g; # don't forget to strip out the commas.
	    $total_flag = 1;

	#Match if CDS line does not include total
        } elsif (/ {12}CDS\s+[:]{0,2}\s+([\d]+\,?[\d]+)/ ){
	    ( $expected_cds_count = $1 ) =~ s/,//g; # don't forget to strip out the commas.

	#Add pseudo gene count for non total CDS count
	} elsif ( / {12}Pseudo Genes\s+[:]{0,2}\s+([\d]+\,?[\d]+)/ ) { 

	    unless($total_flag){
		( my $pseudos = $1 ) =~ s/,//g;
		$expected_cds_count += $pseudos;
		last;
	    }
        }

    }

    return $expected_cds_count;

}


sub count_cds {

    my ( $file ) = @_;

    my $cds_count = 0;

    open( my $cds_fh, '<', $file ) || die "Error opening $file!\n$!\n";

    while ( <$cds_fh> ) {
        $cds_count++ if /^ {5}CDS/;
    }

    return $cds_count;

}


sub get_date {

    my @now = localtime();

    my $timestamp = sprintf( "%04d%02d%02d", $now[5] + 1900, $now[4] + 1, $now[3] );
    
    return $timestamp;

}


sub print_parsed_row {
# Useful for debugging the input rows.

    my ( $biosample, $infraspecific_name, $isolate, $version_status, $assembly_level, $ftp_link ) = @_;

    print<<EOF
BIOSAMPLE: $biosample
INFRASPECIFIC NAME: $infraspecific_name
ISOLATE: $isolate
VERSION STATUS: $version_status
ASSEMBLY LEVEL: $assembly_level
FTP LINK: $ftp_link
EOF

}


sub load_input_ids {

    my ( $id_file ) = @_;

    my @ids;

    open( my $ifh, '<', $id_file ) || die "Can't open $id_file: $!\n";

    while ( <$ifh> ) {

        chomp;
        $_ =~ s/(.*\S)\s+/$1/;
        push @ids, $_;

    }

    return \@ids;

}


sub cleanup_id {
# for when you give it something you don't really mean as the id.

    my ( $id_ref ) = @_;

    $$id_ref =~ s/[^[:alnum:]]//g;

}


sub _log {

    my ( $msg, $level ) = @_;

    $level = 0 unless $level;

    print $lfh $msg if ( $level >= $opts{ loglevel } );

}


sub _warn {

    my ( $msg, $level ) = @_;

    $level = 0 unless $level;

    # warn AND log!
    if ( $level >= $opts{ loglevel } ) {
        warn $msg;
        _log( $msg, $level );
    }

}


sub check_options {

    my $errors = '';

    if ( $opts{ both } ) {

        $opts{ fasta }++;
        $opts{ gb }++;

    } 

    # Default to --gb
    if ( !( $opts{ both } || $opts{ fasta } ) ) {
        $opts{ gb }++;
    }

    $errors .= "Use only one of --cg or --wgs.  To download both, use neither.\n" if ( $opts{ cg } && $opts{ wgs } );

    unless ( $opts{ no_download } ) {

        $errors .= "Need --fasta, --gb, or --both\n" unless ( $opts{ fasta } || $opts{ gb } );

    }

    if ( $opts{ download_only } ) {

        $errors .= "When --download_only is used, please give an input --download_file\n" unless ( $opts{ download_file } );

    } elsif ( ! ( $opts{ organism_search_term } || $opts{ biosample_list } || $opts{ accession_list } || $opts{ name_list } ) ) {

        $errors .= "Please supply one of: --organism_search_term, --biosample_list, or --accession_list\n";
    }

    # set up some default values
    $opts{ id_length }      = $opts{ id_length }        // $DEFAULT_ID_LENGTH;
    $opts{ loglevel }       = $opts{ loglevel }         // $DEFAULT_LOG_LEVEL;
    $working_dir            = $opts{ working_dir }      // $working_dir;
    $opts{ kingdom }        = $opts{ kingdom }          // $DEFAULT_KINGDOM;
    unless ( $opts{ kingdom } =~ /^(bacteria|archaea|fungi|invertebrate|plant|protozoa|vertebrate_mammalian|vertebrate_other|viral)$/ ) {
        $errors .= "--kingdom MUST be one of the following:\n\tbacteria (default)\n\tarchaea\n\tfungi\n\tplant\n\tprotozoa\n\tvertebrate_mammalian\n\tvertebrate_other\n\tviral\n";
    } 
    $opts{ section }        = $opts{ section }          // $DEFAULT_SECTION;
    unless ( $opts{ section } =~ /^(refseq|genbank|both)$/ ) {
        $errors .= "--section MUST be either 'refseq' or 'genbank' or 'both'\n";
    } 
    $opts{ id_type }        = $opts{ id_type }          // $DEFAULT_ID_TYPE;
    unless ( $opts{ id_type } =~ /^(normal|biosample|biosamples)$/ ) {
        $errors .= "--id_type MUST be either 'normal' or 'biosample'\n";
    }
    $opts{ id_type } = 'biosample' if $opts{ id_type } eq 'biosamples';

    if ( $opts{ exclusion_list } ) {

        $errors .= "Can't find --exclusion_list $opts{ exclusion_list } (or it's empty)\n"
            unless ( -s $opts{ exclusion_list } );

    }

    $opts{ output_prefix }  = $opts{ output_prefix }    // "$opts{ section }_$DEFAULT_OUTPUT_BASENAME.$TODAY";


    die $errors if $errors;

}
