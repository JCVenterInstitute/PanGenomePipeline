#!/usr/bin/env perl

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

use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Cwd;
use Data::Dumper;
use File::Slurp;
use File::Path;
use FindBin;
use lib "$FindBin::Bin/../lib";
use TIGR::FASTArecord;
use TIGR::FASTAreader;
use TIGR::FASTAwriter;
use TIGR::Foundation;

my %opts;
GetOptions(\%opts,'fasta_file|f=s',
    'max_num|n=i',
    'output|o=s',
    'extension|e=s',
    'help|h') || die "Couldn't get options! $!\n";

my $tf_object = new TIGR::Foundation;
my @errors;

my $seq_count = 0;
my $file_count = 1;
my $output;
my $extension;

&check_params;

my $fr = new TIGR::FASTAreader( $tf_object, \@errors, $opts{fasta_file} );
my $fw = new TIGR::FASTAwriter( $tf_object, \@errors );
my $filename = "$output/split_fasta.$file_count";
$filename .= $extension if $extension;
$fw->open( $filename );

while ( $fr->hasNext() ) {

    my $seq_obj = $fr->next();

    if( $seq_count < $opts{max_num} ) {

        $fw->write( $seq_obj );
        $seq_count++;

    } else {

        #Close and open new file
        $file_count++;
        $seq_count = 0;
        $filename = "$output/split_fasta.$file_count";
        $filename .= $extension if $extension;
        $fw->open( $filename );
        $fw->write( $seq_obj );
        $seq_count++;

    }

}

exit(0);

sub check_params {

    my $usage = "\nUsage: ./split_fasta.pl --fasta_file <fasta file to split> --max_num <max number of sequences per file> --output <location> --extension <new fasta extension>\n";
    
    if ( !($opts{fasta_file}) || !($opts{max_num}) || $opts{help} ) {

        die "$usage\n";

    }

    $extension = $opts{ extension } // undef;

    $output = $opts{ output } // cwd;
    
    mkpath( $output ) unless( -d $output );

}
