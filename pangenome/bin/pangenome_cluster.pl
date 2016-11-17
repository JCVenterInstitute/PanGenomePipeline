#!/usr/bin/env perl
# $Id: pangenome_cluster.pl 38722 2013-06-20 13:19:24Z ekelsey $

##---------------------------------------------------------------------------##
##  File: pangenome_cluster.pl
##   Author:  Malay <malay@bioinformatics.org>
#******************************************************************************
#* Copyright (C) 2010 Malay K Basu <malay@bioinformatics.org>
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

pangenome_cluster.pl - This script is a runner script that runs different clustering methods on the files provided.

=head1 SYNOPSIS

pangenome_cluster.pl -m inparanoid -b all_vs_all_blast.bla genome1.fas genome2.fas ...


=head1 DESCRIPTION

The script accepts a file containing all vs. all blast result file and a set of fasta files to create the clusters using the method of choice. The result are dumped onto standard output.


=head1 ARGUMENTS 

The program arguments are a set of fasta files. You can even use a wildcard like *.fas.

=head1 OPTIONS

=over 4

=item B<--blast_file|-b all_vs_all.bla>

The file containing all_vs_all blast result.

=item B<--percent_id|-i 35>

percent_id used as blast cutoff score in method algorithm

=item B<--strict|-s >

level of strictness in method algorithm

=item B<--method|-m method>

Method can be either inparanoid or orthomcl or panoct. 

=item B<--working_dir|-w working_dir >

The working dir where the files will be saved. If ommited the current dir is taken as working dir.

=item B<--keep_files|-k>

Keep the temporary files. By default the software will delete all the temporary files created.

=back


=head1 COPYRIGHT

Copyright (c) 2010 Malay K Basu <malay@bioinformatics.org>

=head1 AUTHORS

Malay K Basu <malay@bioinformatics.org>

=cut

##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../CPAN/lib";
use Getopt::Long;
use Pod::Usage;
use Pangenome::Config;
use Cwd;
use File::Path qw(make_path remove_tree);
use Pangenome::MethodFactory;
use Pangenome::ParserFactory;
use Carp;
use strict;
use warnings;

##---------------------------------------------------------------------------##
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
##---------------------------------------------------------------------------##

my %options = ();    # this hash will have the options

#
# Get the supplied command line options, and set flags
#


GetOptions( \%options, 
	    'help|?',
	    'method|m=s',
	    'blast|b=s',
	    'working_dir|w=s',
	    'crib|c=s',
	    'gene_att|g=s',
	    'percent_id|i=s',
	    'strict|s=s',
	    'keep_files|k' ) || pod2usage( -verbose => 0 );


my @fasta_files = @ARGV;
_check_params( \%options );


my $config = Pangenome::Config->new();
my $current_dir = cwd();
$config->set_starting_dir($current_dir);
my $working_dir = $options{working_dir} ? $options{working_dir} : cwd();
my $temp_working_dir = File::Spec->catdir($working_dir, "results");
make_path($temp_working_dir);
$config->set_working_dir ($temp_working_dir);
my $method = Pangenome::MethodFactory->new($options{method});
croak "Could not create object for method ". $options{method}. "\n" unless $method;
$method->set_fasta_files (@fasta_files);
$method->set_blast_file($options{blast});
$method->set_crib_file($options{crib});
$method->set_gene_att_file($options{gene_att});
$method->set_percent_id($options{percent_id});
$method->set_strict($options{strict});
$method->set_config($config);
$method->execute();
my $result_file = $method->result_file();

#print STDERR $result_file, "\n";
my $parser = Pangenome::ParserFactory->new($options{method});
croak "Could not create parser" unless $parser;


$parser->parse_file ($result_file);
$parser->dump(); 

unless (exists $options{keep_files}) {
    #remove_tree ($temp_working_dir,0,1);
}

exit(0);

######################## S U B R O U T I N E S ############################

sub _check_params {
    my $opts = shift;
    pod2usage( -verbose => 2 ) if ( $opts->{help} || $opts->{'?'} );
    pod2usage( -verbose => 1 ) unless ( $opts->{'blast'} );
    unless ($opts->{method}) {print STDERR "Method not set"; croak;}
    if ($opts->{method} =~ /panoct/i) {
	unless ($opts->{crib} && $opts->{gene_att}) {
	    warn 'Running "PanOCT" requires setting crib and gene attribute file', "\n";
	    pod2usage(-verbose=>1);
	}
    }
}
