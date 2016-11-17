# $Id: Run.pm 38073 2012-12-19 20:35:23Z ekelsey $
# Perl module for Pangenome::Method::PanOCT::Run
# Author: Jason Inman <jinman@jcvi.org>
# Copyright (c) 2010 by JCVI. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

Pangenome::Method::PanOCT::Run  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jason Inman <jinman@jcvi.org>


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package Pangenome::Method::PanOCT::Run;
use vars qw(@ISA);
@ISA = qw(Pangenome::Method::RunI);
use Pangenome::Method::RunI;
use Carp;
use File::Copy;
use File::Path;
use File::Spec;
use File::Basename;
use strict;

##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	return $self;
}

# _init is where the heavy stuff will happen when new is called

sub _init {
	my ( $self, @args ) = @_;
	$self->{panoct_dir} =
		File::Spec->catdir( $FindBin::Bin, File::Spec->updir(), "ext",
							"PanOCT" );

}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

=head2 execute()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub execute {
	my $self = shift;
	
	unless ( $self->get_config()->get_current_dir() eq
			 $self->get_config()->get_working_dir() )
	{
		$self->_copy_fasta_files();
		$self->_copy_blast_file();
	}
	
	$self->_copy_gene_att_file();
	$self->_copy_crib_file();
	#warn "before running BTM";
	$self->_run_panoct();

	chdir( $self->get_config()->get_starting_dir() );

}

=head2 result_file()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub result_file {
	my ( $self, @args ) = @_;
	unless (exists $self->{result_file}) {
		croak "Could not find result file\nDid you run the execute method?";
		
	}
	return $self->{result_file};
}

sub get_panoct_dir {
	my $self = shift;
	return $self->{panoct_dir};
}

=head1 PRIVATE METHODS

=cut

sub _copy_fasta_files {
    my $self = shift;

    #Copy individual pep files
    foreach my $f ( $self->get_fasta_files() ) {
	
	copy( $f, $self->get_config()->get_working_dir())
	    || croak "Failed to copy $f to working dir: $!\n";
	
    }

    #Copy combined.fasta
    my $startdir = $self->get_config->get_starting_dir();
    my $combined_fasta = "$startdir/fastas/combined.fasta";

    if(-s $combined_fasta){
	copy($combined_fasta,$self->get_config()->get_working_dir()) ||
	    croak "Failed to copy $combined_fasta to working dir: $!\n";
    }
}

sub _copy_crib_file {
	my $self = shift;
	
	copy( $self->get_crib_file(), $self->get_config()->get_working_dir() )
	   || croak "Failed to copy crib file: ".$self->get_crib_file()." to working directory: ".
	                   $self->get_config()->get_working_dir()."\n$!\n";
	 
}

sub _copy_gene_att_file {
	my $self = shift;
	
	copy( $self->get_gene_att_file(), $self->get_config()->get_working_dir() ) 
	   || croak "Failed to copy gene_att file: ".$self->get_gene_att_file()." to working directory: ".
	                   $self->get_config()->get_working_dir()."\n$!\n";
}


sub _copy_blast_file {
	my $self = shift;
	
	copy( $self->get_blast_file(), $self->get_config()->get_working_dir() )
	   || croak "Failed to copy BLAST file: ".$self->get_blast_file()." to working directory: ".
	                   $self->get_config()->get_working_dir()."\n$!\n";
}

sub _get_combined_fasta {
 
    my $self = shift;
    my $fastas_ref = shift;
    my $workdir = $self->get_config->get_working_dir();
    
    my $combined_fasta = "$workdir/combined.fasta";
 
    if (!(-s $combined_fasta) ) {
        
        open (my $cfh, ">$combined_fasta") || croak "Can't open $combined_fasta: $!\n";
            
        foreach my $fasta_file (@{$fastas_ref}) {

            my $file_name = basename( $fasta_file );
            open(my $ffh, "<$workdir/$file_name") || croak "Can't open $workdir/$file_name: $!\n";
            print $cfh $_ while <$ffh>;	
            
        }

    }
    return ('combined.fasta');
}

sub _run_panoct {

    my $self = shift;
    #warn "run btm";
    my $workdir = $self->get_config()->get_working_dir();
    my $blast_file = basename( $self->get_blast_file() );
    my $crib_file = basename( $self->get_crib_file() );
    my $gene_att_file = basename( $self->get_gene_att_file() );
    my $percent_id = $self->get_percent_id();
    my $strict = $self->get_strict();

    my @fasta_files = $self->get_fasta_files();
    my $combined_fasta = $self->_get_combined_fasta(\@fasta_files);
    
    chdir( $workdir );
    
    #Malay
    my $panoct_exe = $self->get_config()->find_bin("panoct.pl");
    
    #Malay: Changing output redirection to a log file
    my $logfile = File::Spec->catfile ($workdir,'panoct.log');
            
    my $command = "perl $panoct_exe ";
    $command .= "-b $workdir ";
    $command .= "-t $blast_file ";
    $command .= "-f $crib_file ";
    $command .= "-g $gene_att_file ";
    $command .= "-P $combined_fasta ";
    $command .= "-i $percent_id " if $percent_id;
    $command .= "-S $strict " if $strict;
    $command .= "-L 1 -M Y -H Y -V Y -N Y -F 1.33 -G y -c \"0,25,50,75,100\" -T";

    #print "Running: $command\n";exit; 
    system("$command 2>$logfile")== 0 || croak "Error in running PanOCT";

    my $result_file = $workdir . '/' . 'matchtable.txt';
    my $frameshift_file = $workdir . '/' . 'frameshifts.txt';
    
    $self->{result_file} = ((-s $result_file) ? $result_file : undef);
    
}


=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2010 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;
