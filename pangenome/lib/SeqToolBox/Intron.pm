# $Id: Intron.pm 199 2008-03-07 20:29:26Z malay $
# Perl module for CDSFileParser
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2006 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

CDSFileParser  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Malay <malay@bioinformatics.org>


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##


package SeqToolBox::Intron;
use Carp;
use strict;


##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {   
	my $class = shift;
	my $self = {};
	bless $self, ref($class) || $class;
	my $retval = $self->_init(@_);
	#print STDERR "$retval\n";
	if ($retval ){
		return $self;
	}else{
		return;
	}
}  


# _init is where the heavy stuff will happen when new is called

sub _init {
	my($self,%opts) = @_;
	#print STDERR $opts{ftstring},"\n";
	croak "FTString not given" unless $opts{ftstring};
	$self->{strand} = $opts{strand} ? $opts{strand} : '+';
	my $string;
	
	if ($self->{strand} eq '-'){
		$string = $self->get_reverse_cds($opts{ftstring});	
	}else{
		$string = $opts{ftstring};
	}
	
	my $retval = $self->_parse($string);
	 #print STDERR "In init $retval\n";
	if (!$retval) {
		print STDERR "Error: parsing ", $opts{ftstring},"\n";
		return;
	}else {
		return $retval;
	}

}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

##-------------------------------------------------------------------------##

=head2 num()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub num {
	my $self = shift;
	my @num = keys %{$self->{coords}};
	#carp "@num\n";
	return scalar(@num);
}

##-------------------------------------------------------------------------##

=head2 get_normalized_exons_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_normalized_exons {
	my $self= shift;
	#croak "Index not given" unless $index;
	croak "No cds found" unless defined ($self->{coords});
	my @values = @{$self->{coords}};
	if ($self->{strand} eq "-") {
		return @values;
	}else{
		my $start = $values[0];
		return map {($_ - $start) + 1} @values;
	}
}

##-------------------------------------------------------------------------##

=head2 get_normalized_exon_start_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_normalized_exon_start {
	my $self = shift;
	#croak "Index not given" unless $index;
	croak "No cds found" unless defined ($self->{coords});
	my @values = $self->get_normalized_exons();
	my @return;
	for (my $i = 0; $i < @values; $i++){
		push @return, $values[$i];
		$i++;
	}
	return @return;
	
}

=head2 get_normalized_exon_end_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_normalized_exon_end {
	my $self = shift;
#	croak "Index not given" unless $index;
	croak "No cds found" unless defined ($self->{coords});
	my @values = $self->get_normalized_exons();
	my @return;
	for (my $i = 0; $i < @values; $i++){
		push @return, $values[$i+1];
		$i++;
	}
	return @return;
	
}


##-------------------------------------------------------------------------##

=head2 get_intron_phases_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_intron_phases {
	my $self = shift;
	#croak "Index not given" unless $index;
	croak "No cds found" unless defined ($self->{coords});
	#get_intron_pos_subtracting_intron_by_index
	#my @values = $self->get_normalized_exon_end_by_index($index);
	my @values = $self->get_intron_pos_subtracting_intron();
	#pop @values;
	if (@values){
	return map {$_ % 3} @values;
	}else{
		return;
	}
}
#* Add your methods here

##-------------------------------------------------------------------------##

=head2 get_intron_postions_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_intron_positions{
	my ($self,$index) = @_;
	#croak "Index not given" unless $index;
	croak "No cds found" unless defined ($self->{coords});
	my @values = $self->get_normalized_exon_end();
	pop @values;
	if (@values){
	return map {$_ + 1} @values;
	}else{
		return;
	}
	
}

##-------------------------------------------------------------------------##

=head2 get_intron_pos_subtracting_intron_by_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub get_intron_pos_subtracting_intron {
	my $self = shift;
	my @exons = $self->get_normalized_exons();
	my @return;
	my $deduct =0;
	if (@exons < 4){
		return;
	}
	for (my $i = 0; $i < @exons; $i++){
		if (!defined $exons[$i+2]){
		 last;
		}
		if ($i == 0){
			push @return, $exons[$i + 1];
			$i++;
			next;
		}
		$deduct += ($exons [$i] - $exons[$i - 1])  -1;
		push @return, ($exons[$i+1]) - $deduct;
		$i++;
		
	}
	#print STDERR "@return\n";
	return @return;
}

sub get_intron_lengths {
	my $self = shift;
	my @start = $self->get_normalized_exon_start();
	my @end = $self->get_normalized_exon_end();
	my @lengths;
	if (@start <2 || @end < 2){
		return;
	}else {
		for (my $i = 0; $i < @start - 1; $i++){
			push (@lengths, ($start[$i+1] - $end[$i]) - 1) ;
		}
	}
	return @lengths;
}

sub get_cds_length {
	my $self = shift;
	my @end = $self->get_normalized_exon_end();
	#print STDERR "@end\n";
	my $end = $end[$#end];
	#print STDERR "$end\n";
	my @start = $self->get_normalized_exon_start();
	#print STDERR "@start\n";
	my $start = $start[0];
	my $length = ($end - $start) + 1;
	#print STDERR $length,"\n";
	my @introns = $self->get_intron_lengths;
	#print STDERR "@introns\n";
	foreach my $i(@introns) {
		$length -= $i;
	}
	return $length;
}
##-------------------------------------------------------------------------##


=head1 PRIVATE METHODS

=cut

##-------------------------------------------------------------------------##

=head2 _parse()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut


sub _parse {
	my ($self,$ftstring) = @_;
	#print STDERR $ftstring,"\n";
	#open my $infile, "$file" || croak "Can't open $file";
	my $count = 1;
	#while (my $line = <$infile>){
	#	chomp $line;
	#	$line =~ s/\,$//;
		#$line =~ s/\,$//;
		#print STDERR $line,"\n";
	my @temp = split (/\,/,$ftstring);
		#print STDERR "@f\n";
	#	shift (@temp);
		#pop (@f);
		#print STDERR "@f\n";
#		foreach my $i (@f){
#			print "data","\t", $i,"\n";
#		}
		my @f;
		foreach my $i (@temp){
			push (@f, split(/\.\./, $i));
		}
		#print STDERR scalar(@f),"\n";
		if (@f % 2){
			print STDERR "Exon coords are not in pairs\n";
			print STDERR "Line: $ftstring\n";
			print STDERR "Split data: @f\n";
			return;
		}
		

		
		$self->{coords} = \@f;
		#$count++;
	#}
	#close $infile;
	#print STDERR "Returning 1\n";
	return 1;
}

sub get_reverse_cds {
	my ($self,$ft_string) = @_;
	my @start;
	my @end;
	my @ft     = split( /\,/, $ft_string );
	my $number = scalar(@ft);
	my $last_s = 0;
	my $last_e = 0;
	foreach my $f (@ft) {
		my ( $s, $e ) = split( /\.\./, $f );
		if ( $s > $e )     { die "Something wrong $ft_string\n"; }
		if ( $last_s > $s ) { die "Ordering problem in start\n"; }
		if ( $last_e > $e ) { die "Ordering problem in end\n"; }
		$last_s = $s;
		$last_e = $e;
		push @start, $s;
		push @end,   $e;
	}
	if ( scalar(@start) != scalar(@end) ) {
		die "Start and end should be equal in number\n";
	}
	my @reverse_start = reverse(@start);
	my @reverse_end   = reverse(@end);
	my $end_point     = $end[$#end];
	my @s;
	for ( my $i = 0 ; $i < @reverse_start ; $i++ ) {
		my $temp = ( $end_point - $reverse_end[$i]) + 1;
		$temp .= '..';
		$temp .= ( $end_point - $reverse_start[$i] ) + 1;

		push @s, $temp;

	}
	#print STDERR join( ",", $number, @s ),"\n";
	return join( ",", @s );

}


##-------------------------------------------------------------------------##




=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2006 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;  

