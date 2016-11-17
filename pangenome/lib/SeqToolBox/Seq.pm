#$Id: Seq.pm 227 2008-03-28 18:32:57Z malay $
# Perl module for SeqToolBox::Seq
# Author: Malay < malay@bioinformatics.org >
# Copyright (c) 2006 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

=head1 NAME

SeqToolBox::Seq - Class for creating Seq object

=head1 SYNOPSIS

	my $seq = SeqToolBox::Seq->new(-id => "myseqid", -des=> "mydes", -seq=>"ATGC");
	my $lenght = $seq->length();
	...


=head1 DESCRIPTION

The class contains three fields:
	id 	= Seq id
	des = Seq description line
	seq = Sequence string

=cut

package SeqToolBox::Seq;
@ISA = qw(SeqToolBox::Root);
use SeqToolBox::Root;
use Carp;
use strict;
use fields qw(id des seq length);

=head1 CONSTRUCTOR

	Usage: my $seq = SeqToolBox::Seq->new(-id => "myseqid", -des=> "mydes", -seq=>"ATGC");
	
	Args :  
		id 	= A seq id
		des	= A description of the sequence
		seq	= The sequence string
	
	Returns: A seq object

=cut

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);
	return $self;
}

sub _init {
	my ( $self, @args ) = @_;
	my ( $id, $des, $seq )
	  = $self->_rearrange( [ "ID", "DES", "SEQ" ], @args );

	#print STDERR "DEBUG:$des\n";
	$self->{id}  = $id  || "";
	$self->{des} = $des || "";
	$self->{seq} = $seq || "";
}

#sub _parse_file {
#    my $self = shift;
#    my $file = shift;
#    open (FILE, $file) || croak "Can't open $file\n";
#    my %seq;
#    my $id = undef;
#    my $s = "";
#    while (my $line = <FILE>) {
#	if ($line =~ /^\>(\S+)/) {
#	    if ($id) {
#		$seq{$id} = $s;
#		$s = "";
#		$id = $1;
#	    }else {
#
#
#		$id = $1;
#	    }
#	    next;
#	}
#	else {
#
#	}
#    }
#}

=head1 METHODS


=head2 get_Id()

	Returns the id of the sequence

	Usage   : my $id = $seq_obj->get_id()
  	Args    : None
  	Returns : String
  	
=cut

sub get_id {
	my ( $self, @args ) = @_;
	return $self->{id} || "";
}
##-------------------------------------------------------------------------##

=head2 set_id()

	Sets the id for the sequence.

	Usage   : $seq_obj->set_id($id)
  	Args    : A string
  	Returns : None
  	
=cut

sub set_id {
	my ( $self, @args ) = @_;
	if ( defined $args[0] ) {
		$self->{id} = $args[0];
	}
}
##-------------------------------------------------------------------------##

=head2 get_seq()

Returns the sequence string

	Usage   : my $seq = $seq->get_seq()
  	Args    : None
  	Returns : The sequence string
  	
=cut

sub get_seq {
	my ( $self, @args ) = @_;
	return $self->{seq} || "";
}
##-------------------------------------------------------------------------##

=head2 set_seq()

Sets the sequence string

	Usage   : $seq = $seq->set_seq("ATGC")
  	Args    : The sequence string
  	Returns : None
  	
=cut

sub set_seq {
	my ( $self, @args ) = @_;
	if ( defined $args[0] ) {
		$self->{seq} = $args[0];
	}
}
##-------------------------------------------------------------------------##

=head2 get_desc()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_desc {
	my ( $self, @args ) = @_;
	return $self->{des} || "";
}
##-------------------------------------------------------------------------##

=head2 set_desc()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub set_desc {
	my ( $self, @args ) = @_;
	if ( defined $args[0] ) {
		$self->{des} = $args[0];
	}
}
##-------------------------------------------------------------------------##

=head2 length()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub length {
	my ( $self, @args ) = @_;

	#if ($self->{length}) {return $self->{length};}
	croak "ERROR: Length of an undefined sequence can't be calculated!"
	  unless $self->{seq};
	return length( $self->{seq} );
}
##-------------------------------------------------------------------------##

=head2 contains_gap()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub contains_gap {
	my $self = shift;
	my $seq  = $self->get_seq();
	if ( $seq =~ /\-/ ) {
		return 1;
	} else {
		return 0;
	}
}
##-------------------------------------------------------------------------##

=head2 get_ungapped_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_ungapped_seq {
	my $self = shift;
	my $seq  = $self->get_seq();
	$seq =~ s/\-//g;
	return $seq;
}
##-------------------------------------------------------------------------##

=head2 contains_term_codon()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub contains_term_codon {
	my $self = shift;
	my $seq  = $self->get_seq();
	if ( $seq =~ /\*/ ) {
		return 1;
	} else {
		return 0;
	}
}
##-------------------------------------------------------------------------##

=head2 get_unterm_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_unterm_seq {
	my $self = shift;
	if ( $self->contains_term_codon ) {
		my $seq = $self->get_seq();
		return $seq =~ s/\*//g;
	} else {
		return $self->get_seq();
	}
}
##-------------------------------------------------------------------------##

=head2 get_cleaned_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_cleaned_seq {
	my $self = shift;
	my $seq  = $self->get_seq();
	$seq =~ s/\*//g;
	$seq =~ s/\-//g;
	return $seq;
}

=head2 get_pretty_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_pretty_seq {
	my ( $self, @args ) = @_;
	my $seq = $self->get_seq();
	my $s   = "";
	my @s   = split( //, $seq );
	for ( my $i = 1 ; $i < @s + 1 ; $i++ ) {

		if ( $i % 60 ) {
			$s .= $s[ $i - 1 ];
		} else {
			$s .= $s[ $i - 1 ];
			$s .= "\n";
		}
	}
	$s .= "\n";
	return $s;
}

=head2 get_fasta()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_fasta {
	my ( $self, @args ) = @_;
	my $s = ">" . $self->get_id() . " " . $self->get_desc() . "\n";
	$s .= $self->get_pretty_seq();
	return $s;
}

sub get_aln_subseq {
	my ( $self, @args ) = @_;
	my ( $start, $end ) = $self->_rearrange( [ "START", "END" ], @args );
	if ( !$start || !$end ) { croak "Start and end is mandetory option\n"; }
	if ( $end < $start ) { croak "Start cannot be more that end\n"; }
	if ( $start < 1 || $end > $self->length ) {
		croak "Invalid coordinates for sequence extraction\n";
	}
	my $seq   = $self->get_seq;
	my @array = split( //, $seq );
	my @s     = @array[ ( $start - 1 ) .. ( $end - 1 ) ];
	return join( "", @s );
}

=head2 get_subseq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_subseq {
	my ( $self, @args ) = @_;
	my ( $start, $end ) = $self->_rearrange( [ "START", "END" ], @args );
	if ( !$start || !$end ) { croak "Start and end is mandetory option\n"; }
	if ( $end < $start ) { croak "Start cannot be more that end\n"; }
	if ( $start < 1 || $end > $self->length ) {
		croak "Invalid coordinates for sequence extraction\n";
	}
	my $seq = $self->get_cleaned_seq;
	my @array = split( //, $seq );
	if ( $start == $end ) {
		return $array[ $start - 1 ];
	} else {
		my @s = @array[ ( $start - 1 ) .. ( $end - 1 ) ];
		return join( "", @s );
	}
}

sub get_gi {
	my ($self) = shift;
	my $id = $self->get_id();
	if ( $id =~ /gi\|(\d+)\|/ || $id =~ /gi\|(\d+)/ ) {
		return $1;
	} else {
		return;
	}
}

sub revcom {
	my $self = shift;
	my $seq  = $self->get_cleaned_seq();
	if ( $seq =~ /[^ATGCNatgcn]/ ) {
		croak "Only a DNA sequence can be reversed and complemented\n";
	}
	my @char          = split //, $seq;
	my $reversed      = reverse(@char);
	my @reversed_char = split( //, $reversed );
	my $r_seq         = "";
	foreach my $i (@reversed_char) {
		if ( $i eq "A" || $i eq "a" ) {
			$r_seq .= "T";
		} elsif ( $i eq "T" || $i eq "t" ) {
			$r_seq .= "A";

		} elsif ( $i eq "G" || $i eq "g" ) {
			$r_seq .= "C";
		} elsif ( $i eq "C" || $i eq "c" ) {
			$r_seq .= "G";
		} elsif ( $i eq "N" || $i eq "n" ) {
			$r_seq .= "N";
		} else {
			croak "Only A,T,G,C,N are allowed in a DNA sequence\n";
		}
	}

	return $r_seq;

}

sub mark_position {
	my ( $self, $required ) = @_;

	my $seq    = $self->get_seq();
	my @s      = split( //, $seq );
	my $return = "";
	for ( my $i = 0 ; $i < @s ; $i++ ) {
		if ( exists $required->{ $i} ) {
				$return .= $required->{$i};
			}
	
			$return .= $s[$i];
	
		}
		my $newseq
		  = SeqToolBox::Seq->new( -id => $self->get_id(), -seq => $return );
		return $newseq;

}

=head1 SEE ALSO

<l>SeqToolBox::SeqDB</l>

=head1 COPYRIGHTS

Copyright (c) 2003 by Malay < malay@bioinformatics.org >. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

1;
