# $Id: SeqDB.pm 663 2011-10-18 22:55:03Z malay $
# Perl module for SeqToolBox::SeqDB
# Author: Malay <malay@bioinformatics.org>
# Copyright (c) 2006, 2007 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

SeqToolBox::SeqDB  - DESCRIPTION of Object

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

package SeqToolBox::SeqDB;
@ISA = qw(SeqToolBox::Root);
use SeqToolBox::Root;
use SeqToolBox::Seq;
use IO::File;
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
	my $self  = {};
	bless $self, ref($class) || $class;
	$self->_init(@_);

	return $self;
}

# _init is where the heavy stuff will happen when new is called

sub _init {
	my ( $self, @args ) = @_;
	my ( $fh, $file, $format, $index )
		= $self->_rearrange( [ "FH", "FILE", "FORMAT", "INDEX" ], @args );

	if ( !$format ) {
		$self->{'format'} = 'FASTA';
	} elsif ( $format !~ /^FASTA$/i ) {
		croak('Only FASTA format is supported');

	} else {

	}

	if ( $fh && $fh->isa("IO::File") ) {
		$self->{fh} = $fh;
	} elsif ($fh) {
		croak "File Handle must be a IO::File object\n";
	} elsif ($file) {
		$self->{file} = $file;
		$self->_open_file($file);
	} else {

		#croak "Either a IO::File object or a filename has to be given\n";
	}

	if ($index) {
		$self->{file_index} = $index;
	}

	#open (my $f , $file) || croak "Can't open $file\n";
	#$self->{fh} = $f;
	if ( exists $self->{fh} ) {
		$self->_parse_file();
	}
	$self->{_counter} = 0;

	return $self;
}

##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##

=head1 PUBLIC METHODS

=cut

##-------------------------------------------------------------------------##

=head2 has_next()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub has_next {
	my $self = shift;

	#	my $number = $#{$self->{seq}};
	#	if ( $self->{_counter} <= $number ) {
	#		return 1;
	#	} else {
	#		return undef;
	#	}
	#
	if ( $self->_in_memory ) {

		if ( $self->{_counter} < $self->seq_count ) {
			return 1;
		} else {
			return;
		}
	}

	if ( defined( $self->{next_id} ) ) {
		return 1;
	} else {
		return;
	}
}

=head2 reset_iterator()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub reset_iterator {
	my ( $self, @args ) = @_;

	#$self->{fh}->close();
	#$self->_open_file($self->{file});

	#open (my $f , $file) || croak "Can't open $file\n";
	#$self->{fh} = $f;
	unless ( $self->_in_memory() ) {
		my $fh = $self->{fh};
		$fh->seek( 0, 0 );
		$self->_parse_file();
	}

	#	$self->{_counter} = 0;

	#return $self;

}
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 next_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub next_seq {
	my $self = shift;

	if ( !exists $self->{fh} && exists $self->{seq} ) {
		my $present = $self->{_counter};

		if ( defined $self->{seq}->[$present] ) {
			$self->{_counter}++;
			return $self->{seq}->[$present];
		} else {
			$self->{_counter} = 0;
			return;
		}

	}

	if ( $self->{next_id} ) {

		$self->{current_id}  = $self->{next_id};
		$self->{current_des} = $self->{next_des};
		$self->{next_id}     = undef;
		$self->{next_des}    = undef;
		my $fh  = $self->{fh};
		my $seq = "";

		while ( my $line = <$fh> ) {
			chomp $line;
			next unless $line;

			#print STDERR $line,"\n";
			if ( $line =~ /^>(\S+)/ ) {
				$self->{next_id} = $1;

				if ( $line =~ /^>(\S+)\s+(.*)/ ) {
					$self->{next_des} = $2;
				} else {
					$self->{next_des} = "";
				}

				#print STDERR $self->{current_id},"\n";
				return
					SeqToolBox::Seq->new( -id  => $self->{current_id},
										  -des => $self->{current_des},
										  -seq => $seq
					);
			} else {
				$line =~ s/\s+//g;
				$seq .= $line;
			}
		}

		return
			SeqToolBox::Seq->new( -id  => $self->{current_id},
								  -des => $self->{current_des},
								  -seq => $seq
			);

	} else {
		$self->reset_iterator();
		return;
	}
}
##-------------------------------------------------------------------------##

=head2 seq_count()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub seq_count {
	my $self = shift;

	if ( defined $self->{seq_count} ) {
		return $self->{seq_count};
	}

	if ( !exists $self->{fh} && exists $self->{seq} ) {
		$self->{seq_count} = scalar( @{ $self->{seq} } );
		return $self->{seq_count};
	}

	my $fh  = $self->{fh};
	my $pos = $fh->getpos;

	my $count = 0;
	$fh->seek( 0, 0 );

	while ( my $line = <$fh> ) {
		if ( $line =~ /^\>/ ) {
			$count++;
		}
	}
	$fh->setpos($pos);
	return $count;

	#return $#{ $self->{seq} } ? ( $#{ $self->{seq} } + 1 ) : 0;
}

##-------------------------------------------------------------------------##

=head2 get_seq_by_index()

Get the sequence using sequence number. The first sequence is 0.

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_seq_by_index {
	my ( $self, $index ) = @_;

	if ( $self->_in_memory ) {
		if ( defined $self->{seq}->[$index] ) {
			return $self->{seq}->[$index];
		} else {
			croak "The seq with $index not found\n";
		}
	}

	#print STDERR "seqdb:", $index,"\n";

	#	if (defined $self->{seq}->[$index]) {
	#		return $self->{seq}->[$index];
	#	}else {
	#		return undef;
	#	}
	my $fh          = $self->{fh};
	my $current_pos = $fh->getpos();
	$fh->seek( 0, 0 );

	my $last_seq = "";
	my $last_id  = "";
	my $last_des = "";
	my $count    = 0;

	#my $fh = $self->{fh};
	while ( my $line = <$fh> ) {
		chomp $line;

		#print STDERR "$line\n";
		if ( $line =~ /^>(\S+)/ ) {

			#print STDERR "defile found\n";
			$count++;

			if ($last_id) {
				if ( $index == ( $count - 2 ) ) {
					$fh->setpos($current_pos);
					return
						SeqToolBox::Seq->new( -id  => $last_id,
											  -des => $last_des,
											  -seq => $last_seq
						);
				}
				$last_id  = $1;
				$last_seq = "";

				if ( $line =~ /^>(\S+)\s+(.*)/ ) {
					$last_des = $2;
				} else {
					$last_des = "";
				}
			} else {
				$last_id  = $1;
				$last_seq = "";

				if ( $line =~ /^>(\S+)\s+(.*)/ ) {
					$last_des = $2;
				} else {
					$last_des = "";
				}
			}

		} else {
			$last_seq .= $line;

		}

	}

	if ( $index == ( $count - 1 ) ) {
		$fh->setpos($current_pos);
		return
			SeqToolBox::Seq->new( -id  => $last_id,
								  -des => $last_des,
								  -seq => $last_seq
			);
	} else {
		return undef;
	}
}

=head1 PRIVATE METHODS

=cut

##-------------------------------------------------------------------------##

=head2 _open_file()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub _open_file {
	my ( $self, @args ) = @_;
	my $fh = IO::File->new( $args[0], "r" );

	if ( defined $fh ) {

		#open (FILE, $args[0]) || croak "Can't open $args[0]\n";
		$self->{fh} = $fh;
	} else {
		croak "Could not open ", $args[0], "\n";
	}
}

##-------------------------------------------------------------------------##

=head2 _parse_file()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub _parse_file {
	my ( $self, @args ) = @_;

	#print STDERR "parse_file called\n";
	#open (FILE, $args[0]) || croak "Can't open $args[0]\n";
	my $last_seq = "";
	my $last_id  = "";
	my $last_des = "";
	my $fh       = $self->{fh};

	while ( my $line = <$fh> ) {
		chomp $line;

		#print STDERR "$line\n";
		if ( $line =~ /^>(\S+)/ ) {

#print STDERR "defile found\n";
#if ($last_id) {
#$self->add_seq(SeqToolBox::Seq->new(-id=>$last_id,-des=>$last_des,-seq=>$last_seq));

			$last_id = $1;

			#$last_seq ="";
			if ( $line =~ /^>(\S+)\s+(.*)/ ) {
				$last_des = $2;

			} else {
				$last_des = "";
			}

			#}else{
			#	$last_id = $1;
			#	$last_seq = "";
			#	if ($line =~ /^>(\S+)\s+(.*)/ ){
			#	$last_des = $2;
			#}else{
			#	$last_des = "";
			#}
			#}

			last;
		}    #else{
		     #$last_seq .= $line;

		#}
	}

#$self->add_seq(SeqToolBox::Seq->new(-id=>$last_id,-des=>$last_des,-seq=>$last_seq));
#close (FILE);
	$self->{next_id}  = $last_id;
	$self->{next_des} = $last_des;
}

##-------------------------------------------------------------------------##

=head2 add_seq()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub add_seq {
	my ( $self, $seq ) = @_;

	#print STDERR "Add_seq called\n";
	if ( exists $self->{seq} ) {

		#my @array = @{$self->{seq}};
		push @{ $self->{seq} }, $seq;

		#$self->{seq} = \@array;
	} else {
		my @array;
		push @array, $seq;
		$self->{seq} = \@array;
	}
}

sub get_longest_seq {
	my $self           = shift;
	my $longest_index  = 0;
	my $longest_length = 0;

	unless ( $self->{fh} ) {
		croak "File handle not defined\n";
	}
	my $current_pos = $self->{fh}->getpos();
	$self->reset_iterator();
	my $counter = 0;

	while ( my $seq = $self->next_seq() ) {
		$counter++;
		my $len = $seq->length();

		if ( $longest_length < $len ) {
			$longest_length = $len;
			$longest_index  = $counter;
		}
	}
	$self->{fh}->setpos($current_pos);
	return $self->get_seq_by_index( $longest_index - 1 );

}
##-------------------------------------------------------------------------##

=head2 _in_memory()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub _in_memory {
	my ( $self, @args ) = @_;

	if ( !exists $self->{fh} && exists $self->{seq} ) {

		#$self->{seq_count} = scalar(@{$self->{seq}});
		return 1;
	} else {
		return;
	}
}

=head2 create_index()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub create_index {
	my ( $self, $index_file ) = @_;

	# Indexing an in memory file does not make any sense
	if ( $self->_in_memory ) {
		croak "Indexing an in memory SeqDB object is plain wrong\n";
	}

	#indexing can be in memory hash object or a file to use later.
	my %index;

	if ($index_file) {
		require SDBM_File || die "Can't load SDBM_module\n";
		tie( %index, 'SDBM_File', $index_file, O_RDWR | O_CREAT, 0666 )
			|| die "Can't create tied hash\n";
	}

	my $fh = $self->{fh};

	# the file pointer can be at any location. We will first get the current
	# location and then restore it later

	my $current_pos = $fh->getpos();
	$fh->seek( 0, 0 );

	#	my $last_seq = "";
	#	my $last_id  = "";
	#	my $last_des = "";
	#	my $count    = 0;

	#my $fh = $self->{fh};

# Length function does not tell the byte location. To get around this problem we will use
# a lexically scope use bytes statement
	{
		use bytes;

		while ( my $line = <$fh> ) {

			#chomp $line;

			#print STDERR "$line\n";
			if ( $line =~ /^>(\S+)/ ) {

				#print STDERR "defile found\n";
				#			$count++;
				my $id = $1;
#				print STDERR "Position: ", tell($fh), "\n";
				my $pos = ( tell($fh) - length($line) );
				if (exists $index{$id}) {
					croak "Duplicate ID found $id\n";
				}else {
				$index{$id} = $pos;
				}
			}
		}
	}

	$fh->setpos($current_pos);

	if ( !$index_file ) {
		$self->{seq_index} = \%index;
	}

}

=head2 get_seq_by_id()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub get_seq_by_id {
	my ( $self, $id ) = @_;
	my %index;

	my $fh              = $self->{fh};
	my $current_pos     = $fh->getpos();
	my $index_file_name;
	if (exists $self->{file}) {
	$index_file_name = $self->{file} . '.idx';
	}

	#first check for presense of in-memory index;
	if ( exists $self->{seq_index} ) {
		%index = %{ $self->{seq_index} };

	}

	# then check whether an index file has been supplied
	elsif ( exists $self->{file_index} ) {
		require SDBM_File || die "Can't load SDBM_File\n";
		tie( %index, 'SDBM_File', $self->{file_index}, O_RDONLY, 0666 )
			|| die "Can't read index file\n";
	}

	# See whether we can read the index file from disk
	elsif ( exists $self->{file} &&  -s $index_file_name . '.pag' ) {
#		print STDERR "Id $id called\n";
#		print STDERR "Index file name: $index_file_name\n";
		require SDBM_File
			|| die "Can't load SDBM_File\n";
		tie( %index, 'SDBM_File', $index_file_name, O_RDONLY, 0666 )
			|| die "Can't read index file\n";
#		print STDERR "After tie\n";
	} else {
#		print STDERR "Manually seeking\n";
		$fh->seek( 0, 0 );
		$self->_parse_file();

		while ( my $seq = $self->next_seq ) {
			if ( $seq->get_id eq $id ) {
				$fh->setpos($current_pos);
				return $seq;
			}
		}
	}

	if ( %index && exists $index{$id} ) {
#		print STDERR $index{$id}, "\n";
		$fh->seek( $index{$id}, 0 ) || die "Can't set to required position";
		$self->_parse_file();

		if ( my $seq = $self->next_seq ) {
			if ( $seq->get_id eq $id ) {
				$fh->setpos($current_pos);
				return $seq;
			} else {
				croak "Seq id in the index does not match\n";
			}
		}

	} else {
		return;
	}

}

=head2 DESTROY()

Describe your function here

	Usage   :
  	Args    :
  	Returns : 
  	
=cut

sub DESTROY {
	my ( $self, @args ) = @_;

	#close($self->{fh});
	if ( $self->{file} ) {
		$self->{fh}->close;
	}
}

=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2006 by Malay <malay@bioinformatics.org>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

=head1 APPENDIX

=cut

1;
