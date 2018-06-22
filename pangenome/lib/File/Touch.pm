package File::Touch;
$File::Touch::VERSION = '0.09';
use 5.006;
use warnings;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(touch);

use Carp;
use IO::File;
use File::stat;
use Fcntl;

my $SYSOPEN_MODE = O_WRONLY|O_CREAT;
eval {
    $SYSOPEN_MODE |= O_NONBLOCK;
};
if($@) {
    # OK, we don't have O_NONBLOCK:
    # probably running on Windows.
}
eval {
    $SYSOPEN_MODE |= O_NOCTTY;
};
if($@) {
    # OK, we don't have O_NOCTTY:
    # probably running on Windows.
}

sub new
{
    my ($caller, %arg) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = bless{}, $class;

    my $atime_only  = $arg{atime_only} || 0; # If nonzero, change only the access time.
    my $mtime_only  = $arg{mtime_only} || 0; # If nonzero, change only the modification time.
    my $no_create   = $arg{no_create}  || 0; # If nonzero, don't create if not already there.
    my $reference   = $arg{reference};       # If defined, use this file's times instead of current time.
    my $time        = $arg{time};            # If defined, use this time instead of current time.
    my $atime       = $arg{atime};           # If defined, use this time for access time instead of current time.
    my $mtime       = $arg{mtime};           # If defined, use this time for modification time instead of current time.

    if ($atime_only && $mtime_only){
        croak("Incorrect usage: 'atime_only' and 'mtime_only' are both set - they are mutually exclusive.");
    }

    if (defined $time) {
        if ((defined $atime) || (defined $mtime)) {
            croak("Incorrect usage: 'time' should not be used with either ",
                  "'atime' or 'mtime' - ambiguous.");
        }
        $atime = $time unless $mtime_only;
        $mtime = $time unless $atime_only;
    }

    if (defined $reference) {
        if ((defined $time) || (defined $atime) || (defined $mtime)) {
            croak("Incorrect usage: 'reference' should not be used with 'time', 'atime' or 'mtime' - ambiguous.");
        }
        if (-e $reference) {
            my $sb = stat($reference) or croak("Could not stat ($reference): $!");
            $atime = $sb->atime unless $mtime_only;
            $mtime = $sb->mtime unless $atime_only;
        }
        else {
            croak("Reference file ($reference) does not exist");
        }
    }

    $self->{_atime}      = $atime;
    $self->{_mtime}      = $mtime;
    $self->{_no_create}  = $no_create;
    $self->{_atime_only} = $atime_only;
    $self->{_mtime_only} = $mtime_only;

    return $self;
}

sub touch
{
    my ($caller, @files) = @_;
    my $caller_is_obj = ref($caller);
    my $self;

    if ($caller_is_obj){
        $self = $caller;
    }
    else {
        unshift @files, $caller;
        $self->{_atime}      = undef;
        $self->{_mtime}      = undef;
        $self->{_no_create}  = 0;
        $self->{_atime_only} = 0;
        $self->{_mtime_only} = 0;
    }

    my $count = 0;

    foreach my $file (@files) {
        my $time = time();
        my ($atime,$mtime);
        
        if (-e $file) {
            my $sb = stat($file) or croak("Could not stat ($file): $!");
            $atime = $sb->atime;
            $mtime = $sb->mtime;
        }
        else {
            unless ($self->{_no_create}) {
                sysopen my $fh,$file,$SYSOPEN_MODE or croak("Can't create $file : $!");
                close $fh or croak("Can't close $file : $!");
                $atime = $time;
                $mtime = $time;
            }
        }
        unless ($self->{_mtime_only}) {
            $atime = $time;
            $atime = $self->{_atime} if (defined $self->{_atime});
        }
        unless ($self->{_atime_only}) {
            $mtime = $time;
            $mtime = $self->{_mtime} if (defined $self->{_mtime});
        }
        if (utime($atime,$mtime,$file)) {
            $count++;
        }
    }
    return $count;
}

1;

__END__

=head1 NAME

File::Touch - update file access and modification times, optionally creating files if needed

=head1 SYNOPSIS

 use File::Touch;
 @file_list = ('one.txt','../two.doc');
 $count = touch(@file_list);

 use File::Touch;
 $reference_file = '/etc/passwd';
 $touch_obj = File::Touch->new(
                  reference => $reference_file,
                  no_create => 1
              );
 @file_list = ('one.txt','../two.doc');
 $count     = $touch_obj->touch(@file_list);

=head1 DESCRIPTION

Here's a list of arguments that can be used with the object-oriented contruction:

=over 4

=item atime_only => [0|1]

If nonzero, change only the access time of files. Default is zero.

=item mtime_only => [0|1]

If nonzero, change only the modification time of files. Default is zero.

=item no_create => [0|1]

If nonzero, do not create new files. Default is zero.

=item reference => $reference_file

If defined, use timestamps from this file instead of current time. Default is undefined.

=item time => $time

If defined, then this value will be used for both access time and modification time,
whichever of those are set. This time is overridden by the C<atime> and C<mtime> arguments,
if you use them.

=item atime => $time

If defined, use this time (in epoch seconds) instead of current time for access time.

=item mtime => $time

If defined, use this time (in epoch seconds) instead of current time for modification time.

=back

=head1 Examples

=head2 Update access and modification times, creating nonexistent files

 use File::Touch;
 my @files = ('one','two','three');
 my $count = touch(@files);
 print "$count files updated\n";

=head2 Set access time forward, leave modification time unchanged

 use File::Touch;
 my @files = ('one','two','three');
 my $day = 24*60*60;
 my $time = time() + 30 * $day;
 my $ref = File::Touch->new( atime_only => 1, time => $time );
 my $count = $ref->touch(@files);
 print "$count files updated\n";

=head2 Set modification time back, update access time, do not create nonexistent files

 use File::Touch;
 my @files = ('one','two','three');
 my $day = 24*60*60;
 my $time = time() - 30 * $day;
 my $ref = File::Touch->new( mtime => $time, no_create => 1 );
 my $count = $ref->touch(@files);
 print "$count files updated\n";

=head1 AUTHOR

Nigel Wetters Gourlay (nwetters@cpan.org)

=head1 COPYRIGHT

Copyright (c) 2001,2007,2009 Nigel Wetters Gourlay. All Rights Reserved.
This module is free software. It may be used, redistributed
and/or modified under the same terms as Perl itself.

