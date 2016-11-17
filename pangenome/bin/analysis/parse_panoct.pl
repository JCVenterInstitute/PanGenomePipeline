#!/usr/bin/env perl
# $Id: parse_panoct.pl 38648 2013-06-04 16:56:48Z ekelsey $

##---------------------------------------------------------------------------##
##  File: parse_panoct.pl
##       
##  Author:
##        Malay <malay@bioinformatics.org>
##
##  Description:
##     
#******************************************************************************
#* Copyright (C) 2010 Malay K Basu <malay@bioinformatics.org> 
#* This work is distributed under the license of Perl iteself.
###############################################################################

=head1 NAME

parse_panoct.pl - One line description.

=head1 SYNOPSIS

parse_panoct.pl [options] -o <option>


=head1 DESCRIPTION

Write a description of your prgram. 


=head1 ARGUMENTS 

=over 4

=item B<--option|-o>

First option.



=back

=head1 OPTIONS

Something here.


=head1 SEE ALSO

=head1 COPYRIGHT

Copyright (c) 2010 Malay K Basu <malay@bioinformatics.org>

=head1 AUTHORS

Malay K Basu <malay@bioinformatics.org>

=cut



##---------------------------------------------------------------------------##
## Module dependencies
##---------------------------------------------------------------------------##

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;


my $file = shift;
open ( my $fh, $file ) || die "Can't open $file\n";
my $count = 0;

while ( my $line = <$fh> ) {

    chomp $line;
    my @f = split ( /\t/, $line );
    my @real_f;
    
    foreach my $i (@f) {  
	
        next if $i =~ /^\-{2,}/;

        push @real_f , $i;

    }

    print join ("\t", @real_f), "\n";

}
