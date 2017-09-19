#!/usr/bin/env perl

###############################################################################
#                                                                             #
#       Copyright (C) 2010 J. Craig Venter Institute (JCVI).                  #
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
