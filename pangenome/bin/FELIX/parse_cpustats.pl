#!/usr/bin/env perl
#Copyright (C) 2017-2022 The J. Craig Venter Institute (JCVI).  All rights reserved #This program is free software: you can redistribute it and/or modify #it under the terms of the GNU General Public License as published by #the Free Software Foundation, either version 3 of the License, or #(at your option) any later version.

#This program is distributed in the hope that it will be useful, #but WITHOUT ANY WARRANTY; without even the implied warranty of #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the #GNU General Public License for more details.

#You should have received a copy of the GNU General Public License #along with this program.  If not, see <http://www.gnu.org/licenses/>.

use FileHandle;
use Getopt::Long;
use Carp;
use strict;
use warnings;
use List::Util qw[min max];

my @cpus = (); # cpu usage
my @rams = (); # ram usage
my $total_cpu = 0;
my $max_cpu = 0;
my $mean_cpu;
my $median_cpu;
my $total_ram = 0;
my $max_ram = 0;
my $mean_ram;
my $median_ram;

GetOptions('help' => \my $help,
	'debug' => \my $debug,
	'list=s' => \my $list);

if ($help) {
   system("clear");
   print STDERR <<_EOB_;
GetOptions('help' => help,
	'debug' => debug,
	'list=s' => list);
_EOB_
    exit(0);
}
print "Sample\tRAM(Mbytes)\tCPU(seconds)\n";
# read file which specifies the output name for each sample, and the cpu_stats file for each sample
my $count = 0;
open (my $infile, "<", $list) || die ("ERROR: cannot open file $list\n");
while (my $line = <$infile>)  {
    chomp $line;
    (my $samplename, my $filename) = split(/\t/, $line);  # split the scalar $line on tab

    # read in cpustats
    my $cpufile;
    unless (open ($cpufile, "<", $filename) )  {
	die ("cannot open cpustats file: $filename!\n");
    }
    my $cpu;
    my $user;
    my $system;
    my $ram;
    while (my $line2 = <$cpufile>) {
	chomp($line);
	if ($line2 =~ /\s+User time \(seconds\): (\d+\.\d+)/) {
	    $user = $1;
	} elsif ($line2 =~ /\s+System time \(seconds\): (\d+\.\d+)/) {
	    $system = $1;
	} elsif ($line2 =~ /\s+Maximum resident set size \(kbytes\): (\d+)/) {
	    $ram = $1;
	}
    }
    $cpu = $user + $system;
    $cpus[$count] = $cpu;
    if ($cpu > $max_cpu) {
	$max_cpu = $cpu;
    }
    $total_cpu += $cpu;
    $ram /= 1000;
    $rams[$count] = $ram;
    if ($ram > $max_ram) {
	$max_ram = $ram;
    }
    $total_ram += $ram;
    print "$samplename\t$ram\t$cpu\n";
    close($cpufile);
    $count++;
}
close($infile);
$mean_cpu = $total_cpu / $count;
@cpus = sort {$a <=> $b} @cpus;
$median_cpu = ($count % 2) ? $cpus[($count / 2)] : (($cpus[(($count / 2) - 1)] + $cpus[($count / 2)]) / 2);
$mean_ram = $total_ram / $count;
@rams = sort {$a <=> $b} @rams;
$median_ram = ($count % 2) ? $rams[($count / 2)] : (($rams[(($count / 2) - 1)] + $rams[($count / 2)]) / 2);
print "Mean\t$mean_ram\t$mean_cpu\n";
print "Median\t$median_ram\t$median_cpu\n";
print "Total\t(Max) $max_ram\t$total_cpu\n";
print "Max\t$max_ram\t$max_cpu\n";

exit(0);
