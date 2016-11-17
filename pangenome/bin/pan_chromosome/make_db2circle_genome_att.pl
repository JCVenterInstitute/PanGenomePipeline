#!/usr/bin/env perl
#
#Copy (C) 2016 The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

# created  April 08, 2014
# last modified 03/02/2016
# script to create genome.att, config.file and data.file files for db2circle (flat file independent version) using:
# 1) pan-genome core genes "assembly" Core.attfGI file (created by gene_order.pl)
# 2) panOCT shared_clusters.txt file
# 3) pan-genome role_id_lookup.txt file (optional now)

my $prog = $0;
$prog =~ s/.*\///;
my $invocation = $prog . " @ARGV";

use strict;
use Cwd;
use Getopt::Std;
getopts ('a:b:c:dhr:V');

############## Declare global variables #################
our ($opt_a,$opt_b,$opt_c,$opt_d,$opt_h,$opt_r,$opt_V);
my ($basedir, $att_file, $clusters_file, $roleid_file, $cluster, $end5, $end3, $role_id, $main_role, $locus, $com_name, $feat_type, $length, $prefix, $inserts, $genlen, $DEBUG);
my $cwd = &getcwd();
my @a = ();
my @b = ();
my @c = ();
my @d = ();
my @roles = ();
my $line = "";
my %role_id_hash = (); # key = role_id, value = main role
my %cluster2id_hash = (); # key = cluster_id, value = role_id
my $version = "1.0";

############## parse options ###############
if ($opt_h) { &option_help; }
if ($opt_V) {die "$version\n";}
if ($opt_b) {$basedir = $opt_b;}  else { $basedir = $cwd; }
if ($opt_a) {$att_file = $opt_a;} else { &option_help; }
if ($opt_c) {$clusters_file = $opt_c;} else { &option_help; }
if ($opt_r) {$roleid_file = $opt_r;} else { undef $roleid_file; }
if ($opt_d) {$DEBUG = 1;} else {$DEBUG = 0;} # Debug mode.

###############
# Subroutines #
###############
sub option_help {
  print <<_EOB_;
    $prog - To create a genome.att file for db2circle (flat file independent version)
      Usage: $prog <options>
    Switch: -h for help\
  Option:
    -b: base directory path [default = PWD]
    -a: att file from pan-genome core genes assembly with gaps for largest fGIs (Core.attfGI [required]
    -c: shared_clusters.txt file from PanOCT [required]
    -r: role_id_lookup.txt file from Pan-genome pipeline [optional; default embedded]
    -V: print version information
    -d: DEBUG MODE (default = 0)
  Format: 
  Output:
_EOB_
  exit;
}

sub print_config_file {
  my ($path,$len,$DEBUG) = @_;
  print "Pan-genome length is <$len> bp.\n" if ($DEBUG);
  open (CONFIG, ">$basedir/config.file") || die "can't open file $basedir/config.file\n";
  print CONFIG "LENGTH:\t$len\n";
  print CONFIG "TOP_BORDER\t1\tFORWARD\tTOP_BORDER\n";
  print CONFIG "BKGRD_L\t1\tFORWARD\tBCKGRND_ISLAND\tGRAY88\n";
  print CONFIG "L_FGI_50\t1\tBOTH\tISLAND\tRED\n";
  print CONFIG "L_FGI_50_1-3\t1\tBOTH\tQUARTER_ISLAND\tRED\n";
  print CONFIG "L_FGI_50_4-6\t1\tBOTH\tHALF_ISLAND\tRED\n";
  print CONFIG "L_FGI_50_7-9\t1\tBOTH\tTHREEQUARTER_ISLAND\tRED\n";
  print CONFIG "L_FGI_50_10+\t1\tBOTH\tISLAND\tRED\n";
  print CONFIG "L_FGI_25\t1\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "L_FGI_25_1-3\t1\tBOTH\tQUARTER_ISLAND\tORANGE\n";
  print CONFIG "L_FGI_25_4-6\t1\tBOTH\tHALF_ISLAND\tORANGE\n";
  print CONFIG "L_FGI_25_7-9\t1\tBOTH\tTHREEQUARTER_ISLAND\tORANGE\n";
  print CONFIG "L_FGI_25_10+\t1\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "L_FGI_10\t1\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "L_FGI_10_1-3\t1\tBOTH\tQUARTER_ISLAND\tGREEN\n";
  print CONFIG "L_FGI_10_4-6\t1\tBOTH\tHALF_ISLAND\tGREEN\n";
  print CONFIG "L_FGI_10_7-9\t1\tBOTH\tTHREEQUARTER_ISLAND\tGREEN\n";
  print CONFIG "L_FGI_10_10+\t1\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "L_FGI_1\t1\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "L_FGI_1_1-3\t1\tBOTH\tQUARTER_ISLAND\tBLUE\n";
  print CONFIG "L_FGI_1_4-6\t1\tBOTH\tHALF_ISLAND\tBLUE\n";
  print CONFIG "L_FGI_1_7-9\t1\tBOTH\tTHREEQUARTER_ISLAND\tBLUE\n";
  print CONFIG "L_FGI_1_10+\t1\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "BKGRD_M\t2\tFORWARD\tBCKGRND_ISLAND\tGRAY88\n";
  print CONFIG "M_FGI_50\t2\tBOTH\tISLAND\tRED\n";
  print CONFIG "M_FGI_50_1-3\t2\tBOTH\tQUARTER_ISLAND\tRED\n";
  print CONFIG "M_FGI_50_4-6\t2\tBOTH\tHALF_ISLAND\tRED\n";
  print CONFIG "M_FGI_50_7-9\t2\tBOTH\tTHREEQUARTER_ISLAND\tRED\n";
  print CONFIG "M_FGI_50_10+\t2\tBOTH\tISLAND\tRED\n";
  print CONFIG "M_FGI_25\t2\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "M_FGI_25_1-3\t2\tBOTH\tQUARTER_ISLAND\tORANGE\n";
  print CONFIG "M_FGI_25_4-6\t2\tBOTH\tHALF_ISLAND\tORANGE\n";
  print CONFIG "M_FGI_25_7-9\t2\tBOTH\tTHREEQUARTER_ISLAND\tORANGE\n";
  print CONFIG "M_FGI_25_10+\t2\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "M_FGI_10\t2\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "M_FGI_10_1-3\t2\tBOTH\tQUARTER_ISLAND\tGREEN\n";
  print CONFIG "M_FGI_10_4-6\t2\tBOTH\tHALF_ISLAND\tGREEN\n";
  print CONFIG "M_FGI_10_7-9\t2\tBOTH\tTHREEQUARTER_ISLAND\tGREEN\n";
  print CONFIG "M_FGI_10_10+\t2\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "M_FGI_1\t2\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "M_FGI_1_1-3\t2\tBOTH\tQUARTER_ISLAND\tBLUE\n";
  print CONFIG "M_FGI_1_4-6\t2\tBOTH\tHALF_ISLAND\tBLUE\n";
  print CONFIG "M_FGI_1_7-9\t2\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "M_FGI_1_10+\t2\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "BKGRD_S\t3\tFORWARD\tBCKGRND_ISLAND\tGRAY88\n";
  print CONFIG "S_FGI_50\t3\tBOTH\tISLAND\tRED\n";
  print CONFIG "S_FGI_50_1-3\t3\tBOTH\tQUARTER_ISLAND\tRED\n";
  print CONFIG "S_FGI_50_4-6\t3\tBOTH\tHALF_ISLAND\tRED\n";
  print CONFIG "S_FGI_50_7-9\t3\tBOTH\tTHREEQUARTER_ISLAND\tRED\n";
  print CONFIG "S_FGI_50_10+\t3\tBOTH\tISLAND\tRED\n";
  print CONFIG "S_FGI_25\t3\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "S_FGI_25_1-3\t3\tBOTH\tQUARTER_ISLAND\tORANGE\n";
  print CONFIG "S_FGI_25_4-6\t3\tBOTH\tHALF_ISLAND\tORANGE\n";
  print CONFIG "S_FGI_25_7-9\t3\tBOTH\tTHREEQUARTER_ISLAND\tORANGE\n";
  print CONFIG "S_FGI_25_10+\t3\tBOTH\tISLAND\tORANGE\n";
  print CONFIG "S_FGI_10\t3\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "S_FGI_10_1-3\t3\tBOTH\tQUARTER_ISLAND\tGREEN\n";
  print CONFIG "S_FGI_10_4-6\t3\tBOTH\tHALF_ISLAND\tGREEN\n";
  print CONFIG "S_FGI_10_7-9\t3\tBOTH\tTHREEQUARTER_ISLAND\tGREEN\n";
  print CONFIG "S_FGI_10_10+\t3\tBOTH\tISLAND\tGREEN\n";
  print CONFIG "S_FGI_1\t3\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "S_FGI_1_1-3\t3\tBOTH\tQUARTER_ISLAND\tBLUE\n";
  print CONFIG "S_FGI_1_4-6\t3\tBOTH\tHALF_ISLAND\tBLUE\n";
  print CONFIG "S_FGI_1_7-9\t3\tBOTH\tTHREEQUARTER_ISLAND\tBLUE\n";
  print CONFIG "S_FGI_1_10+\t3\tBOTH\tISLAND\tBLUE\n";
  print CONFIG "ORF\t4\tFORWARD\tISLAND\n";
  print CONFIG "ORF\t5\tREVERSE\tISLAND\n";
  print CONFIG "BOTTOM_BORDER\t5\tFORWARD\tBOTTOM_BORDER\n";
  close (CONFIG);
}
###############
#   M A I N   #
###############
  # read in and store role_id to main_role lookup info
  if ($roleid_file) {
    print "Using user-supplied role_id lookup table\n" if ($DEBUG);
    open (ROLEIDFILE, "<$basedir/$roleid_file");
    while (<ROLEIDFILE>)  {
      chomp;
      @a=split(/\t/);
      print "ROLE_ID2MAIN_ROLE: <$a[0]>\t<$a[1]>\n" if ($DEBUG);
      $role_id_hash{$a[0]} = $a[1];
    }
    close (ROLEIDFILE);
  }
  else { # These are the JCVI/TIGR hard-coded main role categories, which are default
    print "Using built-in JCVI/TIGR main role ids\n" if ($DEBUG);
    %role_id_hash = (0 => 'NONE',
		     46 => 'metabolism',
		     69 => 'Amino acid biosynthesis',
		     70 => 'Amino acid biosynthesis',
		     71 => 'Amino acid biosynthesis',
		     73 => 'Amino acid biosynthesis',
		     74 => 'Amino acid biosynthesis',
		     75 => 'Amino acid biosynthesis',
		     76 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     77 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     78 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     79 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     80 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     81 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     82 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     83 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     84 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     85 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     86 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     88 => 'Cell envelope',
		     89 => 'Cell envelope',
		     90 => 'Cell envelope',
		     91 => 'Cell envelope',
		     92 => 'Cellular processes',
		     93 => 'Cellular processes',
		     94 => 'Cellular processes',
		     95 => 'Protein fate',
		     96 => 'Cellular processes',
		     97 => 'Protein fate',
		     98 => 'Cellular processes',
		     100 => 'Central intermediary metabolism',
		     102 => 'Central intermediary metabolism',
		     103 => 'Central intermediary metabolism',
		     104 => 'Central intermediary metabolism',
		     105 => 'Energy metabolism',
		     106 => 'Central intermediary metabolism',
		     108 => 'Energy metabolism',
		     109 => 'Energy metabolism',
		     110 => 'Energy metabolism',
		     111 => 'Energy metabolism',
		     112 => 'Energy metabolism',
		     113 => 'Energy metabolism',
		     114 => 'Energy metabolism',
		     116 => 'Energy metabolism',
		     117 => 'Energy metabolism',
		     118 => 'Energy metabolism',
		     119 => 'Energy metabolism',
		     120 => 'Energy metabolism',
		     121 => 'Fatty acid and phospholipid metabolism',
		     122 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     123 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     124 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     125 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     126 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     127 => 'Purines, pyrimidines, nucleosides, and nucleotides',
		     129 => 'Regulatory functions',
		     130 => 'DNA metabolism',
		     131 => 'DNA metabolism',
		     132 => 'DNA metabolism',
		     133 => 'Transcription',
		     134 => 'Transcription',
		     135 => 'Transcription',
		     136 => 'Protein synthesis',
		     137 => 'Protein synthesis',
		     138 => 'Protein fate',
		     140 => 'Protein fate',
		     141 => 'Transport and binding proteins',
		     142 => 'Transport and binding proteins',
		     143 => 'Transport and binding proteins',
		     144 => 'Transport and binding proteins',
		     145 => 'Transport and binding proteins',
		     146 => 'Transport and binding proteins',
		     147 => 'Transport and binding proteins',
		     149 => 'Cellular processes',
		     152 => 'Mobile and extrachromosomal element functions',
		     154 => 'Mobile and extrachromosomal element functions',
		     156 => 'Hypothetical proteins',
		     157 => 'Unknown function',
		     158 => 'Protein synthesis',
		     159 => 'Energy metabolism',
		     160 => 'Central intermediary metabolism',
		     161 => 'Amino acid biosynthesis',
		     162 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     163 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     164 => 'Energy metabolism',
		     165 => 'Transcription',
		     166 => 'Transcription',
		     168 => 'Protein synthesis',
		     169 => 'Protein synthesis',
		     170 => 'DNA metabolism',
		     175 => 'Viral functions',
		     176 => 'Fatty acid and phospholipid metabolism',
		     177 => 'Fatty acid and phospholipid metabolism',
		     179 => 'Central intermediary metabolism',
		     182 => 'Transport and binding proteins',
		     183 => 'DNA metabolism',
		     184 => 'Energy metabolism',
		     185 => 'Unclassified',
		     186 => 'Mobile and extrachromosomal element functions',
		     187 => 'Cellular processes',
		     188 => 'Cellular processes',
		     189 => 'Protein fate',
		     191 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     261 => 'Regulatory functions',
		     262 => 'Regulatory functions',
		     263 => 'Regulatory functions',
		     264 => 'Regulatory functions',
		     698 => 'Central intermediary metabolism',
		     699 => 'Signal transduction',
		     700 => 'Signal transduction',
		     701 => 'Cellular processes',
		     702 => 'Cellular processes',
		     703 => 'Unknown function',
		     704 => 'Hypothetical proteins',
		     705 => 'Cellular processes',
		     706 => 'Cellular processes',
		     707 => 'Biosynthesis of cofactors, prosthetic groups, and carriers',
		     708 => 'Mobile and extrachromosomal element functions',
		     710 => 'Signal transduction',
		     719 => 'Cellular processes');
  }
# READ in and store cluster_id to role_id info
# take just first role_id since there can be many
open (CLUSTFILE, "<$basedir/$clusters_file");
while (<CLUSTFILE>)  {
  chomp;
  if (/^cluster/) { next; } # skip the header line
  @b=split(/\t/);
  @roles = split(/,/, $b[4]);
  $role_id = $roles[0]; # take first role_id
  if ($role_id == "") { $role_id = "NONE"; } # for proteins with unassigned role categories
  print "CLUSTER2ROLE_ID: <$b[0]>\t<$b[4]>\t<$role_id>\n" if ($DEBUG);
  $cluster2id_hash{$b[0]} = $role_id;
}
close (CLUSTFILE);
open (GENOATT, ">$basedir/genome.att"); # create genome.att file so we can draw ORFs (option -G in db2circle_rainbow_flatfile.spl)
open (DATAFILE, ">$basedir/data.file"); # create the data file so we can decorate the fGIs on the circle (option -F in db2circle_rainbow_flatfile.spl)
open (ATTFILE, "<$basedir/$att_file");
while (<ATTFILE>)  {
  chomp;
  @c=split(/\t/);
  if (($c[1] =~ /CONTEXT/) || ($c[0] != "1")) { next; } # don't print the contexts or non-chromosome cycles
  $locus = $c[1];
  $end5 = $c[2];
  $end3 = $c[3];
  $com_name = $c[4];
  if ($c[5] =~ /CL/) {
    $cluster = $locus;
    $cluster =~ s/CL_//; # remove CL_ prefix from cluster ids
    $main_role = $role_id_hash{$cluster2id_hash{$cluster}};
    print "ATTFILE_ORFs: $cluster, $end5, $end3, $main_role, $locus, $com_name>\n" if ($DEBUG);
    print GENOATT "ORF\t$end5\t$end3\t$main_role\t$locus\t$com_name\n";
  }
  elsif (($c[5] eq "fGI")) {
    $length = abs($end5-$end3) + 1; # calculate the span of the fGI
    if ($length > "20000") {
      $prefix = "L";
    }
    elsif (($length > "10000") && ($length <= "20000")) {
      $prefix = "M";
    }
    elsif (($length >= "1000") && ($length <= "10000")) {
      $prefix = "S";
    }
    elsif ($length < "1000") {
      next;
    }
    if ($com_name =~ /;/) { # if number insertions has been calculated
      @d=split(/;/, $com_name);
      if ($d[1] >= 10) {
	$inserts = "10+";
      }
      elsif (($d[1] >= 7) && ($d[1] <= 9)) {
	$inserts = "7-9";
      }
      elsif (($d[1] >= 4) && ($d[1] <= 6)) {
	$inserts = "4-6";
      }
      elsif (($d[1] >= 1) && ($d[1] <= 3)) {
	$inserts = "1-3";
      }
      $feat_type = $prefix . "_" . uc($d[0]) . "_" . $inserts;
    }
    else {
      $feat_type = $prefix . "_" . uc($com_name);
    }
    print "ATTFILE_fGIs: <$locus>, <$end5>, <$end3>, <$length bp>, <$feat_type>\n" if ($DEBUG);
    print DATAFILE "$locus\t$end5\t$end3\t$feat_type\n";
  }
}
if ($end3 > $end5) {
  $genlen = $end3;
}
else {
  $genlen = $end5;
}
&print_config_file($basedir,$genlen,$DEBUG);
print DATAFILE "TBORDER\t1\t$genlen\tTOP_BORDER\n";
print DATAFILE "BBORDER\t1\t$genlen\tBOTTOM_BORDER\n";
print DATAFILE "BACK_L\t1\t$genlen\tBKGRD_L\n";
print DATAFILE "BACK_M\t1\t$genlen\tBKGRD_M\n";
print DATAFILE "BACK_S\t1\t$genlen\tBKGRD_S\n";
close (ATTFILE);
close (GENOATT);
close (DATAFILE);
exit;
