###############################################################################
#                                                                             #
#       Copyright (C) 2016-2017 J. Craig Venter Institute (JCVI).             #
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

# Note, depths have been assigned as follows

# bckgrnd = 10
#     1.00 = 9 = base gene or base island
#     0.80 = 8
#     0.75 = 7 = threequarter
#     0.60 = 6
#     0.50 = 5 = half
#     0.40 = 4
#     0.25 = 3 = quarter
#     0.20 = 2
# features = 1
#  borders = 0

sub draw_elements { # modified by Derrick E. Fouts 01/30/2007 for drawing different heights based on %id
    my($end5, $end3, $l, $c, $id, $type, $comment, $seq_len) = @_;
    my($h1,$h2,$i,$multiplier,$depth,$t);

    $h1 = convert_layer($l);

    if ($id eq "BCKGRND")  {
	$multiplier = "1";
        $depth = "10";
    }
    elsif (($id == "100") | ($id eq "ISLAND") | ($id eq "GENE"))  {
	$multiplier = "1";
        $depth = "9";
    }
    elsif ($id == "98")  {
	$multiplier = "0.8";
	$depth = "8";
    }
    elsif ($id eq "THREEQUARTER")  {
	$multiplier = "0.75";
	$depth = "7";
    }
    elsif ($id == "95")  {
	$multiplier = "0.6";
	$depth = "6";
    }
    elsif ($id eq "HALF")  {
	$multiplier = "0.5";
	$depth = "5";
    }
    elsif ($id == "90")  {
	$multiplier = "0.4";
	$depth = "4";
    }
    elsif ($id eq "QUARTER")  {
	$multiplier = "0.25";
	$depth = "3";
    }
    elsif ($id == "85")  {
	$multiplier = "0.2";
	$depth = "2";
    }
    elsif (($id eq "WIDE") && ($type eq "GENE")) {
	$multiplier = "1";
        $depth = "9";
	$width = "15";
    }
    if ($id eq "POINTER")  {
	$multiplier = "0.1";
        $depth = "1";
    }
     
    $h2 = $h1 + ($gene_height*$multiplier);

    ($end5, $end3) = ($end3, $end5) if ($end5 > $end3);

    if ($type eq "ISLAND")  {
	if (($end3-$end5) < 100)  {  # to correct for case where the region to print is smaller than default size of tick (100 bp)
	    $end3 = $end5+100;
	}

	if ($seq_len < 500000)  {
	    $t = 10;
	}
	else  {
	    $t = 100;
        }

	for($i=$end5;$i<$end3;$i+=$t) {
	    if ($GIF_TOGGLE) {
		gif_spoke($i,$h1,$h2,$c);
	    }
	    else {
		xfig_spoke($depth,$i,$h1,$h2,$c);
	    }
	}
    }
    elsif ($type eq "GENE")  {
	if ($GIF_TOGGLE) {
	    gif_spoke($end5, $h1, $h2, $c, $comment);
	}
	else {
	    xfig_spoke($depth, $end5, $h1, $h2, $c, $width);
	}
    }
}


##################################################################

sub GetAvg {
    my(%BarGraph) = @_; # changed 06/10/2004 because Perl v5.8.1-RC3 built for darwin-thread-multi-2level did not like the local(\%BarGraph) = @_;
    my($avg) = 0;
    my($total_values) = 0;
    my($num_values) = 0;
    foreach $key(keys %BarGraph) {
	$total_values += $BarGraph{$key}->{'value'};
	$num_values++;
    }
    $avg = $total_values/$num_values;
    return($avg);
}

sub avg_value_bar_graph { 
    my($end5, $v, $l, $c,$avg_bar) = @_;
    my($h1,$h2,$i);

    $h1 = convert_layer_avg($l);
    
    if($v > $avg_bar) {
	$h2 = $h1 + (($gene_height * abs($v-$avg_bar)) ) ;
	
    }else{
	$h2 = $h1 - (($gene_height * abs($v-$avg_bar)) ) ;
	
    }

    ($end5, $end3) = ($end3, $end5) if ($end5 > $end3);

    if ($GIF_TOGGLE) {
	gif_spoke($end5,$h1,$h2,$c);
    }
    else {
	xfig_spoke(0, $end5,$h1,$h2,$c);
    }

}

sub convert_layer_avg {
    my($l) = @_;

    if ($l > $NUMBER_OF_LAYERS) {
	print STDERR "MAXIMUM LAYERS EXCEEDED: $l\n";
	print STDERR "CONVERTING TO ONE\n";
	$l = 1;
    }

    return($max_radius - ($layer_height * $l) + ($layer_height/2));
}

################################################################

sub draw_bar_graph { 
    my($end5, $v, $l, $c) = @_;
    my($h1,$h2,$i);

    $h1 = convert_layer($l);

    $h2 = $h1 + ($gene_height * $v);

    ($end5, $end3) = ($end3, $end5) if ($end5 > $end3);

    if ($GIF_TOGGLE) {
	gif_spoke($end5,$h1,$h2,$c);
    }
    else {	
	xfig_spoke(0, $end5,$h1,$h2,$c);
    }

}

sub gif_tick { 
    my($x, $l) = @_;
    my($h1, $h2);

    $h1 = convert_layer($l);

    $h2 = $h1 + $tick_height;

    gif_spoke($x,$h1,$h2,$BLACK);

}

sub xfig_tick { 
    my($x, $l) = @_;
    my($h1, $h2);

    $h1 = convert_layer($l);

    $h2 = $h1 + $tick_height;

    xfig_spoke(0, $x,$h1,$h2,$BLACK);

}

sub xfig_span { 
    my($end5, $end3, $l) = @_;
    my($h1, $h2, $h3);

    $h1 = convert_layer($l);

    $h3 = $h1 + $span_height;
    $h2 = $h1 + ($span_height / 2);

    ($end5, $end3) = ($end3, $end5) if ($end5 < $end3);

    xfig_spoke(0, $end5,$h1,$h3,$BLACK);
    xfig_spoke(0, $end3,$h1,$h3,$BLACK);
    xfig_draw_arc($end5, $end3, $h2, $BLACK);

}

sub xfig_draw_arc {
    my($d1,$d3,$h,$c) = @_;
    
    ($d1, $d3) = ($d3, $d1) if ($d1 < $d3);
    $d2 = (($d3-$d1) / 2) + $d1;

    printf("5 1 0 1 %d 7 0 0 -1 0.000 0 0 0 0 ",$c);
    printf("%lf %lf %d %d %d %d %d %d\n",
	   get_center_x(), get_center_y(),
	   screen_x_trans($d1,$h), screen_y_trans($d1,$h),
	   screen_x_trans($d2,$h), screen_y_trans($d2,$h),
	   screen_x_trans($d3,$h), screen_y_trans($d3,$h));


}

sub place_gif_circle {
    my($l) = @_;
    my($radius);

    $radius = $max_radius - ($layer_height * $l);

    $image->arc(get_center_x(), 
		get_center_y(), 
		2*$radius, 
		2*$radius, 
		0, 360, $GIF_BLACK);

}

sub place_gif_ticks {
    my($increment) = @_;
    my($h,$i,$j,$d,$justification);

    # this puts a tick mark at intervals around the circle
    for($i=0;$i<$seq_len;$i+=$increment) {
        gif_tick($i, 0);

        $d = 360 * ($i / $seq_len);

        $justification = 0;
        if ($d > 180) {
            $justification = 2;
        }

        if ($d > 120 && $d < 240) {
            $h = convert_layer(0) + $tick_height + $text_height;
        }
        else {
            $h = convert_layer(0) + $tick_height + ($text_height / 2);
        }

        $j = $i;
        $j = 1 if ($j == 0);


        # THIS IS WHERE LABELS ARE PLACED:
        ### due to the different locations around the pie, modification is
        ### required 
        $degrees = (360 * ($j / $seq_len));
        $lngth = (length(commas($j)) * 5)+5;
        if (($degrees > 0) && ($degrees <= 90))
        {
            $mod_x = 0;
            $mod_y = -5;  ## height of one of the character
        }
        if (($degrees > 90) && ($degrees <= 180))
        {
            $mod_x = 0;
            $mod_y = 0;
        }
        if (($degrees > 180) && ($degrees <= 270))
        {
            $mod_x = -$lngth;
            $mod_y = 0;
        }
        if (($degrees > 270) && ($degrees <= 360))
        {
            $mod_x = -$lngth;
            $mod_y = -5; ### again the height
        }
        
        $image->string(gdSmallFont,
         screen_x_trans($j,$h)+$mod_x,
         screen_y_trans($j,$h)+$mod_y,
         commas($j),
         $BLACK);
    }
}


sub place_xfig_circle {
    my($l,$t,$gh) = @_;
    my($radius);

    $radius = $max_radius - ($layer_height * $l)+$gh;

    printf("1 3 0 $t -1 7 0 0 -1 0.000 1 0.0000 ");
    printf("%d %d %d %d %d %d %d %d\n",
	   get_center_x(), 
	   get_center_y(),
	   $radius,
	   $radius,
	   get_center_x(), 
	   get_center_y(),
	   get_center_x(), 
	   get_center_y() + $radius);

}

sub place_xfig_ticks {
    my($increment) = @_;
    my($h,$i,$j,$d,$justification);

    # this puts a tick mark at intervals around the circle
    for($i=0;$i<$seq_len;$i+=$increment) {
	xfig_tick($i, 0);

	$d = 360 * ($i / $seq_len);

	$justification = 1;
	if ($d >= 5 && $d < 175) {
	    $justification = 0;
	} elsif ($d >= 185 && $d < 355 ) {
	    $justification = 2;
	}

	if ($d > 120 && $d < 240) {
	    $h = convert_layer(0) + $tick_height + $text_height;
	}
	else {
	    $h = convert_layer(0) + $tick_height + ($text_height / 2);
	}

	$j = $i;
	$j = 1 if ($j == 0);

	# label that tick mark.
	printf("4 %d -1 0 0 16 60 0.0000 4 765 3225 0", $justification);
	printf("%d %d %s\\001\n", 
	       screen_x_trans($j,$h),
	       screen_y_trans($j,$h),
	       commas($j));

    }
}

sub gif_spoke {
    my($d,$l1,$l2,$c,$comment) = @_;
    
    $image->line(screen_x_trans($d,$l1), 
		 screen_y_trans($d,$l1),
		 screen_x_trans($d,$l2), 
		 screen_y_trans($d,$l2),
		 $c);

}

sub xfig_spoke {
    my($depth,$d,$l1,$l2,$c,$width) = @_;
    # width is optional

    if (length($width) == 0) {
	$width = $spoke_size;
    }
    
    printf("2 1 0 %d %d 7 $depth 0 -1 0.000 0 0 -1 0 0 2 0\n",
	   $width,
	   $c);
    printf("%d %d %d %d\n",
	   screen_x_trans($d,$l1), screen_y_trans($d,$l1),
	   screen_x_trans($d,$l2), screen_y_trans($d,$l2));
}

sub calc_ave {
    my($layer, %BarGraph) = @_;
    my($avg) = 0;
    my($total_values) = 0;
    my($num_values) = 0;

    foreach $key(keys %BarGraph) {
	if ($BarGraph{$key}->{'layer'} == $layer) {

	    $min = $BarGraph{$key}->{'value'};
	    $max = $BarGraph{$key}->{'value'};

	    last;
	}
    }

    foreach $key(keys %BarGraph) {
	if ($BarGraph{$key}->{'layer'} == $layer) {

	    $min = min($min,$BarGraph{$key}->{'value'});
	    $max = max($max,$BarGraph{$key}->{'value'});

	    $total_values += $BarGraph{$key}->{'value'};
	    $num_values++;
	}
    }

    $avg = $total_values/$num_values;

    return($min, $max, $avg);
}


########### Don't think this works very well ############
#sub avg_value_bar_graph { 
#    local($end5, $v, $l, $c, $min, $max, $avg) = @_;
#    local($h1,$h2,$q);
#
#    $h1 = convert_layer($l) + ($layer_height / 2);
#
#    $q = ($layer_height / 2) / ($max - $min);
#
#    $FUDGE_FACTOR = 1.5;
#
#    $h2 = $h1 + (($v - $avg) * $q * $FUDGE_FACTOR);
#
#    if ($GIF_TOGGLE) {
#	gif_spoke($end5,$h1,$h2,$c);
#    }
#    else {
#	xfig_spoke($end5,$h1,$h2,$c);
#    }
#
#}


sub convert_layer {
    my($l) = @_;

    if ($l > $NUMBER_OF_LAYERS) {
	print STDERR "MAXIMUM LAYERS EXCEEDED: $l\n";
	print STDERR "CONVERTING TO ONE\n";
	$l = 1;
    }

    return($max_radius - ($layer_height * $l));
}

sub screen_x_trans {
    my($p, $r) = @_;
    my($d,$x);

    $d = 360 * ($p / $seq_len);
    return(($r * sin(deg2rad($d))) + $max_radius + $border);
}

sub screen_y_trans {
    my($p, $l) = @_;
    my($d);

    $d = 360 * ($p / $seq_len);

    return(($l * cos(deg2rad($d)) * -1) + $max_radius + $border);
}

sub rad2deg {
    my ($r) = @_;

    return(((180 / $PI) * ($r)) % 360);
}

sub deg2rad {
    my ($d) = @_;

    return(($d/180) * $PI);
}

sub init_xfig_drawing {

    print "\#FIG 3.2\n";
    print "Landscape\n";
    print "Center\n";
    print "Inches\n";
    print "Letter  \n";
    print "100.00\n";
    print "Single\n";
    print "-2\n";
    print "1200 2\n";

    for($i=$COLOR_START;length($rgb{$i}) != 0;$i++) {
	printf("0 %d \#%s\n",$i,$rgb{$i});
    }

}

sub get_center_x {
    return($max_radius + $border);
}

sub get_center_y {
    return($max_radius + $border);
}

sub commas { 
    my($_) = @_;
    1 while s/(.*\d)(\d\d\d)/$1,$2/;
    $_;
}

sub init_image  {
    my($x,$y) = @_;
    my($image) = new GD::Image($x,$y);
        
    for($i=$COLOR_START;length($rgb{$i}) != 0;$i++) {

	$gif_role_color{$i} = $image->colorAllocate(hex(substr($rgb{$i},0,2)),
						hex(substr($rgb{$i},2,2)),
						hex(substr($rgb{$i},4,2)));

    }

    $BACK_GROUND = $image->colorAllocate(255,255,255);
    $GIF_BLACK = $image->colorAllocate(0,0,0);
    $GIF_WHITE = $image->colorAllocate(255,255,255);

    $image->fill(2,2,$BACK_GROUND);
    return($image);
}

sub write_image_to_disk {
    my($f,$image) = @_;
        
    if(open(IMAGE,"+> $f")) {
        print IMAGE  $image->gif;
        close(IMAGE);
	# might be helpful to do this so the web browser
	#  can read it.
        chmod(0664,$f);
    }
}

sub write_image_to_stdout {
    my($image) = @_;
        
    print $image->gif;

}

sub max {
    my($x,$y) = @_;
    return ($x >= $y) ? $x :$y;
}
        
sub min {
    my($x,$y) = @_;
    return ($x < $y) ? $x :$y;
}

1;
