#!/usr/bin/perl
# use module

# File: visml_to_svg.pl
#
#
# Description:
# Converts VisANT visml files to Scalable Vector Graphic (SVG) files
#
#
# Requires Modules:  XML::Simple
#                    SVG
#
# Usage:
# ./visml_to_svg.pl input_visml_file.xml output_svg_file.svg

use CGI;
use XML::Simple;
use SVG;
use String::Random;

my $query = new CGI;

my $infile = $query->param("vismlfile");

my $upload_file = $query->upload("vismlfile");

open (UPLOADFILE, ">tmp/$upload_file");

while (<$upload_file>) {

    print UPLOADFILE;

}

close UPLOADFILE;

my $file = "tmp/$upload_file";

print $query->header( );


# create an xml object and basic data structures
my $xml = new XML::Simple;
my %node;
my %namelookup;
my $data = $xml->XMLin("$file");
my %gradient;
my $gradient_index = 0;
my $nodecount = $data->{nodecount};

#figure out the dimensions of our SVG file, graphical properties of different nodes
my $maxx;
my $maxy;
my %methodcolor;
my @methods = keys %{$data->{method}};

for (my $i = 0; $i < $nodecount; $i++) {

    my $name = $data->{Nodes}->{VNodes}->[$i]->{data}->{name};
    my $n = "$name";

    if (($name ne "") & ($data->{Nodes}->{VNodes}->[$i]->{counter} > 0)) {

        $namelookup{$n} = $i;
        $node{$i}{label} = $name;
        $node{$i}{x} = $data->{Nodes}->{VNodes}->[$i]->{x};
        $node{$i}{y} = $data->{Nodes}->{VNodes}->[$i]->{y};
	$node{$i}{size} = $data->{Nodes}->{VNodes}->[$i]->{w} / 1.5;
        $maxx = max($node{$i}{x}, $maxx);
        $maxy = max($node{$i}{y}, $maxy);
  
    }

    if ($data->{Nodes}->{VNodes}->[$i]->{children} ne "") {
	$node{$i}{isMeta} = "true";
    }
    else {
	$node{$i}{isMeta} = "false";
    }

    if ($data->{Nodes}->{VNodes}->[$i]->{group} ne "") {	
	$node{$i}{isChild} = "true";
    }
    else {
	$node{$i}{isChild} = "false";
    }


}

#initialize the SVG file
$maxx += 50;
$maxy += 50;
my $svg= SVG->new(
		  width => $maxx,
		  height => $maxy
		  );

# initialize SVG marker types and colors of different edges in the file
foreach my $m (@methods) {
    my $mc = $data->{method}->{$m}->{color};
    my @c = split(/,/, $mc);
    my $mcolor = "rgb($c[0],$c[1],$c[2])";
    $methodcolor{$m} = $mcolor;
    my $arrowhead = $svg->marker(
                                 'id'     => "$m",
                                 'refX'  => "4",
                                 'refY' => "2.9",
#                                 'viewBox' => "0 0 12 12",
                                 'orient' => "auto",
                                 'markerUnits'  => "strokeWidth",
                                 'markerWidth' => "6",
                                 'markerHeight' => "6",
				 'stroke' => 'none',
				 'fill' => $mcolor
                                 );
    my $arrowline = $arrowhead->path(
				     'd' => "M 0 0 L 6 3 L 0 6 z"
				     );
}

my @color = (100, 130, 100);
my $colorname = "rgb($color[0],$color[1],$color[2])";
make_color_gradient("g$gradient_index", $colorname);
$gradient{$color[0]}{$color[1]}{$color[2]} = $gradient_index;
$gradient_index++;

#initialize node coloring
for (my $i = 0; $i < $nodecount; $i++) {
	
    my $ncc = $data->{Nodes}->{VNodes}->[$i]->{ncc};
    if ($ncc ne "") {

	my @color = split(/\s/, $ncc);
	$colorname = "rgb($color[0],$color[1],$color[2])";

	if ( !(exists $gradient{$color[0]}{$color[1]}{$color[2]}) ) {
	    
	    make_color_gradient("g$gradient_index", $colorname);
	    $node{$i}{color} = "g".$gradient_index;
	    $node{$i}{colorflat} = $colorname;
	    $gradient{$color[0]}{$color[1]}{$color[2]} = $gradient_index;
	    $gradient_index++;
	}
	else {

	    $node{$i}{color} = "g".$gradient{$color[0]}{$color[1]}{$color[2]};
	    $node{$i}{colorflat} = $colorname;

	}
    }
    else {

	if ($node{$i}{isMeta} eq "true") {
	    $node{$i}{colorflat} = "rgb(150,150,150)";
	}
	else {
	    $node{$i}{color} = "g0";
	}
    }
    
}


# calculate edge count data for each node
my @edges = @{$data->{Edges}->{VEdge}};
my $edgecount = @edges;
my %linkcount;

for (my $e = 0; $e < $edgecount; $e++) {

    if (exists $linkcount{$data->{Edges}->{VEdge}->[$e]->{from}} ) {
	$linkcount{$data->{Edges}->{VEdge}->[$e]->{from}} += 1;
    }
    else {
	$linkcount{$data->{Edges}->{VEdge}->[$e]->{from}} = 1;
    }


    if (exists $linkcount{$data->{Edges}->{VEdge}->[$e]->{to}} ) {
        $linkcount{$data->{Edges}->{VEdge}->[$e]->{to}} += 1;
    }
    else {
        $linkcount{$data->{Edges}->{VEdge}->[$e]->{to}} = 1;
    }

}




# draw all edges first
for (my $e = 0; $e < $edgecount; $e++) { 

    
    my $fromnode = $data->{Edges}->{VEdge}->[$e]->{from};
    my $tonode = $data->{Edges}->{VEdge}->[$e]->{to};

    if ( (exists $namelookup{$fromnode}) & (exists $namelookup{$tonode}) ) {
	
	my $from = $namelookup{$fromnode};	
	my $node1x = $node{$from}{x};
	my $node1y = $node{$from}{y};
	my $size1 = $node{$from}{size};
	my $size2;
	
	if ($linkcount{$fromnode} == 1) {
	    
	    if ($fromnode eq $tonode) {
		
		my $method_id = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->{method};		
		my $m = substr($method_id, 0, 5);
		draw_self_edge($node1x, $node1y, $size1, $methodcolor{$m});

	    }
	    
	    else {
		my $to = $namelookup{$tonode};
		my $node2x = $node{$to}{x};
		my $node2y = $node{$to}{y};
		$size2 = $node{$to}{size};
		my $method_id = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->{method};
		my $m = substr($method_id, 0, 5);
                my $totype = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->{toType};
                my $fromtype = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->{fromType};
		draw_line($node1x, $node1y, $node2x, $node2y, $size1, $size2, $methodcolor{$m}, $method_id, 1, 1, $totype, $fromtype);

	    }
	
	}
	else {
	    
	    my $to = $namelookup{$tonode};
	    my $node2x = $node{$to}{x};
	    my $node2y = $node{$to}{y};
	    
	    my $j = 0;
	    
	    while (($j < $linkcount{$fromnode})) {
		
		my $findtonode = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->[$j]->{to};
		$size2 = $node{$to}{size};

		if ($findtonode eq $tonode) {
		    
		    my @methods = split(/,/,$data->{Nodes}->{VNodes}->[$from]->{data}->{link}->[$j]->{method});

		    if ($fromnode eq $tonode) {

			my $method_id = $methods[0];			
			my $m = substr($method_id, 0, 5);
			draw_self_edge($node1x, $node1y, $size1, $methodcolor{$m});


		    }

		    else {

			my $nm = @methods;
			my $k = 1;
			
			foreach my $method (@methods) {
			    
			    my $m = substr($method, 0, 5);
			    my $totype = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->[$j]->{toType};
			    my $fromtype = $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->[$j]->{fromType};
			    my $weight = max(0, $data->{Nodes}->{VNodes}->[$from]->{data}->{link}->[$j]->{weight});
			    draw_line($node1x, $node1y, $node2x, $node2y, $size1, $size2, $methodcolor{$m}, $m, $k, $nm, $totype,$fromtype, $weight);
			    $k++;
			    
			}
		    }
		    $j++;
		}
		else {
		    
		    $j++;
		}
	    }
	}
    }	
}

#draw child nodes
for (my $i = 0; $i < $nodecount; $i++) {

    if (exists $node{$i}{label} &  $node{$i}{isChild} eq "true" & $node{$i}{isMeta} eq "false") {

        my $tag = $svg->circle(cx=>$node{$i}{x},
                               cy=>$node{$i}{y},
                               r=>$node{$i}{size},
                               fill=> "url(#$node{$i}{color})",
                                'fill-opacity'=>'0.5',
                               );
    }
}


#draw meta (parent) nodes
for (my $i = 0; $i < $nodecount; $i++) {

    if (exists $node{$i}{label} & $node{$i}{isMeta} eq "true") {

	if ($data->{Nodes}->{VNodes}->[$i]->{childVisible} eq "false") {
	    
	    my $tag = $svg->circle(cx=>$node{$i}{x},
				   cy=>$node{$i}{y},
				   r=>$node{$i}{size},
				   'fill' => $node{$i}{colorflat},
				   'fill-opacity' => '0.5',
				   'stroke' => $node{$i}{colorflat},
				   'stroke-width' => 4,
				   'stroke-opacity' => '1'
				   );
	}
	else {

	    my $w = $data->{Nodes}->{VNodes}->[$i]->{w};
            my $h = $data->{Nodes}->{VNodes}->[$i]->{h};


            my $tag = $svg->rectangle(
				      x=>$node{$i}{x} - ($w/2),
				      y=>$node{$i}{y} - ($h/2),
				      width=>$w,
				      height=>$h,
				      rx=>10,
				      ry=>10,
				      'fill' => $node{$i}{colorflat},
				      'fill-opacity' => '0.3',
				      'stroke' => $node{$i}{colorflat},
				      'stroke-width' => 4,
				      'stroke-opacity' => '0.5'
				      );

	}


    }
}


#draw remaining free (parentless) nodes
for (my $i = 0; $i < $nodecount; $i++) {

    if (exists $node{$i}{label} ) {

	if ( $node{$i}{isMeta} ne "true" & $node{$i}{isChild} ne "true") {

	    my $tag = $svg->circle(cx=>$node{$i}{x},
				   cy=>$node{$i}{y},
				   r=>$node{$i}{size},
				   fill=> "url(#$node{$i}{color})",
				   'fill-opacity'=>'0.5',
				   );
	}
    }
}
	   

for (my $i = 0; $i < $nodecount; $i++) {

    if ((exists $node{$i}{label}) & ( $data->{Nodes}->{VNodes}->[$i]->{counter} != 0)) {

        if ( $data->{Nodes}->{VNodes}->[$i]->{labelOn} eq "true" ) {
            my $alias =  $data->{Nodes}->{VNodes}->[$i]->{vlabel};
	    my $nodesize = $node{$i}{size};
            my $text = $svg->text(
                                  'x'     => $node{$i}{x} - ($nodesize),
                                  'y'     => $node{$i}{y} + ($nodesize/3),
                                  'font-size' => $data->{Nodes}->{VNodes}->[$i]->{labelSize},
				  'font'  => "Vera"
                                  )->cdata($alias);
        }

    }
}




# use svg namespace and generate a document with its own DTD


my $string = new String::Random;
my $outfile = $string->randpattern("cccnnn");
open(OUT, ">>tmp/$outfile.svg");
print OUT $svg->xmlify;

`rsvg -x .2 -y .2 -f png tmp/$outfile.svg tmp/$outfile.png`;

print <<END_HTML;
<HTML>
<HEAD>
<TITLE>
Converter
</TITLE>
</HEAD>
<BODY>
<P>
You converted <b>$infile</b> to 
<A HREF="tmp/$outfile.svg" TARGET="_new">$outfile.svg</A>, a scaled-down version is shown below.
Right-click on link to save SVG.
</P>
<P>
<IMG SRC="tmp/$outfile.png" />
</P>
</BODY>
</HTML>
END_HTML

###############
# SUBROUTINES##
###############


sub max {
    if ($_[0] > $_[1]) {
	return $_[0];
    }
    else {
	return $_[1];
    }
}


sub min {
    if ($_[0] < $_[1]) {
        return $_[0];
    }
    else {
        return $_[1];
    }
}

# draw arbitrary line with arbitrary direction and marker type
sub draw_line {
    
    my $node1x = $_[0] + 0.00001;
    my $node1y = $_[1] + 0.00001;
    my $node2x = $_[2];
    my $node2y = $_[3];
    my $size1 = $_[4];
    my $size2 = $_[5];
    my $color = $_[6];
    my $marker = $_[7];
    my $k = $_[8];
    my $nm = $_[9];
    my $totype = $_[10];
    my $fromtype = $_[11];
    my $weight = ($_[12]*9) + 1;

    my $m = ($node2y - $node1y + .00000001) /  ($node2x - $node1x + 0.0000001) ;
    
    my ($x1, $x2, $y1, $y2);

    if ($node1x < $node2x & $node1y < $node2y) {
	$x1 = $node1x + ($size1 / sqrt(1 + ($m**2)));
        $x2 = $node2x - ($size2 / sqrt(1 + ($m**2)));
        $y1 = $node1y + ($m *($size1 / sqrt(1 + ($m**2))));
        $y2 = $node2y - ($m *($size2 / sqrt(1 + ($m**2))));
    }
    elsif ($node1x < $node2x & $node1y > $node2y )  {
        $x1 = $node1x + ($size1 / sqrt(1 + ($m**2)));
        $x2 = $node2x - ($size2 / sqrt(1 + ($m**2)));
        $y1 = $node1y + ($m *($size1 / sqrt(1 + ($m**2))));
        $y2 = $node2y - ($m *($size2 / sqrt(1 + ($m**2))));

    }
    elsif ($node1x > $node2x & $node1y < $node2y) {
        $x1 = $node1x - ($size1 / sqrt(1 + ($m**2)));
        $x2 = $node2x + ($size2 / sqrt(1 + ($m**2)));
        $y1 = $node1y - ($m *($size1 / sqrt(1 + ($m**2))));
        $y2 = $node2y + ($m *($size2 / sqrt(1 + ($m**2))));

    }
    elsif ($node1x > $node2x & $node1y > $node2y) {
        $x1 = $node1x - ($size1 / sqrt(1 + ($m**2)));
        $x2 = $node2x + ($size2 / sqrt(1 + ($m**2)));
        $y1 = $node1y - ($m *($size1 / sqrt(1 + ($m**2))));
        $y2 = $node2y + ($m *($size2 / sqrt(1 + ($m**2))));
    }


    my $starty = $y1 + (($y2 - $y1) * (($k-1)/$nm));
    my $stopy = $y1 + (($y2 - $y1) * ($k/$nm));

    my $startx = $x1 + (($x2 - $x1) * (($k-1)/$nm));
    my $stopx = $x1 + (($x2 - $x1) * (($k)/$nm));

    my $d = sqrt( ($stopy - $starty)**2  + ($stopx - $startx)**2 );

    if ($totype == 1) {

	if ($k == $nm) {	    

	    my $tag = $svg->line(
				 x1=>round($startx), 
				 y1=>round($starty),
				 x2=>round($stopx), 
				 y2=>round($stopy),
				 style=>{
				     'stroke'=> $color,
				     'stroke-width'=>$weight
				     }
				 );
	    my $tag = $svg->line(
				 x1=>round($startx),
				 y1=>round($starty),
				 x2=>round($stopx),
				 y2=>round($stopy),
				 style=>{
                                     'marker-end' => "url(#$marker)",
                                     'stroke'=> $color,
                                     'stroke-width'=>$weight
                                     }
                                 );
	}
	else {
	    my $tag = $svg->line(
				 x1=>round($startx),
				 y1=>round($starty),
				 x2=>round($stopx),
				 y2=>round($stopy),
				 style=>{
				     'stroke'=> $color,
				     'stroke-width'=>$weight
				     }
				 );
	}
    }

    elsif ($fromtype == 1) {


	if ($k == 1) {	
    
	    my $tag = $svg->line(
				 x2=>round($startx), 
				 y2=>round($starty),
				 x1=>round($stopx), 
				 y1=>round($stopy),
				 style=>{
				     'stroke'=> $color,
				     'stroke-width'=>$weight
				     }
				 );
	    my $tag = $svg->line(
				 x2=>round($startx),
				 y2=>round($starty),
				 x1=>round($stopx),
				 y1=>round($stopy),
				 style=>{
                                     'marker-end' => "url(#$marker)",
                                     'stroke'=> $color,
                                     'stroke-width'=>$weight
                                     }
                                 );
	}
	else {
	    my $tag = $svg->line(
				 x1=>round($startx),
				 y1=>round($starty),
				 x2=>round($stopx),
				 y2=>round($stopy),
				 style=>{
				     'stroke'=> $color,
				     'stroke-width'=>$weight
				     }
				 );
	}
    }
   else {

       my $tag = $svg->line(
			    x2=>round($startx),
			    y2=>round($starty),
			    x1=>round($stopx),
			    y1=>round($stopy),
			    style=>{
				'stroke'=> $color,
				'stroke-width'=>$weight
				}
			    );
   }

}

# defines a color gradient
sub make_color_gradient {

    my $gradient = $svg->gradient(
                                  '-type' => "radial",
                                  'id'    => "$_[0]",
                                  'cx'    => "35%",
                                  'cy'    => "35%",
                                  'r'     => "50%"
                                  );
    my $pattern = $gradient->stop(
                                  'offset'     => "5%",
                                  'stop-color' => 'rgb(255,255,255)'
                                  );
    my $pattern = $gradient->stop(
                                  'offset'     => "95%",
                                  'stop-color' => $_[1]
                                  );
}

# draws a looping edge from a node to itself
sub draw_self_edge  {

    my $x = $_[0];
    my $y = $_[1];

    my $size = $_[2];
    my $color = $_[3];
    
    my $startx = round($x - ($size/1.4));
    my $starty = round($y - ($size/1.4));

    my $stopx = round($x + ($size/1.4));
    my $stopy = round($y - ($size/1.4));

    my $rx = round($size);
    my $ry = round($size**.9);
    
    my $tag = $svg->path(
		      'd' => "M$startx,$starty A$rx $ry 180 1 1 $stopx $stopy",
		      style => {
			  'stroke'=> $color,
			  'fill' => 'none',
			  'stroke-width'   => '1',
			  }
		      );
    
}


sub round {
    my $r = sprintf("%.1f", $_[0]);
    return $r;
}