package HTML::TreeStructured;
use 5.006;
use strict;
use warnings;

use Carp qw(croak);
use HTML::Template 2.6;

=head1 NAME

HTML::TreeStructured - Perl extension for generating tree structured HTML

=head1 SYNOPSIS

	use HTML::TreeStructured;

	### Describe tree via collection of Node and its properties

	### Method 1: Via ArrayRef
	###
	### Node can be a string or '/' concatenated strings to show ancestry
	### Properties are name/value pairs

	my $tree1 = [
		['/aaa', 	color => 'green'],
		['/aaa/bbb'	mouseover => 'This is addl info'],
		['/aaa/ccc',	color => 'red', active => 0]
	];

	### Method 2: Via Hashref

	my $tree2 = {
		aaa => {
			color => 'green',
			bbb   => {
				mouseover => 'This is addl info',
			},
			ccc   => {
				color	=> 'red',
				active	=> 0,
			},
	};

	Interpreted Node Properties:

	color		= Color of the node name
	mouseover	= Mouse Over text for the node (Info image is displayed next to node)
	active		= 0 would cause strike thru on node
	highlight	= color code used for marker highlight of node
	url		= URL to hyperlink the node
	tooltip		= popup when mouse is over the link (together with url) (See HTML::Tooltip::Javascript)
	closed		= 1 if node be closed on default display (default, all nodes are open)
	comment		= Text to display next to node in bold
	weight		= A numeric value on node which will be used for sorting node position in at sibling level
			  (Default, nodes are sorted in ascending order per dictionary order)


	### Now get HTML equivalent for the tree
	### The associated JavaScript for nodes close/open and ExpandAll/CollapseAll is generated alongside

	$tree_html = HTML::TreeStructured->new(
		name         => 'tree_name',
		image_path   => '/images/',
		data         => $tree1,
		title        => "My Tree",
		title_width  => 300,
		level        => {},     ### If scalar, close BEYOND this depth. Depth start at 0.
					### If Hash, close for depths specified in keys
	)->output;

	### The same module can be used to generate FAQ - see "examples/faq.cgi"

=cut

our $VERSION = "1.01";

sub new 
{
    my $pkg = shift;
    
    # setup defaults and get parameters
    my $self = bless({ 

			title_width	=> 300,
			child_indent	=> 20,
			image_path	=> ".",
			level	     => {},	### Default open up the full tree - This would contain 
		       				### Scalar => depth BEYOND which nodes are closed for tree display
		       				### Hash   => depths matching keys are closed for tree display
						### NB: Depth starts at ZERO
						### E.g. Value of
						### level=2 means Tree nodes are closed at depth 3 and more (all others are open)
						### level={2=>1} means Tree node at depth 2 are closed (all others are open)
                       @_,
                     }, $pkg);
    
    # fix up image_path to always end in a /
    $self->{image_path} .= "/" unless $self->{image_path} =~ m!/$!;

    # check required params
    foreach my $req (qw(name data title)) {
        croak("Missing required parameter '$req'") unless exists $self->{$req};
    }

    if (ref($self->{data}) eq 'ARRAY') {
    	$self->{data} = process_arrayref($self->{data});
    } else {
	my $res;
    	my $data = $self->{data};
	my @kkk = grep { ref($data->{$_}) eq 'HASH' } keys %$data;
	if (1 == @kkk) {
		$res = process_hashref($kkk[0], $data->{$kkk[0]});
	} else {
		$res = process_hashref('ROOT', $data);
	}
	$self->{data} = $res;
    }

    #use Data::Dumper;
    #print '<pre>', Dumper($self->{data}), '</pre>';

    return $self;
}

sub output 
{
    my $self = shift;
    our $TEMPLATE_SRC;
    my $template = HTML::Template->new(scalarref          => \$TEMPLATE_SRC,
                                       die_on_bad_params => 0,
                                       global_vars       => 1,
                                      );

    # build node loop
    my @loop;
    $self->_output_node(node   => $self->{data},
                        loop   => \@loop,
			depth  => 1,
			level  => $self->{level},
                       );
    my @parents;			### Collect all nodes with children - for use in ExpandAll/CollapseAll
    map { push(@parents, {id => $_->{id}}) if ($_->{has_children}) } @loop;
    # setup template parameters
    $template->param(loop => \@loop);
    $template->param(parents => \@parents);
    $template->param(map { ($_, $self->{$_}) } qw(name title title_width child_indent image_path));
    # get output for the widget
    my $output = $template->output;

    return $output;
}

# recursively add nodes to the output loop
sub _output_node 
{
    my ($self, %arg) = @_;
    my $node = $arg{node};
    my $depth = $arg{depth};
    my $level = $arg{level};

    #use Data::Dumper;
    #print "<pre>", Dumper($node), "</pre>";

    my $id = next_id();
    push @{$arg{loop}}, { label       => $node->{label},		### Label to appear in tree
                          value       => $node->{value},		### Hidden Value (whats the use?, but good to have)
                          id          => $id,				### Unique Id (no immediate use, but good to have)
                          open        => display_closed_node($depth, $level, $node->{closed}) ? 0 : 1,	
			  						### During Display, whether to close/open
			  url	      => $node->{url},			### Url to link 
			  mouseover   => $node->{mouseover},		### mouseover message
			  tooltip     => $node->{tooltip},		### Tooltip popup box (See HTML::Tooltip::Javascript)
			  active      => ((defined($node->{active}) and $node->{active} == 0) ? 0 : 1),	
			  						### Is this node active? If not, strike thru during display
			  color       => $node->{color} || 'black',	### Color of the label
			  highlight   => $node->{highlight},		### Color to use for marker highlight
			  comment     => $node->{comment},		### Comment in bold within parentheses next to label
                        };
    
    if ($node->{children} and @{$node->{children}}) {
        $arg{loop}[-1]{has_children} = 1;
        for my $child (@{$node->{children}}) {
            $self->_output_node(node   => $child,
                                loop   => $arg{loop},
				depth  => $depth + 1,
				level  => $level,
                               );
        }
        push @{$arg{loop}}, { end_block => 1 };
    }
    
}

sub display_closed_node
{
	my $ddd = shift; # Depth Info
	my $lll = shift; # Level Info 
			 # scalar ==> Close all nodes BEYOND this depth
			 # hash   ==> Close all nodes for specified keys
	my $closed = shift;

	if (defined($closed)) {
		return $closed;
	}

	if (ref($lll) eq 'HASH') {
		return ($lll->{$ddd} ? 1 : 0);
	} else {
		return ($ddd > $lll ? 1 : 0);
	}
}

{ 
    my $id = 1;
    sub next_id { $id++ }
}

our $TEMPLATE_SRC = <<END;
<style type="text/css">
<!--

  /* title bar style.  The width here will define a minimum width for
     the widget. */
  .hpts-title {
     padding:          2px;
     margin-bottom:    4px;     
     font-size:        large;
     color:            #ffffff;
     background-color: #666666;
     width:            <tmpl_var title_width>px;
  }

  /* style of a block of child nodes - indents them under their parent
     and starts them hidden */
  .hpts-block {
     margin-left:      <tmpl_var child_indent>px;
     display:          none;
  }

  /* style for selected labels */
  .hpts-label-selected {
     background:       #98ccfe;
  }

  /* style for labels after being unselected */
  .hpts-label-unselected {
     background:       #ffffff;
  }

-->
</style>

<script language="javascript">

  /* expand or collapse a sub-tree */
  function <tmpl_var name>_toggle_expand(id) {
     var obj = document.getElementById("<tmpl_var name>-desc-" + id);
     var plus = document.getElementById("<tmpl_var name>-plus-" + id);
     var node = document.getElementById("<tmpl_var name>-node-" + id);
     if (obj.style.display != 'block') {
        obj.style.display = 'block';
        plus.src = "<tmpl_var image_path>minus.png";
        node.src = "<tmpl_var image_path>open_node.png";
     } else {
        obj.style.display = 'none';
        plus.src = "<tmpl_var image_path>plus.png";
        node.src = "<tmpl_var image_path>closed_node.png";
     }
  }

  /* expand or collapse a sub-tree */
  function <tmpl_var name>_Expand(id) {
     var obj = document.getElementById("<tmpl_var name>-desc-" + id);
     var plus = document.getElementById("<tmpl_var name>-plus-" + id);
     var node = document.getElementById("<tmpl_var name>-node-" + id);
        obj.style.display = 'block';
        plus.src = "<tmpl_var image_path>minus.png";
        node.src = "<tmpl_var image_path>open_node.png";
  }

  /* expand or collapse a sub-tree */
  function <tmpl_var name>_Collapse(id) {
     var obj = document.getElementById("<tmpl_var name>-desc-" + id);
     var plus = document.getElementById("<tmpl_var name>-plus-" + id);
     var node = document.getElementById("<tmpl_var name>-node-" + id);

        obj.style.display = 'none';
        plus.src = "<tmpl_var image_path>plus.png";
        node.src = "<tmpl_var image_path>closed_node.png";
  }

  function <tmpl_var name>_CollapseAll() {
	<tmpl_loop parents>
		<tmpl_var name>_Collapse(<tmpl_var id>);
	</tmpl_loop>
  }

  function <tmpl_var name>_ExpandAll() {
	<tmpl_loop parents>
		<tmpl_var name>_Expand(<tmpl_var id>);
	</tmpl_loop>
  }

</script>


<div id="<tmpl_var name>-outer" style="display: block">
  <div class="hpts-title" id="<tmpl_var name>-title>"><tmpl_var title></div>
  <a href="JavaScript:<tmpl_var name>_ExpandAll()">ExpandAll</a> | <a href="JavaScript:<tmpl_var name>_CollapseAll()">CollapseAll</a>
  <div>
  <tmpl_loop loop>
    <tmpl_unless end_block>
       <div>
          <tmpl_if has_children>
              <img id="<tmpl_var name>-plus-<tmpl_var id>" width=16 height=16 src="<tmpl_var image_path><tmpl_if open>minus<tmpl_else>plus</tmpl_if>.png" alt="image" onclick="<tmpl_var name>_toggle_expand(<tmpl_var id>)"><span id="<tmpl_var name>-line-<tmpl_var id>" ondblclick="<tmpl_var name>_toggle_expand(<tmpl_var id>)">
	  <tmpl_else>
              <img width=16 height=16 src="<tmpl_var image_path>L.png" alt="image"><span id="<tmpl_var name>-line-<tmpl_var id>">
          </tmpl_if>
          <img id="<tmpl_var name>-node-<tmpl_var id>" width=16 height=16 src="<tmpl_var image_path><tmpl_if has_children><tmpl_if open>open<tmpl_else>closed</tmpl_if><tmpl_else>file</tmpl_if>_node.png" alt=<tmpl_if has_children>"Double click to toggle (Open/Close)"<tmpl_else>"image"</tmpl_if>>

	  <tmpl_if form_input_file>
	  <tmpl_unless has_children>
	  	<input type=file size=50 name="<tmpl_var field_name>">
	  	<input type=hidden name="<tmpl_var field_name>_">
	  </tmpl_unless>
	  </tmpl_if>

	  <tmpl_if form_input_text>
	  <tmpl_unless has_children>
	  	<input type=text size=50 name="<tmpl_var field_name>" onblur="document.<tmpl_var name>['<tmpl_var field_name>_'].value=document.<tmpl_var name>['<tmpl_var field_name>'].value" onchange="document.<tmpl_var name>['<tmpl_var field_name>_'].value=document.<tmpl_var name>['<tmpl_var field_name>'].value">
	  	<input type=hidden name="<tmpl_var field_name>_">
	  </tmpl_unless>
	  </tmpl_if>

	  <tmpl_if url><a href="<tmpl_var url>" <tmpl_if tooltip><tmpl_var tooltip><tmpl_else> title="Click to view details for '<tmpl_var label>'" </tmpl_if> target=_new></tmpl_if>
	  <tmpl_if mouseover><a href="JavaScript:alert('<tmpl_var label> : <tmpl_var mouseover>');"  style="text-decoration: none; cursor:help;" title="<tmpl_var mouseover>" OnMouseOver="window.status='<tmpl_var mouseover>'; return true;" OnMouseOut="window.status=''; return true;">
	  </tmpl_if>
	  <tmpl_if highlight><b style="background-color:<tmpl_var highlight>"></tmpl_if><tmpl_unless active><strike></tmpl_unless><font color="<tmpl_var color>"><tmpl_var label></font><tmpl_unless active></strike></tmpl_unless><tmpl_if highlight></b></tmpl_if><tmpl_if url></a></tmpl_if>
	  <tmpl_if mouseover><img src="<tmpl_var image_path>tip.png" width="14" height="16" border="0" alt="<tmpl_var mouseover>"></a></tmpl_if>
	  <tmpl_if comment><b>(<tmpl_var comment>)</b></tmpl_if>
          </span>
       </div>
       <tmpl_if has_children>
          <div id="<tmpl_var name>-desc-<tmpl_var id>" class="hpts-block" <tmpl_if open>style="display: block"<tmpl_else>style="display: none"</tmpl_if>>
       </tmpl_if>
    <tmpl_else>
      </div>
    </tmpl_unless>
  </tmpl_loop>

  </div>
</div>

END

### Process file having tree specification
sub process_tree_file
{
	my $in = shift;
	my $separator = shift || "\t";

	open(FFF, $in);
	my @lines = <FFF>;
	close(FFF);
	chomp(@lines);

	return process_arrayref([ map { [split($separator, $_)] } grep { ! /^\s*#/ } @lines ]);
}

### Recursively process path and create XML-structure
sub process_arrayref
{
	my $in = shift;

	my $data = {};

	for my $i (@$in) {
		my @iii = @$i;
		my $head = shift(@iii);
		my @a = grep {$_} split('/', $head);
		my $a = join('', map { "{'" . $_ . "'}" } @a);
		my %prop = @iii;
		if (%prop) {
			while (my ($k,$v) = each %prop) {
				my $x = $v;
				$x =~ s/\\/\\\\/g;
				$x =~ s/'/\\'/g;
				#print "<hr><pre>v=$v</pre><hr><pre>x=$x</pre><hr>\n";
				eval "\$data->$a" . "{$k} = '$x'";
			}
		} else {
			eval "\$data->$a = {}";
		}
	}

	my $res;
	my @kkk = keys %$data;
	if (1 == @kkk) {
		$res = process_hashref($kkk[0], $data->{$kkk[0]});
	} else {
		$res = process_hashref('ROOT', $data);
	}
	return $res;
}

## Recursively take hashref and construct XML-structure (merely transforms the raw structure created in process_path_arrayref)
sub process_hashref
{
	my $k = shift;
	my $v = shift;

	my $res = { label => $k, children => [] };
	map {
		if (ref($v->{$_}) eq 'HASH') {
			push(@{$res->{children}}, process_hashref($_, $v->{$_}));
		} else {
			$res->{$_} = $v->{$_};
		}
	# } sort {$a <=> $b} keys %$v;
	} sort { (ref($v->{$a}) eq 'HASH' and defined($v->{$a}{weight}) and ref($v->{$b}) eq 'HASH' and defined($v->{$b}{weight})) ? $v->{$a}{weight} <=> $v->{$b}{weight} : $a cmp $b } keys %$v;
	return $res;
}

=head1 AUTHOR

Ramana Mokkapati, <mvr707@yahoo.com> 10 May 2004

I have been using HTML tables for structuring HTML presentation.
After seeing HTML::PopupTreeSelect from Sam Tregar <sam@tregar.com>
I liked the idea of stylesheets to indent HTML and adapted the same.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=cut

1;
