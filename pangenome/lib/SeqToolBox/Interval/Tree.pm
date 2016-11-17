# $Id: Tree.pm 550 2010-04-30 17:50:48Z malay $
# Perl module for SeqToolBox::Interval::Tree
# Author: Malay <malaykbasu@gmail.com>
# Copyright (c) 2009 by Malay. All rights reserved.
# You may distribute this module under the same terms as perl itself


##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##


=head1 NAME

SeqToolBox::Interval::Tree  - DESCRIPTION of Object

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


package SeqToolBox::Interval::Tree;

use vars qw(@ISA);
@ISA = qw();
@EXPORT_OK = qw();
use strict;


##-------------------------------------------------------------------------##
## Constructors
##-------------------------------------------------------------------------##

=head1 CONSTRUCTOR

=head2 new()

=cut

sub new {
    my $args  = shift;
    my $class = ref $args || $args;
    my $self  = {};
    bless $self, $class;
    $self->_init();
    return $self;

}

sub _init {
    my $self = shift;
    $self->{ranges} = ();
}



##-------------------------------------------------------------------------##
## METHODS
##-------------------------------------------------------------------------##


=head1 PUBLIC METHODS

=cut

sub getroot {
    my $self = shift;
    if ( defined $self->{root} ) {
        return $self->{root};
    }
    else {
        return undef;
    }
}

sub setroot {
    my ( $self, $node ) = @_;

    delete( $node->{left} );
    $self->{root} = $node;

#print "Current root: ", $self->getroot()->to_string, "and it's left is: ", $self->getroot()->{left}, "\n";
}

#sub insert {
#    my ($self, $range) = @_;
#    if ($self->{ranges}) {

#    my @array = @{$self->{ranges}};
#    $array[@array] = $range;
#    $self->{ranges} = \@array;
#}
#    else {
#	my @array;
#	$array[0] = $range;
#	$self->{ranges} = \@array;
#    }
#}

sub get_tree {
    my $self = shift;

    # print "print tree\n";
    my $root = $self->getroot;

    # print "Consesus Tree:";
    my $s = "";
    while ($root) {

        #print "loop\n";
        $s .= $root->{name} . '[' . $root->{start} . '..' . $root->{end} . '];';

        $root = $root->right();

        # print "end\n";
    }
    $s =~ s/\;$//;
    return $s;

    #print "print tree end\n";
}

sub get_names {
    my $self = shift;

    # print "print tree\n";
    my $root = $self->getroot;

    # print "Consesus Tree:";
    my @s;
    while ($root) {

        #print "loop\n";
        $s[@s] = $root->{name};

        $root = $root->right();

        # print "end\n";
    }

    return @s;

    #print "print tree end\n";
}

sub name_exists {
	my ($self, $name) = @_;
	my $root = $self->getroot;
	while ($root) {
		if ($root->{name} eq $name) {
			return 1;
		}
		$root = $root->right();
	}
	return 0;
}

sub get_node_by_name {
	my ($self, $name) = @_;
	my $root = $self->getroot;
	while ($root) {
		if ($root->{name} eq $name) {
			return $root;
		}
		$root = $root->right();
	}
	return;
}

#sub get_tree {
#    my $self = shift;
##    my $tree = Tree->new();
##    foreach my $node(@{$self->{ranges}}) {
##	$tree->insert_node($node);
##    }
##    $tree->print_tree();
##}

sub insert {

    my ( $self, $range ) = @_;

    #$self->print_tree();
    #print "\n";

    print STDERR "inserting " . $range->to_string() . "\n";

    if ( !$self->getroot() ) {
        $self->setroot($range);
        $range->setasroot($self);
        return $self;
    }

    my $parent = undef;

    my $node = $self->getroot;

    # my $newrange = Range->new();
    my @intersected;

    while ($node) {
        $parent = $node;

        #print "while loop\n";

        if ( $node->intersects($range) ) {

            #$node->merge($range);
            # $newrange->set_start($node->start());
            #	    $newrange->set_end($node->end());
            #	    $range = $newrange;
            #$range = $node;
            #$node= $self->getroot;
            #print $range->to_string().' intersects '.$node->to_string."\n";
            $intersected[@intersected] = $node;

            #print "Right node:", $node->right()->to_string;
            $node = $node->right();

        }
        else {

            if ( $range->proximalto($node) ) {

                #$node = $node->left();
                # next;
                #print "node proximal\n";
                last;
            }
            else {
                $node = $node->right();

                #next;
                #print "node right\n";
            }
        }
    }
    if ( !@intersected ) {

        #print "not intersected\n";
        #print "Parent: ".$parent->to_string , "\n";
        #print $parent->left();
        if ( $range->proximalto($parent) ) {

            #print "Proximal to paret\n";
            if ( $parent->isroot() ) {

                $self->setroot($range);
                $range->setasroot($self);
                $parent->{isroot} = 0;

                #print "Inserted ",$range->to_string(),"as root\n";
            }

            my $leftnode = $parent->left();
            $range->{right}    = $parent;
            $parent->{left}    = $range;
            $range->{left}     = $leftnode;
            $leftnode->{right} = $range;

            #print "link set\n";

        }
        else {
            my $rightnode = $parent->right();
            if ($rightnode) {
                $range->{right}    = $rightnode;
                $rightnode->{left} = $range;
            }

            $parent->{right} = $range;
            $range->{left}   = $parent;
        }
    }
    else {

        #$self->rearrange( $range, @intersected );
        #return 1;
    }

}

sub rearrange {
    print STDERR "Rearrage called\n";
    my ( $self, $range, @nodes ) = @_;
    my $first_node = $nodes[0];
    my $last_node  = $nodes[$#nodes];

    #    my $left_start;
    #    my $right_start;
    if ( @nodes > 1 ) {
        if ( $range->{name} =~ /^smart/i ) {

            $first_node->replace_with( $range, $self );
            $range->{right} = $last_node->{right};

            #$first_node->{end}   = $last_node->{end};
            #$first_node->{right} = $last_node->{right};
        }
    }
    else {
        print STDERR "One node intersected\n";
        if ( $nodes[0]->{name} =~ /^smart/i ) {

        }

        else {
            print STDERR "new node is smart\n";

            $nodes[0]->replace_with( $range, $self );
        }
    }
}

sub print_tree {
    print STDERR "Print tree called\n";
    my $self = shift;

    # print "print tree\n";
    my $root = $self->getroot;
    print "Root", $root->to_string(), "\n";
    print $root->{right}->to_string(), "\n";

    # print "Consesus Tree:";
    while ($root) {

        #print "loop\n";
        print $root->{start} . '..' . $root->{end} . ';';

        $root = $root->right();

        # print "end\n";
    }

    #print "print tree end\n";
}




=head1 PRIVATE METHODS

=cut




=head1 SEE ALSO

=head1 COPYRIGHTS

Copyright (c) 2009 by Malay <malaykbasu@gmail.com>. All rights reserved.
This program is free software; you can redistribute it and/or modify it under
the same terms as Perl itself.

=cut


=head1 APPENDIX

=cut

1;