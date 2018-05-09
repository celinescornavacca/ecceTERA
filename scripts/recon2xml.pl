#!/usr/bin/perl

# Convert a TERA reconciliation file to phyloxml format.
# Branch lengths can be read if the TERA gene newick file is given
# as the second argument.

use strict;

my $spacing = "  ";


my $reconFile = shift @ARGV;
if( !$reconFile || $reconFile eq '-h' ) {
    print "USAGE: <reconFile> [geneNewickFile]$/";
    print "Convert a TERA reconciliation file to phyloxml format.$/";
    print "Branch lengths can be read if the TERA gene newick file is given$/";
    print "as the second argument.$/";
    exit(1);
}

my ($root,$eventListsRef,$geneParentsRef,$geneChildrenRef) 
        = readRecon( $reconFile );

my %gBranchLengths;
if( @ARGV ) {
    my $newickFile = shift @ARGV;
    readBranchLengthsFromNewick( $newickFile );
}

open XML, ">$reconFile.xml" or die "Could not create xml file";

    print XML <<EOF;
<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
<phylogeny rooted="true" rerootable="true">
EOF

printClades( $root, $spacing );
print XML "</phylogeny></phyloxml>$/";
close XML;

exit(1);

sub printClades {
    my ($id, $spaces) = @_;

    my $sp = $spaces . $spacing;
    print XML "$sp<clade>$/";
    if( !$geneChildrenRef->{$id} ) {
        print XML "$sp$spacing<name>$id<\/name>$/";
    }
    if( $gBranchLengths{$id} ) {
        print XML "$sp$spacing<branch_length>$gBranchLengths{$id}"
                  . "<\/branch_length>$/";
    }
#print XML "$sp$spacing<confidence type=\"bootstrap\">100<\/confidence>$/";

    my $reverseClades = 0; # switch order of clades if true
    foreach my $event (@{$eventListsRef->{$id}}) {
        $reverseClades = printEvent( $event, $sp.$spacing );
    }

    if( $geneChildrenRef->{$id} ) {
        my @clades = @{$geneChildrenRef->{$id}};
        if( $reverseClades ) {
            @clades = reverse @clades;
        }
        foreach my $childId (@clades) {
            printClades( $childId, $sp );
        }
    }
    print XML "$sp<\/clade>$/";
}

sub printEvent {
    my ($event, $spaces) = @_;

    my $reverseClades = 0;
    print XML "$spaces<rec:event>$/";
    if( $event =~ /(.*),([ST]L\w*),(.*)@(.*)/ ) {
        my $lossLoc = $3;
        my $name = $2;
        if( $lossLoc =~ /(.*),(.*)/ ) {
            $lossLoc = $1; # new format
        }
        if( $name eq 'SL' ) {
            printSpeciation( $1, $spaces.$spacing );
        } elsif( $name eq 'TL' ) {
            $lossLoc = $1;
            printTransfer( $1, $lossLoc, $spaces.$spacing );
        } elsif( $name eq 'TLFD' ) {
            $lossLoc = 'dead';
            printTransfer( 'dead', $lossLoc, $spaces.$spacing );
        } elsif( $name eq 'TLTD' ) {
            $lossLoc = $1;
            printTransfer( $1, 'dead', $spaces.$spacing );
        } else {
            die "Unknown event: $event";
        }
        printLoss( $lossLoc, $spaces.$spacing )
    } elsif( $event =~ /(.*),([SDT]\w*),(.*),(.*)@(.*)/ ) {
        if( $2 eq 'S' ) {
            printSpeciation( $1, $spaces.$spacing );
        } elsif( $2 eq 'D' ) {
            printDuplication( $1, $spaces.$spacing );
        } elsif( $2 eq 'DD' ) {
            printDuplication( 'dead', $spaces.$spacing );
        } elsif( $2 eq 'T' || $2 eq 'TFD' ) {
            my $toId = $4;
            if( $1 eq $toId ) {
                $toId = $3; 
                $reverseClades = 1;
            }
            printTransfer( $1, $toId, $spaces.$spacing );
        } elsif( $2 eq 'TTD' ) {
            printTransfer( $1, 'dead', $spaces.$spacing );
        } elsif( $2 eq 'TFD' ) {
            my $toId = $4;
            $toId = $3 if( $4 eq '-1' );
            printTransfer( 'dead', $toId, $spaces.$spacing );
        } else {
            die "Unknown event: $event";
        }
    } else {
        die "Unknown event format: $event";
    }
    print XML "$spaces<\/rec:event>$/";
    $reverseClades;
}

sub printSpeciation {
    my ($loc, $spaces) = @_;
    print XML "$spaces<rec:speciation>$/";
    print XML "$spaces$spacing<rec:species_location>$loc"
             . "<\/rec:species_location>$/";
    print XML "$spaces<\/rec:speciation>$/";
}

sub printDuplication {
    my ($loc, $spaces) = @_;
    print XML "$spaces<rec:duplication>$/";
    print XML "$spaces$spacing<rec:species_location>$loc"
             . "<\/rec:species_location>$/";
    # time slice
    print XML "$spaces<\/rec:duplication>$/";
}

sub printTransfer {
    my ($origin, $recipient, $spaces) = @_;
    print XML "$spaces<rec:transfer>$/";
    print XML "$spaces$spacing<rec:species_origin>$origin"
             . "<\/rec:species_origin>$/";
    print XML "$spaces$spacing<rec:species_recipient>$recipient"
             . "<\/rec:species_recipient>$/";
    # time slice
    print XML "$spaces<\/rec:transfer>$/";
}

sub printLoss {
    my ($loc, $spaces) = @_;
    print XML "$spaces<rec:loss>$/";
    print XML "$spaces$spacing<rec:species_location>$loc"
             . "<\/rec:species_location>$/";
    print XML "$spaces<\/rec:loss>$/";
}


sub getBranchLengths {
    my ($newickStr) = @_;

    my @pairs = $newickStr =~ /[^(),]+/g;
    foreach my $p (@pairs) {
        if( $p =~ /(.*):(.*)/ ) {
            $gBranchLengths{$1} = $2;
        } else {
            die "Problem near $p";
        }
    }
}

sub readBranchLengthsFromNewick {
    my ($newickFile) = @_;

    my $description;
    open NEWICK, $newickFile or die "Could not open $newickFile";
    while( <NEWICK> ) {
        chomp;
        next if /^\s*$/;
        if( /^(.*;)\s*$/ ) {
            $description .= $1;
            last; 
        } elsif( /^(.*)\s*$/ ) {
            $description .= $1;
        }
    }
    close NEWICK;

    if( $description != /^\((.*)\)\d+;$/ ) {
        getBranchLengths( $1 );
    } else {
        die "No newick string found in $newickFile with the form (...)#;";
    }


}

# read the TERA recon format
#55:71,S,15,70@1;70,SL,32@1;69,SL,41@1;68,S,42,67@1:11,54
sub readRecon {
    my ($reconFile) = @_;

    my $root;
    my %eventLists;
    my %geneParents;
    my %geneChildren;
    open RECON, $reconFile or die "Could not open $reconFile";
    while( <RECON> ) {
        chomp;
        next if /^\s*$/;
        if( /^(.+):(.+):(.+),(.+)\s*$/ ) {
            my $geneId = $1;
            my $eventStr = $2;
            my $childL = $3;
            my $childR = $4;
            $eventStr =~ s/'//g;
            my @events = split ';', $eventStr;
            shift @events if $root; #first event is repeated previous last
            $eventLists{$geneId} = \@events;
            push @{$geneChildren{$geneId}}, $childL;
            push @{$geneChildren{$geneId}}, $childR;
            $geneParents{$childL} = $geneId;
            $geneParents{$childR} = $geneId;
            $root = $geneId if !$root;
        } elsif( /^(.+):(.+)\s*$/ ) {
            my $geneId = $1;
            my $eventStr = $2;
            $eventStr =~ s/'//g;
            my @events = split ';', $eventStr;
            shift @events if( $root ); #first event is repeated lasted event
            $eventLists{$geneId} = \@events;
            $root = $geneId if !$root;
        } else {
            die "Bad line: $_";
        }
    }
    close RECON;

    ($root, \%eventLists, \%geneParents, \%geneChildren);
}
