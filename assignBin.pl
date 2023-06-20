#!/usr/bin/env perl

use strict;
use warnings;

my $kraken2out = $ARGV[0];
my $ncbiTaxDir = $ARGV[1];
my $minFraction = $ARGV[2];
my $minNumber = $ARGV[3];
my $outfile = $ARGV[4];
my ($tax2parent, $merged, $tax2rank, $tax2name, $deleted) = load_ncbi_taxonomy($ncbiTaxDir);

my $assignmentCounts = {};

open (KRKN, '<', $kraken2out) or die $!;
while (my $line = <KRKN>) {
    my @elems = split(/\t/, $line);
    # kraken uses 'C', metabuli '1' for classified sequences
    if (($elems[0] eq 'C') || ($elems[0] eq '1')) {
        my $taxid = int($elems[2]);
        $taxid = $merged->{$taxid} if defined($merged->{$taxid});
        if (!defined($deleted->{$taxid})) {
            $assignmentCounts->{$taxid}++;
            while ($taxid != $tax2parent->{$taxid}) {
                #die "no parent for $taxid" unless defined($tax2parent->{$taxid});
                $taxid = $tax2parent->{$taxid};
                $taxid = $merged->{$taxid} if defined($merged->{$taxid});
                $assignmentCounts->{$taxid}++;
            }
        }
    }
}
close(KRKN);

foreach my $rank (qw(species genus family order class phylum superkingdom)) {
    my $rankProfile = {};
    my $max=0;
    my $sum=0;
    my $maxTaxId;
    foreach my $tid (keys %$assignmentCounts) {
        # extract tax. profile at given rank
        die "no rank for $tid" unless defined($tax2rank->{$tid});
        if ($tax2rank->{$tid} eq $rank and $assignmentCounts->{$tid} >= $minNumber) {
            $rankProfile->{$tid} = $assignmentCounts->{$tid};
            if ($rankProfile->{$tid} > $max) {
                $max = $rankProfile->{$tid};
                $maxTaxId = $tid; 
            }
            $sum += $rankProfile->{$tid};
        }
    }

    #print STDERR "profile at $rank has ".$sum." entries, dominant fraction is ".$max/$sum."\n"; 

    if ($sum > 0 and $max/$sum >= $minFraction) {
        #print "assigment at rank ".$rank."\n";
        my $lin = lineage($maxTaxId, $tax2parent,$merged, $tax2name, $tax2rank);
        chop($lin); chop($lin);
        open(OUT, '>', $outfile) or die $!;
        print OUT $lin . "\n";
        close(OUT);
        exit 0;
    }
}

open(OUT, '>', $outfile) or die $!;
print OUT "root\n";
close(OUT);

sub isMajorRank {
    my ($tid, $tax2rank) = @_;
    foreach my $rank (qw(species genus family order class phylum superkingdom)) {
        return 1 if $tax2rank->{$tid} eq $rank;
    }
    return 0;
}

sub lineage {
    my ($tid, $t2p, $merged, $tax2name, $tax2rank) = @_;
    if ($tid == 1) {
        return $tax2name->{$tid}.'; ';
    }
    $tid = $merged->{$tid} if exists($merged->{$tid});
    my $parent = $t2p->{$tid};
    die "no parent for $tid\n" unless defined($parent);
    return undef unless defined($parent);
    my $lin = lineage($parent, $t2p, $merged, $tax2name, $tax2rank);
    return undef unless defined($lin);
    return $lin . $tax2name->{$tid}.'; ' if isMajorRank($tid, $tax2rank);
    return $lin;
}

sub load_ncbi_taxonomy {
    my $taxDir = shift;
    my $data = {};
    my $tax2rank = {};
    my ($taxid, $parent, $rank);
    open(my $namesfile, '<', $taxDir.'/nodes.dmp') or die 'Cannot open ".$taxDir."/nodes.dmp';
    while (<$namesfile>) {
        chop; chop; chop;  # remove '\t|\n'
        ($taxid, $parent, $rank) = split(/\t\|\t/o, $_);
        $data->{int($taxid)} = int($parent);
        $tax2rank->{int($taxid)} = $rank;
    }
    close($namesfile);


    my $merged = {};
    open(my $mergedfile, '<', $taxDir.'/merged.dmp') or die 'Cannot open ".$taxDir."/merged.dmp';
    my $renamed;
    while (<$mergedfile>) {
        chop; chop; chop;  # remove '\t|\n'
        ($taxid, $renamed) = split(/\t\|\t/o, $_);
        $merged->{$taxid} = $renamed;
    }
    close($mergedfile);

    my $tax2name = {};
    open(my $names, '<', $taxDir.'/names.dmp') or die 'Cannot open ".$taxDir."/names.dmp';
    while (<$names>) {
        chop; chop; chop;  # remove '\t|\n'
        my @elems = split(/\t\|\t/o, $_);
        if ($elems[3] eq 'scientific name') {
            $tax2name->{$elems[0]} = $elems[1];
        }
    }
    close($names);

    my $deleted = {};
    open(my $delnodes, '<', $taxDir.'/delnodes.dmp') or die 'Cannot open ".$taxDir."/delnodes.dmp';
    while (<$delnodes>) {
        chop; chop; chop;  # remove '\t|\n'
        $deleted->{$_} = 1;
    }
    close($delnodes);


    return ($data, $merged, $tax2rank, $tax2name, $deleted);
}

