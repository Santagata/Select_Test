#!/usr/bin/perl

use strict;
use warnings; 
use Data::Dumper;

MAIN: {
    my $c_file = $ARGV[0] or die "usage: $0 CDS_FILE ALN_FILE\n";
    my $a_file = $ARGV[1] or die "usage: $0 CDS_FILE ALN_FILE\n";

    my @c_files = ();
    open IN, $c_file or die "cannot open $c_file:$!";
    while (my $line = <IN>) {
        chomp $line;
        next if ($line =~ m/^\s*$/);
        $line =~ s/\s+$//;        
        push @c_files, $line;
    }

    my @a_files = ();
    open IN, $a_file or die "cannot open $a_file:$!";
    while (my $line = <IN>) {
        chomp $line;
        next if ($line =~ m/^\s*$/);
        $line =~ s/\s+$//;
        push @a_files, $line;
    }

#    my %data = ();
    for (my $i = 0; $i < @c_files; $i++) {
        print "perl pal2nal.pl $a_files[$i] $c_files[$i] -output paml -nogap -codontable 1 > $c_files[$i]_align.phy\n";
        system "perl pal2nal.pl $a_files[$i] $c_files[$i] -output paml -nogap -codontable 1 > $c_files[$i]_align.phy\n";
#        $data{$a_files[$i]} = $c_files[$i];
    }
#    print Dumper \%data;
}



