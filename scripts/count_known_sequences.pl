#!/usr/bin/env perl
use strict;
use warnings;

my ($pattern_file) = @ARGV;
die "usage: count_known_sequences.pl patterns.tsv\n" unless defined $pattern_file;

open my $pfh, "<", $pattern_file or die "cannot open $pattern_file: $!\n";
my %strategy;
my %patterns;
while (my $line = <$pfh>) {
    chomp $line;
    next if $line eq "";
    my ($deletion_id, $strategy, $group_index, $pattern) = split /\t/, $line, 4;
    next unless defined $pattern && $pattern ne "";
    $strategy{$deletion_id} = $strategy;
    push @{ $patterns{$deletion_id}{$group_index} }, $pattern;
}
close $pfh;

my %counts = map { $_ => 0 } keys %strategy;
my $line_number = 0;
while (my $line = <STDIN>) {
    ++$line_number;
    next unless $line_number % 4 == 2;
    chomp $line;
    my $seq = uc($line);
    for my $deletion_id (keys %strategy) {
        my $matched = 1;
        for my $group_index (keys %{ $patterns{$deletion_id} }) {
            my $group_matched = 0;
            for my $pattern (@{ $patterns{$deletion_id}{$group_index} }) {
                if (index($seq, $pattern) >= 0) {
                    $group_matched = 1;
                    last;
                }
            }
            if (!$group_matched) {
                $matched = 0;
                last;
            }
        }
        ++$counts{$deletion_id} if $matched;
    }
}

for my $deletion_id (sort keys %counts) {
    print "$deletion_id\t$counts{$deletion_id}\n";
}
