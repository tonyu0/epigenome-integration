#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 2) {
    die "Usage: $0 <arg1> <arg2>\n<arg1> is a genome reference file, <arg2> is a chromosome index you want)\n";
}

my ($faFile, $outChr) = @ARGV;
open(IN, '<', $faFile) or die "failed to open the specified genom reference file\n" ;
open(OUT, '>', "outChr$outChr.fa") or die "failed to open output file\n";
$\ = "\n"; # change the last char of each print to \n

my $outFlag = 0;

while(<IN>) {
    chomp;
    my $line = $_;
    # Header description is like:
    # 1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
    if($line =~/(^>)(\d+)\s+(.+?)\s+(.+?)\s+/){ 
        if ($outFlag) { # already outputting data, and reached the next chromosome
            last;
        }
        # $header = $2 ;
        # $len = length($header) if($len<length($header));

        my @infoList = split(/:/, $4);
        my $genome = $infoList[1];
        my $chrIdx = $infoList[2];
        my $chrLen = $infoList[4];
        # print OUT "Genome: $genome, ChrIdx: $chrIdx, ChrLen: $chrLen\n"; # use this line to check header info
        if($chrIdx == $outChr) {
            print OUT $line; # start output selected chromosome data
            $outFlag = 1;
        } 
    } else {
        if ($outFlag) {
            print OUT $line;
        }
    }
    # else
    # {
    #     $line =~ tr/a-z/A-Z/;
    #     $seq{$header} = length $line;
    # }
}

close(IN);
close(OUT);