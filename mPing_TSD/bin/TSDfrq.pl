#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt, "input:s" ,"help");


my $help=<<USAGE;
perl $0 --input ../input/mping_flank_600.fa

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

getfastalen($opt{input});

##################################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}

##
#>mping.Chr10_21829433_21829435 mping.rc mping A119_2 Chr10:21829433..21829435 FLANK1:1..600 TSD1:598..600 TE:601..1030 TSD2:1031..1033 FLANK2:1031..1630
#AACTTTCTCCATTTCTCGATGACCGATTCATCGTTGACAAGTGTTTATTACAGTAGAAAAGAGACAAGTTCTTCAGCGTGCATCACGGTTAGTACTCCCATTTATGACCTCAAGCCAGTCGCATCTAAGTTCGATAGGCGTGCAGACACTGAAGTTAACCTT
##
sub getfastalen
{
$/=">";
my %hash;
my $total;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$temp\n";
    my ($strand, $tsds1, $tsde1, $tsds2, $tsde2);
    if ($temp =~/ mping\.(rc|fwd) .*TSD1:(\d+)\.\.(\d+) .*TSD2:(\d+)\.\.(\d+)/){
       $strand = $1 eq 'fwd'  ? 1 : 0;
       $tsds1 = $2;
       $tsde1 = $3;
       $tsds2 = $4;
       $tsde2 = $5;
       $l = $tsde1 - $tsds1 + 1;
       my $tsd = $strand == 1 ? substr($seq, $tsds1-1, $l) : revcom(substr($seq, $tsds1-1, $l));
       $tsd = 'TTA' if $tsd eq 'TAA';
       $tsd = 'TCA' if $tsd eq 'TGA';
       $hash{$tsd} += 1;
       $total += 1;
       #print "$strand\t$tsd\n"; 
    }
    #print "$strand\t$tsds1\t$tsde1\t$tsds2\t$tsde2\n";
}
close IN;
$/="\n";
for (keys %hash){
    my $frq= $hash{$_}/$total;
    my $rev= revcom($_);
    print "$_\t$rev\t$hash{$_}\t$frq\n";
}
return \%hash;
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}

 
