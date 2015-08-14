#!/usr/bin/perl
use strict;


my $input=$ARGV[0];
my $strain=$ARGV[2];
my $output=$ARGV[1];
my $incolumns=$ARGV[3];


open IN, "<",$input or die "could not open input $!";
open OUT, ">",$output or die "could not open output $!";

my %columns;

my @printColumns=split(/,/,$incolumns);
print "$incolumns\n";


my $count=0;
while(<IN>){
    my $line=$_;
    if(index($line,"##")==0){
        
    }elsif(index($line,"#CHROM")==0){
        
        my @col=split(/\t/,$line);
        my $curInd=0;
        foreach my $colName(@col){
            $columns{$colName}=$curInd;
            $curInd++;
        }
    }else{
        my @col=split(/\t/,$line);
        if(index($col[$columns{$strain}],"1/1:")==0 and $col[$columns{"FILTER"}] eq "PASS"){
            foreach my $curCol(@printColumns){
                if($curCol eq "POS"){
                    #my $posPlus=$col[$columns{$curCol}]+length($col[$columns{"ALT"}]);
                    my $posPlus=$col[$columns{$curCol}]+1;
                    print OUT $col[$columns{$curCol}]."\t".$posPlus."\t";
                }elsif($curCol eq "CHROM"){
                    print OUT "chr".$col[$columns{$curCol}]."\t";
                }else{
                    print OUT $col[$columns{$curCol}]."\t";
                }
            }
            print OUT $col[$columns{$strain}]."\n";
        }
    }
    $count++;
}

close IN;
close OUT;