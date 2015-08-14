#!/usr/bin/perl
use strict;


#
#
#   Use to Merge strain specific snpsInProbes.txt files based on strain list
#   Then outputs Mask file
#
# Example
# perl CreateMask.pl Affy /Volumes/Data/mm10/array/430v2/Aligned/Mouse430_2.default.perfectMatches.txt ILS,ISS /Volumes/Data/mm10/array/430v2/SNPs_Probes /Volumes/Data/mm10/array/430v2/Masks/LXS.Aug2013.430v2.mask.txt
#
# This merges Probes with SNPS in ILS and ISS strains (/Volumes/Data/mm10/array/430v2/SNPs_Probes/ILS_snpsInProbes.txt)
# And reads in the 

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


my $type=$ARGV[0]; #type AffyExon Affy etc.
my $multimatch=$ARGV[1]; #file for multiple matches
my $perfectmatch=$ARGV[2]; #file for perfect matches
my $probeFA=$ARGV[3];   #file for all probes
my $snpStrain=$ARGV[4]; #a file prefix for each strain to mask probes with snps
my $snpDir=$ARGV[5]; #path to directory with all stains of snpsInProbes.txt
my $output=$ARGV[6]; #path file name for output file







#create a list of Probes with SNPS
my @snpList=split(/,/,$snpStrain);
#read in all snpsInProbes.txt files and create a masked and non masked list
my %maskedIDHOH;


my %completeIDHOH;


my %noHitIDHOH;



#read in probe fasta file to make a list of all probes
    print "try to open: $probeFA\n";
    open(IN, "<".$probeFA);
    while(<IN>){
        my $line=$_;
        $line=trim($line);
        if(index($line,">")==0){
            my @col=split(/;/,$line);
            my @col2=split(/:/,$col[0]);
            my $probeid=$col2[2].":".$col2[3].":".$col2[4];
            $completeIDHOH{$probeid}={};
            $noHitIDHOH{$probeid}=1;
        }
    }
    close IN;
    
    my @tmpList=keys %completeIDHOH;
    print @tmpList." Total Probes\n";

#read in the multiMatch file and remove any multiple alignments from the complete list
    my %multiHOH;
    my @multiList;
    my $confirm=0;
    my $totalMulti=0;
    print "try to open: $multimatch\n";
    open(IN, "<".$multimatch);
    while(<IN>){
        my $line=$_;
        $line=trim($line);
        if(!exists($multiHOH{$line}) and index($line,"AFFX")==-1){
            $totalMulti++;
            #add to multi list
            $multiHOH{$line}=1;
            push(@multiList,$line);
            #remove from probelist
            delete $completeIDHOH{$line};
            delete $noHitIDHOH{$line};
            if(exists($maskedIDHOH{$line})){
                $confirm++;
            }else{
                $maskedIDHOH{$line}=1;
            }
        }
    }
    close IN;
    my @tmpMask=keys %maskedIDHOH;
    @tmpList=keys %completeIDHOH;
    print "\n\nAfter Multi Hits Masked:".@tmpMask."\n";
    print "After Multi Hits Non-Masked:".@tmpList."\n";

#read in the perfectMatch file and remove any multiple alignments from the nonmaskedHOH list
    print "try to open: $perfectmatch\n";
    open(IN, "<".$perfectmatch);
    while(<IN>){
        my $line=$_;
        $line=trim($line);
        my @col=split(/\t/,$line);
        if(exists($noHitIDHOH{$col[3]})){
            delete $noHitIDHOH{$col[3]};
        }
        $completeIDHOH{$col[3]}{chr}=$col[0];
        $completeIDHOH{$col[3]}{start}=$col[1];
        $completeIDHOH{$col[3]}{stop}=$col[2];
    }
    close IN;
    
    my @tmpnohit=keys %noHitIDHOH;
    print @tmpnohit." no hits to genome for alignment so masking.\n";
    
    foreach my $nohit(@tmpnohit){
        if(exists($completeIDHOH{$nohit})){#remove masked entry
                    delete $completeIDHOH{$nohit};
        }
        #continue setting up mask entry
        if(!exists($maskedIDHOH{$nohit}) and index($nohit,"AFFX")==-1){
            $maskedIDHOH{$nohit}=1;
        }
    }
    
    
    @tmpMask=keys %maskedIDHOH;
    @tmpList=keys %completeIDHOH;
    print "\n\nAfter No Hits Masked:".@tmpMask."\n";
    print "After No Hits Non-Masked:".@tmpList."\n";



foreach my $snpPrefix (@snpList) {
    open(IN, "<".$snpDir."/".$snpPrefix."_snpsInProbes.txt");
    #print "open:".$snpDir."/".$snpPrefix."_snpsInProbes.txt\n";
    my $masked=0;
    my $total=0;
    while(<IN>){
        my $line=$_;
        if(index($line,"random")==-1 and index($line,"chrM")==-1 and index($line,"chrUn")==-1){
            my @col=split(/\t/,$line);
            $total++;            
            if($col[4]==0||index($col[3],"AFFX")>-1){#not masked
            }else{#masked
                $masked++;
                if(exists($completeIDHOH{$col[3]})){#remove masked entry
                    delete $completeIDHOH{$col[3]};
                }
                #continue setting up mask entry
                if(!exists($maskedIDHOH{$col[3]})){
                    
                    $maskedIDHOH{$col[3]}=1;
                }
            }
        }
    }
    print "$masked snp masked of $total in $snpPrefix\n";
    @tmpMask=keys %maskedIDHOH;
    @tmpList=keys %completeIDHOH;
    print "Masked:".@tmpMask."\n";
    print "Non-Masked:".@tmpList."\n";
    close IN;
}

#need to remove probesets with <4 probes
if($type eq "Affy"){
    my %maskedPS;
    my %probesets;
    my @tmpIDs=keys %completeIDHOH;
    foreach my $tmpID(@tmpIDs){
        my @col=split(/:/,$tmpID);
        if(!exists($probesets{$col[0]})){
            $probesets{$col[0]}=1;    
        }else{
            $probesets{$col[0]}=$probesets{$col[0]}+1;
        }
    }
    my $psCount=0;
    my @psIDs=keys %probesets;
    foreach my $psID(@psIDs){
        if($probesets{$psID}<4 and index($psID,"AFFX")==-1){
            $maskedPS{$psID}=1;
            $psCount++;
        }
    }
    my $pCount=0;
    foreach my $tmpID(@tmpIDs){
        my @col=split(/:/,$tmpID);
        if(exists($maskedPS{$col[0]})){
            $maskedIDHOH{$tmpID}=1;
            delete $completeIDHOH{$tmpID};
            $pCount++;
        }
    }
    print "\n$psCount probesets removed w/ <4 probes\n$pCount probes removed after masking probe sets\n\n";
}else{
    
}
    
    @tmpMask=keys %maskedIDHOH;
    @tmpList=keys %completeIDHOH;
    print "After masking Probe sets<4 probes Masked:".@tmpMask."\n";
    print "After masking Probe sets<4 probes Non-Masked:".@tmpList."\n";
    
    #need to remove mismatched probes for each perfect match probe
    if($type eq "Affy"){
        @tmpMask=keys %maskedIDHOH;
        foreach my $maskID(@tmpMask){
            my @col=split(/:/,$maskID);
            my $mmID=$col[0].":".$col[1].":".($col[2]+1);
            if(!exists($maskedIDHOH{$mmID}) and index($mmID,"AFFX")==-1){
                $maskedIDHOH{$mmID}=1;
            }
        }
    }
    
    @tmpMask=keys %maskedIDHOH;
    @tmpList=keys %completeIDHOH;
    my @sorttmpList=sort @tmpList;
    print "\n\nTotal Masked:".@tmpMask."\n";
    print "Total Non-Masked:".@tmpList."\n";
    
    #output Mask based on array type and list of probe locations for non masked probes
    if($type eq "Affy"){
        open(OUT,">".$output);
        foreach my $tmpID (@tmpMask){
            my @col=split(/:/,$tmpID);
                print OUT "$col[0]\t$col[1]\t$col[2]\n";
        }
        close(OUT);
        open(OUT,">".$output.".location.txt");
        foreach my $tmpID (@sorttmpList){
            my @col=split(/:/,$tmpID);
            print OUT $completeIDHOH{$tmpID}{chr}."\t".$completeIDHOH{$tmpID}{start}."\t".$completeIDHOH{$tmpID}{stop}."\t$tmpID\n";
        }
        close(OUT);
    }