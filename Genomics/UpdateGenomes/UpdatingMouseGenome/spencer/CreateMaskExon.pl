#!/usr/bin/perl
use strict;
use DateTime;


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

sub writePGF{
    my($pgfHOHRef,  $header,$output) = @_;
    my %pgfHOH=%$pgfHOHRef;
    open OUT,">".$output;
    print OUT $header;
    my @psList=keys %pgfHOH;
    
    my @sorted=sort @psList;
    foreach my $ps(@sorted){
        print OUT $pgfHOH{$ps}{ID}."\t".$pgfHOH{$ps}{cat};
        print OUT "\n";
        my $tmpRef=$pgfHOH{$ps}{probes};
        my %probeHOH=%$tmpRef;
        my @probeList=keys %$tmpRef;
        my @sorted=sort {
                            $pgfHOH{$ps}{probes}{$a}{pinfo}{atom}<=>$pgfHOH{$ps}{probes}{$b}{pinfo}{atom}
                            or
                            $pgfHOH{$ps}{probes}{$a}{pinfo}{pID}<=>$pgfHOH{$ps}{probes}{$b}{pinfo}{pID}
                         } @probeList;
        my $prevAtom="";
        foreach my $probe(@sorted){
            if($prevAtom eq $pgfHOH{$ps}{probes}{$probe}{pinfo}{atom}){
                
            }else{
                print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{atom}."\n";
                $prevAtom=$pgfHOH{$ps}{probes}{$probe}{pinfo}{atom};
            }
            print OUT "\t\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{pID};
            print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{type};
            print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{gc};
            print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{plen};
            print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{intPos};
            print OUT "\t".$pgfHOH{$ps}{probes}{$probe}{pinfo}{seq}."\n";
        }
    }
    close OUT;
}

my $arrayDir=$ARGV[0];
my $pgf=$ARGV[1]; #pgf filename Should be in /Source/name
my $snpStrain=$ARGV[2]; #a file prefix for each strain to mask probes with snps
my $psPrefix=$ARGV[3];
my $outputFileName=$ARGV[4];
my $ver=$ARGV[5];
my $psSuffix=$ARGV[6];
my $psFile=$ARGV[7];

my $arrayStr=substr($psPrefix,0,rindex($psPrefix,"."));
$arrayStr=substr($arrayStr,0,rindex($arrayStr,".")+1);

$psPrefix=$arrayDir."/Source/".$psPrefix;

my @maskList;
push(@maskList,"core");
push(@maskList,"extended");
push(@maskList,"full");
push(@maskList,"all");





#create a list of Probes with SNPS
my @snpList=split(/,/,$snpStrain);
#read in all snpsInProbes.txt files and create a masked and non masked list

my %maskedP;

my @maskedPS;

my %uniqueP;

my %pgfHOH;
my %pToPS;
my %psToTc;
my %uniqueTC;

my %noHitP;
my $header="";

#read in csv for psToTc
open (IN,"<",$arrayDir."/Source/".$psFile);
print "reading:".$arrayDir."/Source/".$psFile."\n";
while(<IN>){
    my $line=$_;
    $line =~ s/"//g;
    if(index($line,"#")==0 or index($line,"probeset_id")==0){#header line
    }else{
        my @col=split(/,/,$line);
        if(!($col[0] eq "") and !($col[3] eq "") ){
            $psToTc{$col[0]}=$col[3];
            $uniqueTC{$col[3]}{tcprobes}{$col[0]}={};
        }
    }
}

my @tmpTCList=keys %uniqueTC;
my @tmpPSTCList=keys %psToTc;
print "TC Size:".@tmpTCList."\n";
print "CSV PS Size:".@tmpPSTCList."\n";
close IN;
print "read in csv\n";

#read in pgf to make a list of all probes
    print "try to open: $pgf\n";
    open(IN, "<".$arrayDir."/Source/".$pgf);
    my $atom;
    my $ps="";
    my $count=0;
    print "Processing PGF";
    while(<IN>){
        my $line=$_;
        if(index($line,"\t\t")==0){#probe line
            my @col=split(/\t/,$line);
            $pToPS{$col[2]}=$ps;
            my $seq=$col[7];
            $seq=trim($seq);
            $seq=~ s/\r//g;
            #print "probe:".$col[2]."\n";
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{pID}=$col[2];
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{atom}=$atom;
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{type}=$col[3];
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{gc}=$col[4];
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{plen}=$col[5];
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{intPos}=$col[6];
            $pgfHOH{$ps}{probes}{$col[2]}{pinfo}{seq}=$seq;
        }elsif(index($line,"\t")==0){#atom id line
            my @col=split(/\t/,$line);
            $atom=$col[1];
            $atom=trim($atom);
            $atom =~ s/\r//g;
            #print "atom:".$atom."\n";
        }elsif(index($line,"#")==0){#header line
            $line =~ s/\r//g;
            $header=$header.$line;
        }else{##probeset line
            $count++;
            my @col=split(/\t/,$line);
            $ps=$col[0];
            my $cat=trim($col[1]);
            $cat=~ s/\r//g;
            #print "probe set:".$ps."\n";
            $pgfHOH{$col[0]}{ID}=$col[0];
            $pgfHOH{$col[0]}{cat}=$cat;
            $pgfHOH{$col[0]}{probes}={};
            if($count%50000==0){
                print ".";
            }
        }
        
    }
    close IN;
    
    my @tmpList=keys %pToPS;
    print "\nInitial ".@tmpList." Total Probes\n";
    @tmpList=keys %pgfHOH;
    print "Initial ".@tmpList." Total Probe Sets\n";


foreach my $snpPrefix (@snpList) {
    open(IN, "<".$arrayDir."/SNPs_Probes/".$snpPrefix."_snpsInProbes.txt");
    my $masked=0;
    my $total=0;
    while(<IN>){
        my $line=$_;
        if(index($line,"random")==-1 and index($line,"chrM")==-1 and index($line,"chrUn")==-1){
            my @col=split(/\t/,$line);
            $total++;            
            if($col[6]==0){
                if(!exists($maskedP{$col[3]})){#not masked
                    $uniqueP{$col[3]}=1;
                    $pgfHOH{$pToPS{$col[3]}}{probes}{$col[3]}{pinfo}{chr}=$col[0];
                    $pgfHOH{$pToPS{$col[3]}}{probes}{$col[3]}{pinfo}{start}=$col[1];
                    $pgfHOH{$pToPS{$col[3]}}{probes}{$col[3]}{pinfo}{stop}=$col[2];
                    $pgfHOH{$pToPS{$col[3]}}{probes}{$col[3]}{pinfo}{strand}=$col[5];
                }
            }else{#masked
                if(!exists($maskedP{$col[3]}) and ($pgfHOH{$pToPS{$col[3]}}{cat} eq "main" or index($pgfHOH{$pToPS{$col[3]}}{cat},"normgene")==0)){
                    $masked++;
                    if(exists($uniqueP{$col[3]})){#remove masked entry
                        delete $uniqueP{$col[3]};
                    }
                    #continue setting up mask entry
                    $maskedP{$col[3]}=1;
                    delete $pgfHOH{$pToPS{$col[3]}}{probes}{$col[3]};
                    delete $pToPS{$col[3]};
                }
            }
        }
    }
    print "$masked snp masked of $total in $snpPrefix\n";
    my @tmpMask=keys %maskedP;
    my @tmpList=keys %pToPS;
    print "Masked probes:".@tmpMask."\n";
    print "Non-Masked probes:".@tmpList."\n";
    close IN;
}



my @tmpMask=keys %uniqueP;
print "\nUnique location list size: ".@tmpMask."\n";

@tmpList=keys %pgfHOH;
my $psRemoveCount=0;
my $pRemoveCount=0;
my $pNonUnique=0;
my $procCount=0;

my $psRemChr=0;
my $psRemPloc=0;

print "Removing PS/NonUnique";
foreach my $tmpProbeSet (@tmpList) {
    if($pgfHOH{$tmpProbeSet}{cat} eq "main" or index($pgfHOH{$tmpProbeSet}{cat},"normgene")==0){
        #print "probeset:".$tmpProbeSet."\n";
        $procCount++;
        my $tmpRef=$pgfHOH{$tmpProbeSet}{probes};
        my @probeList=keys %$tmpRef;
        my $priorCount=@probeList;
        my $curRem=0;
        my $locationDisagreement=0;
        my $tmpChr="unset";
        my $tmpStrand="";
        my $min=1999999999;
        my $max=0;
        foreach my $tmpProbe (@probeList){
            if(!exists $uniqueP{$tmpProbe}){
                #print "DNE\n";
                $curRem++;
                $pNonUnique++;
                $maskedP{$tmpProbe}=1;
                delete $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe};
            }elsif(exists $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}){
                if($pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{start}<$min){
                    $min=$pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{start};
                }
                if($pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{stop}>$max){
                    $max=$pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{stop};
                }
                if($tmpChr eq $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{chr}
                   and
                   $tmpStrand eq $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{strand}){
                    
                }else{
                    if($tmpChr eq "unset"){
                        if(!(exists $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{chr}) or $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{chr} eq ""
                           or
                           !(exists $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{strand}) or $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{strand} eq ""){
                            
                        }else{
                            $tmpChr=$pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{chr};
                            $tmpStrand=$pgfHOH{$tmpProbeSet}{probes}{$tmpProbe}{pinfo}{strand};
                        }
                    }else{
                        $locationDisagreement=1;
                    }
                }
            }
        }
        $tmpRef=$pgfHOH{$tmpProbeSet}{probes};
        @probeList=keys %$tmpRef;
        if($max-$min>15000){
            $psRemoveCount++;
            $psRemPloc++;
            foreach my $tmpProbe(@probeList){
                $pRemoveCount++;
                $maskedP{$tmpProbe}=1;
                delete $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe};
            }
            delete $pgfHOH{$tmpProbeSet};
            push (@maskedPS,$tmpProbeSet);
        }elsif($locationDisagreement==1){
            $psRemoveCount++;
            $psRemChr++;
            foreach my $tmpProbe(@probeList){
                $pRemoveCount++;
                $maskedP{$tmpProbe}=1;
                delete $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe};
            }
            delete $pgfHOH{$tmpProbeSet};
            push (@maskedPS,$tmpProbeSet);
        }elsif(@probeList<3){
            $psRemoveCount++;
            foreach my $tmpProbe(@probeList){
                $pRemoveCount++;
                $maskedP{$tmpProbe}=1;
                delete $pgfHOH{$tmpProbeSet}{probes}{$tmpProbe};
            }
            delete $pgfHOH{$tmpProbeSet};
            push (@maskedPS,$tmpProbeSet);
        }else{
            $pgfHOH{$tmpProbeSet}{start}=$min;
            $pgfHOH{$tmpProbeSet}{stop}=$max;
            $pgfHOH{$tmpProbeSet}{strand}=$tmpStrand;
            $pgfHOH{$tmpProbeSet}{chr}=$tmpChr;
        }
        if($procCount%10000==0){
            print ".";
        }
    }else{
        
    }
}
    print "\n\n";
    print "$pNonUnique Non-Unique removed\nRemove Probe sets <3 probes:\n";
    print $pRemoveCount." removed Probes\n";
    print $psRemoveCount." removed Probe sets\n";
    print $psRemChr." removed Chr Probe sets\n";
    print $psRemPloc." removed Loc Probe sets\n\n";
    @tmpList=keys %maskedP;
    
    
    print @tmpList." Total masked probes\n";
    print @maskedPS." Total masked probe sets\n\n";
    
    @tmpList=keys %pToPS;
    print @tmpList." Total Probes\n";
    @tmpList=keys %pgfHOH;
    print @tmpList." Total Probe Sets\n";
    
    
    
    my $tmpOutput=$arrayDir."/Masks/".$outputFileName;
    writePGF( \%pgfHOH, $header, $tmpOutput);
    
    foreach my $maskType(@maskList){
        if($maskType eq "all"){
            #read in probeset
            my $psHeader="";
            #print "reading:".$psPrefix.$maskType.".ps\n";
            open IN,"<",$psPrefix."core.ps";
            #read ps file
            #my $curDT=DateTime->now(time_zone=>'local');
            #my $dt_class='DateTime::Format::HTTP';
            my $line="#";
            while(index($line,"#")==0){
                $line=<IN>;
                if(index($line,"#")==0){
                    #if(index($line,"#%create_date=")==0){
                    #    $psHeader=$psHeader."#%create_date=".$dt_class->format_datetime($curDT)."\n";
                    #}else{
                        $psHeader=$psHeader.$line;
                    #}
                }
            }
            if(index($line,"probeset_id")==0){
                $psHeader=$psHeader.$line;   
            }
            my $tmpOut=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".ps";
            my $tmpLoc=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".psloc.bed";
            #print "writing:$tmpOut\n";
            open (OUTPS,">",$tmpOut) or die  "Cannot open $tmpOut: $!";
            open (OUTLOC,">",$tmpLoc) or die  "Cannot open $tmpLoc: $!";
            print OUTPS $psHeader;
            my @psList=keys %pgfHOH;
            foreach my $psID(@psList){
                print OUTPS $psID."\n";
                if(index($pgfHOH{$psID}{chr},"chr")>-1){
                    print OUTLOC $pgfHOH{$psID}{chr}."\t".$pgfHOH{$psID}{start}."\t".$pgfHOH{$psID}{stop}."\t".$psID."\t".$pgfHOH{$psID}{strand}."\n";
                }
            }
            close IN;
            close OUTPS;
            close OUTLOC;
        }else{
            my $psHeader="";
            #print "reading:".$psPrefix.$maskType.".ps\n";
            open IN,"<",$psPrefix.$maskType.".ps";
            #read ps file
            #my $curDT=DateTime->now(time_zone=>'local');
            #my $dt_class='DateTime::Format::HTTP';
            my $line="#";
            while(index($line,"#")==0){
                $line=<IN>;
                if(index($line,"#")==0){
                    #if(index($line,"#%create_date=")==0){
                    #    $psHeader=$psHeader."#%create_date=".$dt_class->format_datetime($curDT)."\n";
                    #}else{
                        $psHeader=$psHeader.$line;
                    #}
                }
            }
            if(index($line,"probeset_id")==0){
                $psHeader=$psHeader.$line;   
            }
            my $tmpOut=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".ps";
            my $tmpLoc=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".psloc.bed";
            #print "writing:$tmpOut\n";
            open (OUTPS,">",$tmpOut) or die  "Cannot open $tmpOut: $!";
            open (OUTLOC,">",$tmpLoc) or die  "Cannot open $tmpLoc: $!";
            print OUTPS $psHeader;
            while(<IN>){
                $line=$_;
                $line=trim($line);
                if(exists $pgfHOH{$line}){
                    print OUTPS $line."\n";
                    if(index($pgfHOH{$line}{chr},"chr")>-1){
                        print OUTLOC $pgfHOH{$line}{chr}."\t".$pgfHOH{$line}{start}."\t".$pgfHOH{$line}{stop}."\t".$line."\t".$pgfHOH{$line}{strand}."\n";
                    }
                }
            }
            close IN;
            close OUTPS;
            close OUTLOC;
        }
    }
    
    
    my @tcList=keys %uniqueTC;
    my $remIndivPSMask=0;
    my $remIndivPSLoc=0;
    my $remIndivPSOutlierLoc=0;
    my $remTCloc=0;
    my $remTClen=0;
    my $remTCps=0;
    my $remTCpsmasked=0;
    
    print "\n\nTC size:".@tcList."\n";
    
    print "Masking TC\n";
    #one loop should check for all at once
    
    foreach my $tcID(@tcList){
        my $tmpRef=$uniqueTC{$tcID}{tcprobes};
        my @psList=keys %$tmpRef;
        my $initPSCount=@psList;
        my %tmpChr;
        my %tmpStrand={'+'=>0,'-'=>0};
        my $finalStrand="";
        my $finalChr="";
        my $min=1999999999;
        my $max=0;
        my $tcRemoved=0;
        my $psToRemove=0;
        my @locations;
        my $tcLog="transcript cluster:".$tcID."\n";
        foreach my $psID(@psList){
            #check masked
            if(!(exists $pgfHOH{$psID})){
                $tcLog=$tcLog.$psID.": MASKED\n";
                delete $uniqueTC{$tcID}{tcprobes}{$psID};
                $remIndivPSMask++;
            }else{
                if(!($pgfHOH{$psID}{strand} eq "") and !($pgfHOH{$psID}{chr} eq "")){
                    $tcLog=$tcLog.$psID."\t".$pgfHOH{$psID}{cat}."\t".$pgfHOH{$psID}{strand}."\t".$pgfHOH{$psID}{chr}."\n";
                    if($pgfHOH{$psID}{start}<$min){
                        $min=$pgfHOH{$psID}{start};
                    }
                    if($pgfHOH{$psID}{stop}>$max){
                        $max=$pgfHOH{$psID}{stop};
                    }
                    push(@locations,$pgfHOH{$psID}{start});
                    
                    $tmpStrand{$pgfHOH{$psID}{strand}}=$tmpStrand{$pgfHOH{$psID}{strand}}+1;
                    
                    if(exists $tmpChr{$pgfHOH{$psID}{chr}}){
                        $tmpChr{$pgfHOH{$psID}{chr}}=$tmpChr{$pgfHOH{$psID}{chr}}+1;
                    }else{
                        $tmpChr{$pgfHOH{$psID}{chr}}=1;
                    }
                }else{
                    print $psID."non masked with no location\t".$pgfHOH{$psID}{cat}."\tst=".$pgfHOH{$psID}{strand}.">\tchr=".$pgfHOH{$psID}{chr}.">\n";
                }
            }
        }
        @psList=keys %$tmpRef;
        if(@psList>0){
            $tcLog=$tcLog."+".$tmpStrand{'+'}." vs -".$tmpStrand{'-'};
            #check location
            if($tmpStrand{'+'}>$tmpStrand{'-'} and $tmpStrand{'-'}>0){
                #print "+".$tmpStrand{'+'}." vs -".$tmpStrand{'-'}."\n";
                $finalStrand="+";
                $psToRemove=1;
            }elsif($tmpStrand{'+'}<$tmpStrand{'-'} and $tmpStrand{'+'}>0){
                #print "+".$tmpStrand{'+'}." vs -".$tmpStrand{'-'}."\n";
                $finalStrand="-";
                $psToRemove=1;
            }elsif($tmpStrand{'+'}==$tmpStrand{'-'}){
                #print "+".$tmpStrand{'+'}." vs -".$tmpStrand{'-'}."\n";
                $tcLog=$tcLog."== REMOVED";
                $tcRemoved=1;
                delete $uniqueTC{$tcID};
                $remTCloc++;
            }elsif($tmpStrand{'+'}>0){
                $finalStrand="+";
            }elsif($tmpStrand{'-'}>0){
                $finalStrand="-";
            }
            $tcLog=$tcLog."\nFinal:".$finalStrand."\t psremove:".$psToRemove."\n";
            if($tcRemoved==0){
                #check chr consensus
                #my @tmpChrList=keys %tmpChr;
                my @sortedKeys=sort {$tmpChr{$b} <=> $tmpChr{$a}} keys %tmpChr;
                $tcLog=$tcLog.@sortedKeys." Chromosomes\n";
                if(@sortedKeys>1){
                    $tcLog=$tcLog."Chr Count: ".$sortedKeys[0]."=".$tmpChr{$sortedKeys[0]}."\t".$sortedKeys[1]."=".$tmpChr{$sortedKeys[1]}."\n";
                    if($tmpChr{$sortedKeys[0]}==$tmpChr{$sortedKeys[1]}){
                        #print "Chr Count: ".$sortedKeys[0]."=".$tmpChr{$sortedKeys[0]}."\t".$sortedKeys[1]."=".$tmpChr{$sortedKeys[1]}."\n";
                        #mask TC
                        $tcRemoved=1;
                        delete $uniqueTC{$tcID};
                        $remTCloc++;
                    }else{
                        $finalChr=$sortedKeys[0];
                        #mask PS not on dominant Chr
                        $psToRemove=1;
                    }
                }else{
                    $tcLog=$tcLog."Chromosome:".$sortedKeys[0]."\n";
                    $finalChr=$sortedKeys[0];
                }
            }
            if($tcRemoved==0 and $psToRemove=1){#if probesets to remove and TC has not been removed
                #go through and remove probesets that don't match the final strand or chromosome
                $tcLog=$tcLog."Remove probesets:\n";
                foreach my $psID(@psList){
                    if($pgfHOH{$psID}{strand} eq $finalStrand and $pgfHOH{$psID}{chr} eq $finalChr){
                        #do nothing so far so good
                    }else{
                        #remove PS
                        $tcLog=$tcLog.$psID.",";
                        delete $uniqueTC{$tcID}{tcprobes}{$psID};
                        $remIndivPSLoc++;
                    }
                }
                $tcLog=$tcLog."\n";
            }
            if($tcRemoved==0){
                #check dist
                if(($max-$min)>1000000){
                    my @low;
                    my @high;
                    my @sortLoc=sort {$a<=>$b}@locations;
                    my $curLow=$sortLoc[0];
                    my $curHigh=$sortLoc[@sortLoc-1];
                    my $rem=0;
                    while(@sortLoc>1 and $curHigh-$curLow>1000000){
                        if($rem==0){
                            push(@low,shift(@sortLoc));
                            $rem=1;
                        }else{
                            push(@high,pop(@sortLoc));
                            $rem=0;
                        }
                        if(@sortLoc>1){
                            $curLow=$sortLoc[0];
                            $curHigh=$sortLoc[@sortLoc-1];
                        }
                    }
                    if(@sortLoc==1 and @locations>2){
                        $tcRemoved=1;
                        delete $uniqueTC{$tcID};
                        $remTClen++;
                    }else{
                        while($curHigh-$curLow<=1000000 and (($rem==1 and @high>0) or ($rem==0 and @low>0)) ){
                            if($rem==1){
                                push(@sortLoc,pop(@high));
                            }else{
                                unshift(@sortLoc,pop(@low));
                            }
                            $curLow=$sortLoc[0];
                            $curHigh=$sortLoc[@sortLoc-1];
                        }
                        if(@low>0 or @high>0){
                            my %toRemove;
                            foreach my $val(@low){
                                $toRemove{$val}=1;
                            }
                            foreach my $val(@high){
                                $toRemove{$val}=1;
                            }
                            @psList=keys %$tmpRef;
                            foreach my $psID(@psList){
                                if(exists $toRemove{$pgfHOH{$psID}{start}}){
                                    #remove PS
                                    delete $uniqueTC{$tcID}{tcprobes}{$psID};
                                    $remIndivPSOutlierLoc++;
                                }
                            }
                        }
                    }
                }
                
            }
            if($tcRemoved==0){
                #check # ps in tc
                $tmpRef=$uniqueTC{$tcID}{tcprobes};
                my @probeList=keys %$tmpRef;
                #my $percLeft=@probeList/$initPSCount*100;
                #if($percLeft<25){# if less than 25% of original PS left
                #    delete $uniqueTC{$tcID};
                #    $remTCps++;
                #}#els
                if(@probeList<1){# if no probeset left
                    delete $uniqueTC{$tcID};
                    $remTCps++;
                    
                    #print $tcLog."\n\n";
                    
                }
                else{# else set the location
                    $uniqueTC{$tcID}{start}=$min;
                    $uniqueTC{$tcID}{stop}=$max;
                    $uniqueTC{$tcID}{strand}=$finalStrand;
                    $uniqueTC{$tcID}{chr}=$finalChr;
                }
            }
        }else{
            $remTCpsmasked++;
            delete $uniqueTC{$tcID};
        }
    }
    
    print $remTClen." removed for length >300kbp\n";
    print $remTCloc." removed for location problems chr/strand\n";
    print $remTCps." removed for lack of probesets after QC\n";
    print $remTCpsmasked." removed all probesets were masked\n";
    
    @tcList=keys %uniqueTC;
    print "TC size:".@tcList."\n";
    
    print $remIndivPSMask." removed masked Probesets\n";
    print $remIndivPSLoc." removed for Chr/Strand\n";
    print $remIndivPSOutlierLoc." removed for outlier location\n";
    
    
    
    #WRITE MPS FILES
    foreach my $maskType(@maskList){
        my $maskedTC=0;
        if($maskType eq "all"){
            #read in probeset
            open (IN,"<",$psPrefix."full.mps") or die  "Cannot open ".$psPrefix."full.mps".": $!";
            #read ps file
            my $psHeader="";
            my $line="#";
            while(index($line,"#")==0){
                $line=<IN>;
                if(index($line,"#")==0){
                        $psHeader=$psHeader.$line;
                }
            }
            if(index($line,"probeset_id")==0){
                $psHeader=$psHeader.$line;   
            }
            
            my $tmpOut=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".mps";
            my $tmpLoc=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".transloc.bed";
            #print "writing:$tmpOut\n";
            open (OUTMPS,">",$tmpOut) or die  "Cannot open $tmpOut: $!";
            open (OUTLOC,">",$tmpLoc) or die  "Cannot open $tmpLoc: $!";
            print OUTMPS $psHeader;
            #my @tcList=keys %uniqueTC;
            print "output tc all:".@tcList."\n";
            foreach my $tcID(@tcList){
                
                    my $tmpRef=$uniqueTC{$tcID}{tcprobes};
                    my @tmpPSList=keys %$tmpRef;
                    my @sorted=sort @tmpPSList;
                    my $tcP=0;
                    my $tcPS=@tmpPSList;
                    my $psList="";
                    foreach my $tmpPS(@sorted){
                        if(exists $uniqueTC{$tcID}{tcprobes}{$tmpPS}){
                            my $tmpRef=$pgfHOH{$tmpPS}{probes};
                            my @tmpList=keys %$tmpRef;
                            $tcP=$tcP+@tmpList;
                            if($psList eq ""){
                                $psList=$tmpPS;
                            }else{
                                $psList=$psList." ".$tmpPS;
                            }
                        }else{
                            $tcPS--;
                        }
                    }
                    print OUTMPS $tcID."\t".$tcID."\t".$psList."\t".$tcP."\n";
                    if(index($uniqueTC{$tcID}{chr},"chr")>-1){
                        print OUTLOC $uniqueTC{$tcID}{chr}."\t".$uniqueTC{$tcID}{start}."\t".$uniqueTC{$tcID}{stop}."\t".$tcID."\t".$uniqueTC{$tcID}{strand}."\n";
                    }
            }
            while(<IN>){
                $line=$_;
                my @col=split(/\t/,$line);
                if($col[1] eq ""){ #check probesets since no transcript cluster
                    if(exists $pgfHOH{$col[0]}){
                        print OUTMPS $col[0]."\t\t".$col[2]."\n";
                        if(index($pgfHOH{$col[0]}{chr},"chr")>-1){
                            print OUTLOC $pgfHOH{$col[0]}{chr}."\t".$pgfHOH{$col[0]}{start}."\t".$pgfHOH{$col[0]}{stop}."\t".$col[0]."\t".$pgfHOH{$col[0]}{strand}."\n";
                        }
                    }
                }
            }
            close IN;
            close OUTMPS;
            close OUTLOC;
        }else{
            open (IN,"<",$psPrefix.$maskType.".mps") or die  "Cannot open ".$psPrefix.$maskType.".mps".": $!";
            #read ps file
            my $psHeader="";
            #$curDT=DateTime->now(time_zone=>'local');
            #$dt_class='DateTime::Format::HTTP';
            my $line="#";
            while(index($line,"#")==0){
                $line=<IN>;
                if(index($line,"#")==0){
                    #if(index($line,"#%create_date=")==0){
                    #    $psHeader=$psHeader."#%create_date=".$dt_class->format_datetime($curDT)."\n";
                    #}else{
                        $psHeader=$psHeader.$line;
                    #}
                }
            }
            if(index($line,"probeset_id")==0){
                $psHeader=$psHeader.$line;   
            }
            
            my $tmpOut=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".mps";
            my $tmpLoc=$arrayDir."/Masks/".$arrayStr.$ver.".".$maskType.".".$psSuffix.".transloc.bed";
            #print "writing:$tmpOut\n";
            open (OUTMPS,">",$tmpOut) or die  "Cannot open $tmpOut: $!";
            open (OUTLOC,">",$tmpLoc) or die  "Cannot open $tmpOut: $!";
            print OUTMPS $psHeader;
            my $tmpTotal=0;
            while(<IN>){
                $line=$_;
                my @col=split(/\t/,$line);
                if($col[1] eq ""){ #check probesets since no transcript cluster
                    if(exists $pgfHOH{$col[0]}){
                        print OUTMPS $col[0]."\t\t".$col[2]."\n";
                        print OUTLOC $pgfHOH{$col[0]}{chr}."\t".$pgfHOH{$col[0]}{start}."\t".$pgfHOH{$col[0]}{stop}."\t".$col[0]."\t".$pgfHOH{$col[0]}{strand}."\n";
                    }
                }else{
                    $tmpTotal++;
                    if(exists $uniqueTC{$col[1]}){#TC not Masked
                        my @tmpPSListMPS=split(/ /,$col[2]);
                        my @sorted=sort @tmpPSListMPS;
                        my $tcP=0;
                        my $tcPS=@tmpPSListMPS;
                        my $psList="";
                        foreach my $tmpPS(@sorted){
                            if(exists $uniqueTC{$col[1]}{tcprobes}{$tmpPS}){
                                my $tmpRef=$pgfHOH{$tmpPS}{probes};
                                my @tmpList=keys %$tmpRef;
                                $tcP=$tcP+@tmpList;
                                if($psList eq ""){
                                    $psList=$tmpPS;
                                }else{
                                    $psList=$psList." ".$tmpPS;
                                }
                            }else{
                                $tcPS--;
                            }
                        }
                        if($tcPS>0){
                            print OUTMPS $col[0]."\t".$col[1]."\t".$psList."\t".$tcP."\n";
                            if(index($uniqueTC{$col[1]}{chr},"chr")>-1){
                                print OUTLOC $uniqueTC{$col[1]}{chr}."\t".$uniqueTC{$col[1]}{start}."\t".$uniqueTC{$col[1]}{stop}."\t".$col[1]."\t".$uniqueTC{$col[1]}{strand}."\n";
                            }
                        }
                    }else{
                        $maskedTC++;
                    }
                }
            }
            print $maskedTC." ".$maskType." masked out of ".$tmpTotal."\n";
            close IN;
            close OUTMPS;
            close OUTLOC;
        }
    }
    