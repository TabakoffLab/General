my $perfectmatch=$ARGV[2]; #file for perfect matches
my $probeFA=$ARGV[3];   #file for all probes


my $output=$ARGV[6]; #path file name for output file

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
    }elsif($type eq "AffyEx"){
        
    }else{
        
    }