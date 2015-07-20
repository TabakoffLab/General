

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


my $pgf="MoEx-1_0-st-v1.r2.LXS.MASKED.perl.pgf";
my $arrayPath="/Volumes/Data/mm10/array/MoEx/Masks/LXS Mask/";

print "try to open: $pgf\n";

my %pToPS;

    
open(IN, "<".$arrayPath.$pgf);
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


my @keyList=keys %pgfHOH;
open(OUT, ">/Users/smahaffey/desktop/LXS.probe_count.txt");
foreach $key(@keyList){
    my $tmp=$pgfHOH{$key}{probes};
    my %tmpHOH=%$tmp;
    my @psList=keys %tmpHOH;
    print OUT $key."\t".@psList."\n";
}
close OUT;