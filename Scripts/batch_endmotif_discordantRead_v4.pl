use strict;
my $list=shift;
my %K6_motif=();
my %K4_motif=();
my $Ends_equal_k4=0;
my $N=0;
my $Ends_equal_k6=0;
my $out="equal_twoends_".$list;
open(LI,$list)or die $!;
open(OUT,">>$out")or die $!;
	print OUT "file\tTotall moleculars\t4-mer_equal proportion\t6-mer_equal proportion\n";
while(my $file=<LI>){
	chomp $file;
    $Ends_equal_k4=0;
    $N=0;
    $Ends_equal_k6=0;
    my %read=();

     open(IN,$file)or die $!;
	while(my $line=<IN>){
		if($line=~/^\@/){
		# print $line;		
		}else{
			my @temp=split /\t/,$line;
			if(!exists($read{$temp[0]})){
				$read{$temp[0]}=$line;
			}else{
				$read{$temp[0]}.=$line;
			}
		}
	}
	close(IN);

	my @read_id=keys(%read);
	
	for(my $i=0;$i<@read_id;$i++){
   	 	my @temp_k6=();
   	 	my @temp_k4=();
   	 	my $k6="";
    	my $k4="";
		my @pair=split /\n/, $read{$read_id[$i]};
		for my $each(@pair){
			my @temp2=split /\t/,$each;
		
		my $test=$temp2[1] & 16;
		if($test == 0){
			$k6=substr($temp2[9],0,6);
     	   $k4=substr($temp2[9],0,4);
		}elsif($test == 16){
			$k6=reverse(substr($temp2[9],-6));
        	my @tempkmer=split //,$k6;
       	 	my $tempk6="";
        	for(@tempkmer){
        	    if($_ eq "A"){
        	    	$tempk6.="T"
         	   }elsif($_ eq "T"){
          		  	$tempk6.="A"
            	}elsif($_ eq "G"){
            		$tempk6.="C"
            	}elsif($_ eq "C"){
            		$tempk6.="G"
            	}else{  # for N
					$tempk6.=$_;
            			print "error...$k6....$tempk6....$temp2[9]...\n";
            	}
        	}
        	$k6=$tempk6;
            	
        	$k4=reverse(substr($temp2[9],-4));
        	my @tempkmer=split //,$k4;
        	my $tempk4="";
        	for(@tempkmer){
            	if($_ eq "A"){
            		$tempk4.="T"
            	}elsif($_ eq "T"){
            		$tempk4.="A"
            	}elsif($_ eq "G"){
            		$tempk4.="C"
            	}elsif($_ eq "C"){
            		$tempk4.="G"
            	}else{  # for N 
            		$tempk4.=$_;
           	 	}
       		 }
            	
        	$k4=$tempk4;
		}else{
					print "error1...\n";
		}
		
   	 	
     	push @temp_k4,$k4;
     	push @temp_k6,$k6;
     } 
     if(@temp_k4==2){
     	my $k3_0_0=substr($temp_k4[0],0,3);
     	my $k3_0_1=substr($temp_k4[0],-3);
     	my $k3_1_0=substr($temp_k4[1],0,3);
     	my $k3_1_1=substr($temp_k4[1],-3);
     	if($temp_k4[0] eq $temp_k4[1] or $k3_0_0 eq $k3_1_1 or $k3_0_1 eq $k3_1_0){
     		$Ends_equal_k4++;
     		if(exists($K4_motif{$k4})){
        		$K4_motif{$k4}++;
               # print $k4."\t".$K4_motif{$k4}."\n";
     		}else{
        		$K4_motif{$k4}=1; 
     		} 
     	}
     	$N++;
     }
     if(@temp_k6==2){
     	my $k5_0_0=substr($temp_k6[0],0,5);
     	my $k5_0_1=substr($temp_k6[0],-5);
     	my $k5_1_0=substr($temp_k6[1],0,5);
     	my $k5_1_1=substr($temp_k6[1],-5);
     	if($temp_k6[0] eq $temp_k6[1] or $k5_0_0 eq $k5_1_1 or $k5_0_1 eq $k5_1_0){
     		$Ends_equal_k6++;
     		if(exists($K6_motif{$k6})){
        		$K6_motif{$k6}++;
    		}else{
       		 	$K6_motif{$k6}=1; 
   			}
     	}
     	
     }
     
       
	}
	open(OUT,">>$out")or die $!;
	print OUT "$file\t$N\t".($Ends_equal_k4/$N)."\t". ($Ends_equal_k6/$N)."\n";
	
	
}

my $kmer="";
my @k6mer=sort{$a<=>$b}keys(%K6_motif);
foreach(@k6mer){   
   	$kmer.=$_."\t".$K6_motif{$_}."\n"; 
}
my $out="k6mer_equal_".$list;
open(OUTT,">$out")or die $!;
print OUTT $kmer;
close(OUTT);
	
my $kmer="";
my @k4mer=sort{$a<=>$b}keys(%K4_motif);
foreach(@k4mer){   
   	$kmer.=$_."\t".$K4_motif{$_}."\n"; 
}
my $out="k4mer_equal_".$list;
open(OUTT,">$out")or die $!;
print OUTT $kmer;
close(OUTT);
close(OUT);