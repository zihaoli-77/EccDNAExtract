=cut
	version 6: extract the identical 6-mer motifs in Type A fragments, and count the number of each 6-mer motifs. 
=cut
use strict;
my $list=shift;
my $ReadLenth=75;

my %split_count=();
my %read=();
my %lenth_dis=();
my %lenth_ecc=();
my $N0=0;
my $N1=0;
my $N2=0;
my $N3=0;
my $N4=0;
my $N5=0;
my $N6=0;
my $N7=0;
my $N8=0;
my $N9=0;
my %dis_total=();
my %dis_op=();
my %dis_same=();
my %K6_motif=();
my %K4_motif=();
my %K6_motif_split=();
my %K4_motif_split=();

my $Ends_equal_k4=0;
my $Ends_equal_k4_same=0;
my $Ends_equal_k4_op=0;
my $Ends_equal_k6=0;
my $Ends_equal_k6_same=0;
my $Ends_equal_k6_op=0;
my $Nequal=0;
my $header="";
my $N7_reads="";

my $statics_f="statics_".$list;
open(OUT3,">$statics_f")or die $!;
	print OUT3 "Filename\ttotal_frag\tcan_not_estimate\tcan_estimate\testimate_frag_length\testimate_ecc_len\testimate_both\tTwoEnds_samechr\tTwoEnds_difOrigin\tTwoEnds_sameOrigin\testimatable_TwoEnds_difOrigin\testimatable_TwoEnds_sameOrigin\tTwoEnds_same4mer_difOrigin\tTwoEnds_same4mer_sameOrigin\tTwoEnds_same6mer_difOrigin\tTwoEnds_same4mer_sameOrigin\n";

open(LI,$list)or die $!;
while(my $file=<LI>){
	chomp $file;
	%split_count=();
	%read=();
	#	%lenth_dis=();
	$N0=0;
	$N1=0;
	$N2=0;
	$N3=0;
	$N4=0;
	$N5=0;
	$N6=0;
	$N7=0;
	$Ends_equal_k4=0;
    $Ends_equal_k4_same=0;
    $Ends_equal_k4_op=0;
    $Ends_equal_k6=0;
    $Ends_equal_k6_same=0;
    $Ends_equal_k6_op=0;
    $Nequal=0;
    $header="";
    
    $N7_reads="";
	# %dis_total=();
	# %dis_op=();
	#	%dis_same=();
	print "Processing....$file\n";
	open(IN,$file)or die $!;
	while(my $line=<IN>){
		if($line=~/^\@/){
		# print $line;	
			$header.=$line;	
		}else{
			my @temp=split /\t/,$line;
			if(!exists($read{$temp[0]})){
				$split_count{$temp[0]}=1;
				$read{$temp[0]}=$line;
			}else{
				$split_count{$temp[0]}++;
			#	print $read{$temp[0]};
				$read{$temp[0]}.=$line;
			}
		}
	}
	close(IN);
	#my $outfile="cluster_".$file;
	#open(OUT,">$outfile")or die $!;
	#	my $process_inf="details_".$file;
	#open(OUT2,">$process_inf")or die $!;
	my $pos5end_f="pos5end_".$file;
	open(OUT1,">$pos5end_f")or die $!;
	my $N7reads_f="N7reads_".$file;
	open(OUT7,">$N7reads_f")or die $!;
    my $cluster_string="";
    my $detail_string="";
    my $end_string="";
    my @read_id=keys(%read);
	my %statics=();
	for(my $i=0;$i<@read_id;$i++){
		if(exists($statics{$split_count{$read_id[$i]}})){
	    	  $statics{$split_count{$read_id[$i]}}++;
		}else{
        	 $statics{$split_count{$read_id[$i]}}=1;
		}
		$cluster_string.=$read{$read_id[$i]};
    	my ($details,$end_inf)=&call_5endto3end($read{$read_id[$i]});
    	$detail_string.=$details;
		$end_string.=$end_inf;
	
	}
	#print OUT $cluster_string;
	print OUT1 $end_string;
	close(OUT1);
	print OUT7  $header.$N7_reads;
	close(OUT7);
	#print OUT2 $detail_string;

	#	print OUT3 "read_fragments_count...count\n";
	#	my @count=sort{$a<=>$b}keys(%statics);
	#foreach(@count){
	#	print OUT3 $_."\t".$statics{$_}."\n";
	#}
	my $output_string=$file;
	$output_string.="\t".($N0+$N1)."\t".($N0/($N0+$N1))."\t".($N1/($N0+$N1))."\t".($N2/($N0+$N1))."\t".($N3/($N0+$N1))."\t".($N4/($N0+$N1))."\t".($N5/($N0+$N1))."\t".($N6/($N0+$N1))."\t".($N7/($N0+$N1))."\t".($N8/($N0+$N1))."\t".($N9/($N0+$N1))."\t".($Ends_equal_k4_same/$N6)."\t".($Ends_equal_k4_op/$N7)."\t".($Ends_equal_k6_same/$N6)."\t".($Ends_equal_k6_op/$N7)."\n";
	print OUT3 $output_string;
	#	print OUT3 "$N0 >150bp:\t".($N0/($N0+$N1))."\ttatolly different fragments origins, out of estimated\n";
	#	print OUT3 "$N1\t:".($N1/($N0+$N1))."\thaving at least one fragment of a same origin, can be partially estimated\n";

}

=cut
	
	my $len_ecc_string="";
	my @len=sort{$a<=>$b}keys(%lenth_ecc);
	foreach(@len){
		$len_ecc_string.=$_."\t".$lenth_ecc{$_}."\n";
	}
	my $out="len_ecc_dis_".$list;
	open(OUTC,">$out")or die $!;
	print OUTC $len_ecc_string;
	close(OUTC);
	

	my $len_string="";
	my @len=sort{$a<=>$b}keys(%lenth_dis);
	foreach(@len){
		$len_string.=$_."\t".$lenth_dis{$_}."\n";
	}
	my $out="len_dis_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $len_string;
	close(OUTT);


	my $end_dis_string="";
	my @end_dis=sort{$a<=>$b}keys(%dis_total);
	foreach(@end_dis){ 
		$end_dis_string.=$_."\t".$dis_total{$_}."\n"; 
	}
	#	print OUT3 $N2."\t:".($N2/($N0+$N1))."\tThe distance between two ends mapped to the same chromsome\n";
	my $out="end_dis_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $end_dis_string;
	close(OUTT);

	my $end_dis_sd="";
	my @end_dis=sort{$a<=>$b}keys(%dis_same);
	foreach(@end_dis){   
    	$end_dis_sd.=$_."\t".$dis_same{$_}."\n"; 
	}
	#	print OUT3 $N3."\t:".($N3/($N0+$N1))."\tThe distance between two ends of same directions mapped to the same chromsome\n";
	my $out="end_dis_sd_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $end_dis_sd;
	close(OUTT);

	
	my $end_dis_op="";
	my @end_dis=sort{$a<=>$b}keys(%dis_op);
	foreach(@end_dis){   
   		$end_dis_op.=$_."\t".$dis_op{$_}."\n"; 
	}
	#	print OUT3 $N4."\t:".($N4/($N0+$N1))."\tThe distance between two ends of opposite directions mapped to the same chromsome\n";
	my $out="end_dis_op_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $end_dis_op;
	close(OUTT);
	
=cut	

	my $kmer="";
	my @k6mer=sort{$a<=>$b}keys(%K6_motif);
	foreach(@k6mer){   
   		$kmer.=$_."\t".$K6_motif{$_}."\n"; 
	}
	my $out="k6mer_TypeA_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);
	
	my $kmer="";
	my @k4mer=sort{$a<=>$b}keys(%K4_motif);
	foreach(@k4mer){   
   		$kmer.=$_."\t".$K4_motif{$_}."\n"; 
   	
	}
	my $out="k4mer_TypeA_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);
	
	my $kmer="";
	my @k6mer=sort{$a<=>$b}keys(%K6_motif_split);
	foreach(@k6mer){   
   		$kmer.=$_."\t".$K6_motif_split{$_}."\n"; 
	}
	my $out="k6mer_split_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);
	
	my $kmer="";
	my @k4mer=sort{$a<=>$b}keys(%K4_motif_split);
	foreach(@k4mer){   
   		$kmer.=$_."\t".$K4_motif_split{$_}."\n"; 
   	
	}
	my $out="k4mer_split_equal_".$list;
	open(OUTT,">$out")or die $!;
	print OUTT $kmer;
	close(OUTT);
	
	
	

close();

sub call_5endto3end{
	my $pair=shift;
	#print $pair;
	my $end5="";
	my %start=();
    my @temp=split /\n/, $pair;
	for(@temp){
		if($_ ne ""){
			my @temp1=split /\t/,$_;
			my $read_flag=$temp1[1] & 192;
#			print "read_flag_".$read_flag."\n";
			if(!exists($start{$read_flag})){
				$start{$read_flag}=$_;
			}else{
				$start{$read_flag}.="\n".$_;
				#	print "Here correct?\n";
				#print $start{$read_flag}."\n";
			}
		}
	}
	my $read_order=0;
	my $start_inf="";
	my $read_id="";
	my @read=keys(%start);
	for(@read){
		#	print $start{$_}."\n";
	     $read_order++;
		my @temp_fragment=split /\n/,$start{$_};
		#print "@temp_fragment\n";
		my $length=0; #length mapped to reference
		my $length2=0; #length of segment
		if(@temp_fragment==1){
		    my @temp2=split /\t/,$temp_fragment[0];
		    
		    if(($temp2[1] & 4) ==4){
		    	$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t*\t*\t*\t$ReadLenth\t*\t$temp2[9];\n";
		    	$end5.="\tR$read_order\t$temp2[2]\t*\t*\t$temp2[9]";
		    	next;
		    }
			
			my @CIGAR=split /[A-Z|=]/,$temp2[5];
			my @CIGAR2=split /[0-9]+/,$temp2[5];
			if($CIGAR2[0] eq ""){
				shift @CIGAR2;
			}	
			#	print "@CIGAR...@CIGAR2\n";
			#  print @CIGAR."...".@CIGAR2."\n";
			$read_id=$temp2[0];
			if(@CIGAR == @CIGAR2){
				my $start_base=0;
				if((($temp2[1] & 16) ==0) && $CIGAR2[0] ne "M"){
					$start_base=$CIGAR[0];
				}
				if((($temp2[1] & 16) ==16) && $CIGAR2[-1] ne "M"){
					$start_base=$CIGAR[-1];
				}

				for (my $i=0;$i<@CIGAR2;$i++){
					if($CIGAR2[$i] eq "M"){
						$length+=$CIGAR[$i];
						$length2+=$CIGAR[$i];
					}
					if($CIGAR2[$i] eq "I"){
						$length2+=$CIGAR[$i];
					}
					if($CIGAR2[$i] eq "D"){
						$length+=$CIGAR[$i];
					}
				}
				#	$end5.="\t".$start_base;
				my $test=$temp2[1] & 16;
				my $binary_number = sprintf("%b", $temp2[1]);
				#				print "test...$binary_number\t".$test."\n";
				if($test == 0){
			    	$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".$temp2[3]."\t".($temp2[3]+$length-1)."\t".$start_base."\t".($length2)."\t+\t$temp2[9];\n";
					if($start_base==0){
						$end5.="\tR$read_order\t$temp2[2]\t$temp2[3]\t+\t$temp2[9]";
					}else{
						$end5.="\tR$read_order\t$temp2[2]\t-\t+\t$temp2[9]";               
					}
					#	print "one frag....+\n";
					#   print $start_inf;
		   		}elsif($test == 16){
			   		$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".($temp2[3]+$length-1)."\t".$temp2[3]."\t".$start_base."\t".($length2)."\t-\t$temp2[9];\n";
					if($start_base==0){
						$end5.="\tR$read_order\t$temp2[2]\t".($temp2[3]+$length)."\t-\t$temp2[9]";
					}else{
						$end5.="\tR$read_order\t$temp2[2]\t-\t-\t$temp2[9]";
					}
				}else{
					print "error1...\n";
				}
				#	print "inf1....$start_inf";
			}else{
					print "error...@CIGAR...@CIGAR2...$temp2[5]\n";
			}
			
		}elsif(@temp_fragment >= 2){
			my %start_frag=();
			my %start_base=();
			for my $each_frag(@temp_fragment){
			    my @temp2=split /\t/,$each_frag;
				if(($temp2[1] & 4) ==4){
					 $length=length($temp2[9]);
		    		$start_inf.=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t*\t*\t*\t$length\t*\t$temp2[9];";
		    		next;
		        }
			    $length=0;
			    $length2=0;
				
                my @CIGAR2=split /[0-9]+/,$temp2[5];
			    my @CIGAR=split /[A-Z|=]/,$temp2[5];
		     	if($CIGAR2[0] eq ""){
			    	shift @CIGAR2;
			    }	
				my $strand="";
				$read_id=$temp2[0];
			   if(@CIGAR == @CIGAR2){
					$start_base{$each_frag}=0;
					my $i=0;
					my $j=0;
					
		   			if(($temp2[1] & 16) ==0){
						$i=0;
						$j=1;
						$strand="+";
               		 }elsif(($temp2[1] & 16) == 16){
						$i=-1;
						$j=-1;
						$strand="-";
					}
					my $test=$temp2[1] & 16;
                    my $binary_number = sprintf("%b", $temp2[1]);
					#   print "test...$binary_number\t".$test."\n";
					#	print "hello1...$i...$j\n";
					while (($i<0 && abs($i)<=@CIGAR2)||($i>=0 && $i<@CIGAR2)){
                        if($CIGAR2[$i] eq "M"){
                           $length+=$CIGAR[$i];
                           $length2+=$CIGAR[$i];
						   $i+=$j;
						   last;
                        }else{
	                       $start_base{$each_frag}+=$CIGAR[$i];
						}
						$i+=$j;
						#	print "hello2...$i...$j\n";
					}
					for ($i;($i<0 && abs($i)<=@CIGAR2)||($i>=0 && $i<@CIGAR2);$i+=$j){
						if($CIGAR2[$i] eq "M"){
							$length+=$CIGAR[$i];
							$length2+=$CIGAR[$i];
						}elsif($CIGAR2[$i] eq "I"){
							$length2+=$CIGAR[$i];
                    	}elsif($CIGAR2[$i] eq "D"){
                        	$length+=$CIGAR[$i];
                    	}
						#	print "hello3...$i...$j\n";
					}
               
					if(($temp2[1] & 16) ==0){
						$start_frag{$each_frag}=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".$temp2[3]."\t".($temp2[3]+$length-1)."\t".$start_base{$each_frag}."\t".($length2)."\t+\t$temp2[9];";
						#	print "two more frag....+\n";
						#	print $start_frag{$each_frag}."\n";
					}elsif(($temp2[1] & 16) == 16){
						$start_frag{$each_frag}=$temp2[0]."\t"."R".$read_order."\t".$temp2[2]."\t".($temp2[3]+$length-1)."\t".$temp2[3]."\t".$start_base{$each_frag}."\t".($length2)."\t-\t$temp2[9];";
					}
					#print "inf2....$start_inf";
               }else{
                        print "error...@CIGAR...@CIGAR2...$temp2[5]\n";
               }
		   }
		   my @sort_frag=sort{$start_base{$a}<=>$start_base{$b}}keys(%start_base);
		   my $k=0;
		   for my $each_frag(@sort_frag){
			   if($k==0){
			   	   my @temp3=split /\t/,$start_frag{$each_frag};
				   my $strand="";
				   if($temp3[-2]=~/(.+)/){
				   		$strand=$1;
			   		}else{
						$strand=$temp3[-2];
					}
					#	print $strand."\n";
					$temp3[-1]=~/(.+);/;
				   if($temp3[-4]==0){
				   		$end5.="\tR".$read_order."\t".$temp3[2]."\t".$temp3[3]."\t".$strand."\t".$1;
			   		}else{
						$end5.="\tR".$read_order."\t".$temp3[2]."\t-\t".$strand."\t".$1;
					}
			   }
			   $k++;
			   $start_inf.=$start_frag{$each_frag};
	  		 }
			$start_inf.="\n";
			#		print "total pair+frag...\n";
			#	print $start_inf;
		}
	}
	$end5=$read_id.$end5;
	#print $end5."\n";
	#	print "total pair+frag...\n";
	my($flag,$len,$con,$cnv)=&estimate_length($start_inf);
	#    print $flag."\t".$len."\t".$con."\n";
	$end5.="\t".$flag."\t".$len."\t".$con."\t".$cnv."\n";
	#print $end5;
	$start_inf=$start_inf.$end5;
	my @temp_end=split /\t/,$end5;
	 my $flag_equal_k4=0;
	 my $flag_equal_k6=0;
	 my @temp_k6=();
   	 my @temp_k4=();
   	 my $k6="";
    my $k4="";
	if(($temp_end[3] ne "*" && $temp_end[8] ne "*") && ($temp_end[3] ne "-") && ($temp_end[8] ne "-")){
	   
       ####kmer of last read
            if($temp_end[9] eq "+"){
            	$k6=substr($temp_end[10],0,6);
            	$k4=substr($temp_end[10],0,4);
            	
            }else{
            	$k6=reverse(substr($temp_end[10],-6));
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
						$tempk6+=$_;
            			#print "error...$k6....$tempk6....$temp_end[5]...\n";
    
            	 
            		}
            	}
            	$k6=$tempk6;
            	
            	$k4=reverse(substr($temp_end[10],-4));
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
            	#if($k4=~/N/){
            	#		print "error...$k4....$tempk4....$temp_end[5]...\n";
    
            	#}
            	$k4=$tempk4;
            	
            }
            push @temp_k4,$k4;
     	    push @temp_k6,$k6;
     	    
     	  ####kmer of former reads
     	  my $k6="";
            my $k4="";
            if($temp_end[4] eq "+"){
            	$k6=substr($temp_end[5],0,6);
            	$k4=substr($temp_end[5],0,4);
            	
            }else{
            	$k6=reverse(substr($temp_end[5],-6));
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
						$tempk6+=$_;
            			#print "error...$k6....$tempk6....$temp_end[5]...\n";
    
            	 
            		}
            	}
            	$k6=$tempk6;
            	
            	$k4=reverse(substr($temp_end[5],-4));
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
            	#if($k4=~/N/){
            	#		print "error...$k4....$tempk4....$temp_end[5]...\n";
    
            	#}
            	$k4=$tempk4;
            	
            }
   	 	    push @temp_k4,$k4;
     	    push @temp_k6,$k6;
     	    my $flag_equal_k4=0;
     	    if(@temp_k4==2){
     	       my $k3_0_0=substr($temp_k4[0],0,3);
     	       my $k3_0_1=substr($temp_k4[0],-3);
     	       my $k3_1_0=substr($temp_k4[1],0,3);
     	       my $k3_1_1=substr($temp_k4[1],-3);
     	       if($temp_k4[0] eq $temp_k4[1] or $k3_0_0 eq $k3_1_1 or $k3_0_1 eq $k3_1_0){
				 $flag_equal_k4=1;
				 for my $k4(@temp_k4){
            	 	if(exists($K4_motif_split{$k4})){
                    	$K4_motif_split{$k4}++;
              			 # print $k4."\t".$K4_motif{$k4}."\n";
              	 	}else{
                	  	$K4_motif_split{$k4}=1; 
               	 	}
           		 }
     	      }
     	      
             }
		  
           if(@temp_k6==2){
     	    my $k5_0_0=substr($temp_k6[0],0,5);
     	    my $k5_0_1=substr($temp_k6[0],-5);
     	    my $k5_1_0=substr($temp_k6[1],0,5);
     	    my $k5_1_1=substr($temp_k6[1],-5);
     	   
     	    if($temp_k6[0] eq $temp_k6[1] or $k5_0_0 eq $k5_1_1 or $k5_0_1 eq $k5_1_0){
			   $flag_equal_k6=1; 
			    for my $k6(@temp_k6){
            	 if(exists($K6_motif_split{$k6})){
                     $K6_motif_split{$k6}++;
                 }else{
                     $K6_motif_split{$k6}=1; 
                 }
           		 }  
     	    }
     	
           }
	
	
	
	}
	
	
	if(($temp_end[2] eq $temp_end[7]) && ($temp_end[3] ne "*" && $temp_end[8] ne "*") && ($temp_end[3] ne "-") && ($temp_end[8] ne "-")){
		my $dis=abs($temp_end[3]-$temp_end[8]);
		$N5++;
		if(exists($dis_total{$dis})){
			$dis_total{$dis}++;
		}else{
			$dis_total{$dis}=1;
		}	
		
		
		if( $flag_equal_k4==1) {
			$Ends_equal_k4++;
				
     		if($temp_end[4] eq $temp_end[9]){
     		     $Ends_equal_k4_same++;
     		}else{
     		     $Ends_equal_k4_op++;
     		}
     		$Nequal++;
		}
		if( $flag_equal_k6==1) {
			$Ends_equal_k6++;
			   
     		if($temp_end[4] eq $temp_end[9]){
     		     $Ends_equal_k6_same++;
     		 }else{
     		     $Ends_equal_k6_op++;
     		 }
		}
		########same orientation from different molecular
		if($temp_end[4] eq $temp_end[9]){
			$N6++;
			if(exists($dis_same{$dis})){
                $dis_same{$dis}++;
            }else{
                $dis_same{$dis}=1; 
            }
			if($flag_equal_k6==1){
            for my $k6(@temp_k6){
            	 if(exists($K6_motif{$k6})){
                     $K6_motif{$k6}++;
                 }else{
                     $K6_motif{$k6}=1; 
                 }
            }
		    }
            if($flag_equal_k4==1){
            for my $k4(@temp_k4){
            	if(exists($K4_motif{$k4})){
                    $K4_motif{$k4}++;
               # print $k4."\t".$K4_motif{$k4}."\n";
                }else{
                    $K4_motif{$k4}=1; 
                }
            }
			}
          
            
            if($flag==1 or $flag==3){
            	$N8++;
            }
		}else{
			$N7++;
			$N7_reads.=$pair;
			if(exists($dis_op{$dis})){
                $dis_op{$dis}++;
            }else{
                $dis_op{$dis}=1;
		    }
		    if($flag == 2){
            	$N9++;
            }

		}
	}
	if($flag > 0){
		$N1++;
		if($flag==1 or $flag ==3){
			$N2++;
			if(!exists($lenth_dis{$len})){
				$lenth_dis{$len}=1;
			}else{
				$lenth_dis{$len}++;
			}
		}
		if($flag ==3){
			$N4++;
		}
		if($flag ==2 or $flag==3){
			$N3++;
			if(!exists($lenth_ecc{$len})){
				$lenth_ecc{$len}=1;
			}else{
				$lenth_ecc{$len}++;
			}
		}
	}else{
		$N0++;
	}
	#	print $start_inf;
	return ($start_inf,$end5);
}



sub estimate_length{
	#fragments has different aligned chromosome, return 150+
	#two or more fragments have the same chromosome, try to estimate length
	my $read_pair=shift;
	my $con_inf="";
	my $cnv_inf="";
	my %flag_R1=();
	my %flag_R2=();
	#print $read_pair;
	my @temp_read=split /\n/,$read_pair;
	my @read_1=split /;/,$temp_read[0];
	my @read_2=split /;/,$temp_read[1];
	if ($read_1[-1] eq ""){
		pop @read_1;
	}
	if ($read_2[-1] eq ""){
        pop @read_2;
     } 
    my %chr=();
	#case 1 read_1_tail_vs_read_2_head
	my @r1t=split /\t/,$read_1[-1];
	my @r2h=split /\t/,$read_2[-1];
	my $i=@read_1-1;
	my $j=@read_2-1;
 	my	$length=0;
	if($r1t[2] eq $r2h[2] && $r2h[3] ne "*" && $r1t[3] ne "*"){
		if($r1t[-2] eq "+"){#$r1t[3]-->$r1t[4]-->$r2h[4]-->$r2h[3]
			if($r2h[3]>$r2h[4] && $r2h[3]>=$r1t[3]){
				my $min=($r2h[4]<$r1t[3])?$r2h[4]:$r1t[3];
				$length+=$r2h[3]-$min+1;	
				$cnv_inf.="$r1t[2]:$min-$r2h[3];";
				$flag_R1{$i}=1;
				$flag_R2{$j}=1;
                #print "C1:R1$i...R2$j...$length\n";
				$con_inf.="R1$i<..<R2$j;";
				#print $cnv_inf."\n";
			}

		}else{#
			if($r2h[3]<$r2h[4] && $r2h[3]<=$r1t[3]){
				my $max=($r2h[4]<$r1t[3])?$r1t[3]:$r2h[4];
				$length+=abs($max-$r2h[3]+1);
				$cnv_inf.="$r1t[2]:$r2h[3]-$max;";
				$flag_R1{$i}=1;
				$flag_R2{$j}=1;
				#	print "C2:R1$i...R2$j...$length\n";
				$con_inf.="R1$i>..>R2$j;";
				#				 print $cnv_inf."\n";
			}
		}
	}
	
	#case 2 read_1_head_vs_read_2_tail
	my @r1h=split /\t/,$read_1[0];
	my @r2t=split /\t/,$read_2[0];
	my $i=0;
    my $j=0;
	if($r1h[2] eq $r2t[2] && $r2h[3] ne "*" && $r1t[3] ne "*"){
		if($r1h[-2] eq "+"){ #$r1h[4]-->$r1h[3]-->$r2t[3]-->$r2t[4]
			if($r2t[3]>$r2t[4] && $r2t[3]<=$r1h[3]){
				if(exists($flag_R1{$i})){
					$length+=$r1h[3]-$r2t[4]+1;
					$cnv_inf.="$r1h[2]:$r2h[4]-$r1h[3];";
		
				}elsif(exists($flag_R2{$j})){
					 $length+=$r1h[4]-$r2t[3]+1;
				     $cnv_inf.="$r1h[2]:$r2t[3]-$r1h[4];";
				 }else{
					$length+=$r1h[4]-$r2t[4]+1;
					$cnv_inf.="$r1h[2]:$r2h[4]-$r1t[4];";
				}
				if(exists($flag_R1{$i})){
					$flag_R1{$i}=3;
				}else{
					$flag_R1{$i}=2;
				}
				if(exists($flag_R2{$i})){
					$flag_R2{$j}=3;
				}else{
					$flag_R2{$j}=2;
				}
				#	print "C3:R1$i...R2$j...$length\n";
				# print $cnv_inf."\n";
				$con_inf.="R2$j<..<R1$i;";
			}

		}else{#$r2t[4]-->$r2t[3]-->$r1h[3]-->$r1h[4]
			if($r2t[3]<$r2t[4] && $r2t[3]>=$r1h[3]){
				if(exists($flag_R1{$i})){
					 $length+=abs($r2t[4]-$r1h[3]+1);
				     $cnv_inf.="$r1h[2]:$r1h[3]-$r2t[4];";
		
				}elsif(exists($flag_R2{$j})){
					 $length+=abs($r2t[3]-$r1h[4]+1);
				     $cnv_inf.="$r1h[2]:$r1h[4]-$r2t[3];";
				}else{
					 $length+=abs($r2t[4]-$r1h[4]+1);
				     $cnv_inf.="$r1h[2]:$r1h[4]-$r2t[4];";
				}
				
				if(exists($flag_R1{$i})){
					$flag_R1{$i}=3;
				}else{
					$flag_R1{$i}=2;
				}
				if(exists($flag_R2{$i})){
					$flag_R2{$j}=3;
				}else{
					$flag_R2{$j}=2;
				}
				# print "C4:R1$i...R2$j...$length\n";
				#print $cnv_inf."\n";
				$con_inf.="R2$j>..>R1$i;";
			}
		}
	}
	my $read_length1=0;
	my $flag=0;
   for (my $i=0;$i<@read_1;$i++){
			my @r1=split /\t/,$read_1[$i];
			$read_length1+=$r1[-3];
			$length+=$r1[-4];
		if(!exists($flag_R1{$i})){
			$length+=$r1[-3];
			if($r1[3]<$r1[4]){
				 $cnv_inf.="$r1[2]:$r1[3]-$r1[4];";
		 	}else{
				$cnv_inf.="$r1[2]:$r1[4]-$r1[3];";
			}
			#	print $cnv_inf."\n";

		}else{
			$flag=$flag_R1{$i}>$flag?$flag_R1{$i}:$flag;
			}
	#	print "length....$flag....$length....".$read_length1."\n";
		
   }
   # print "Other:R1".($read_length1-$ReadLenth)."take from....$length\n";
  $length-=(($read_length1>$ReadLenth)?($read_length1-$ReadLenth):0);

   $read_length1=0;
   for (my $i=0;$i<@read_2;$i++){
			my @r1=split /\t/,$read_2[$i];
			$read_length1+=$r1[-3];
			$length+=$r1[-4];
		if(!exists($flag_R2{$i})){
			$length+=$r1[-3];
			if($r1[3]<$r1[4]){
				 $cnv_inf.="$r1[2]:$r1[3]-$r1[4];";
		 	}else{
				$cnv_inf.="$r1[2]:$r1[4]-$r1[3];";
			}
			
			#		print $cnv_inf."\n";
		}else{
			$flag=$flag_R1{$i}>$flag?$flag_R1{$i}:$flag;
		}
	#	print "length....$flag....$length....".$read_length1."\n";
   }
   # print "Other:R2".($read_length1-$ReadLenth)."take from....$length\n";
   $length-=(($read_length1>$ReadLenth)?($read_length1-$ReadLenth):0);
   if($flag==0){$length=150};
   if($con_inf eq ""){$con_inf.="R10-R1".(@read_1-1)."--R20-R2".(@read_2-1).";";}
	return ($flag,$length,$con_inf,$cnv_inf);
}


