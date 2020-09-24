#!/usr/bin/perl -w
use strict;
use warnings;

#----------------------------------------------------------
sub Usage{
	print "perl shiftReads4ATACseq.pl <in.sam> <shift.sam>\n";
	print "\t<in.sam>		input sam file\n";
	print "\t<shift.sam>		output shifted sam file\n\n";

	print "\tFunction: shift alignments on sense strand 4bp, and shift alignments on anti-sense strand -5bp, to represent the center of of the Tn5 transposon binding event. Alternatively, use \"awk 'BEGIN {FS=\"\\t\";OFS=\"\\t\";} \$9>0 {\$4+=4;\$8=\$8+\$6+4;print \$0;next;} \$9<0 {\$4-=5;\$8=\$8+\$6-5;print \$0;next;}' \"\n";
	print "\tQirui Zhang (qirui.zhang\@med.lu.se)\n\t17-02-2020\n\n";
}

#----------------------------------------------------------
my ($sam, $shift)=@ARGV;
if (@ARGV!=2){
	Usage;
	exit;
}

#----------------------------------------------------------
open (IN, "$sam") or die "$!\n";
open (OUT, "> $shift") or die "$!\n";
while (<IN>){
	chomp;
	if(/^@/){
		print OUT "$_\n";
	}else{
		my @tmp=(split/\s+/,$_);
		my $match=$tmp[5];
		$match=~s/[0-9]{1,}[A-LN-Z]//g;
		$match=~s/[M]$//g;
		my @tmp2=(split/M/,$match);
		my $sum=0;
		foreach my $item (@tmp2){
			$sum+=$item;
		}

		if($tmp[8]>0){
			$tmp[3]+=4;
			$tmp[7]+=$sum;
			$tmp[7]+=4;
		}else{
			$tmp[3]-=5;
			$tmp[7]+=$sum;
			$tmp[7]-=5;
		}

		my $last_item=pop @tmp;
		foreach my $item (@tmp){
			print OUT "$item\t";
		}
		print OUT "$last_item\n";
	}
}
close IN;
close OUT;

