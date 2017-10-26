#!/usr/bin/perl

use strict;

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	chomp($string);
	return $string;
}


open (Fastq , $ARGV[0]) || die "ERROR could not open sam file";

my $l1 = "";
my $lastLine = "";
my @temp;
my $chr = "nope"; 
my $pos = "nope";
my $ref = "nope"; 
my $alt = "nope";
my $sample = "nope"; 
my $counter = 0;
 
while ($l1 = <Fastq>)
{
	my $FC = substr($l1 , 0, 1);
	if ($FC == "#")
	{
		print "$l1";
	}
	else
	{
		@temp = split(/\t/, $l1);	
		$counter = $counter+1;
		chomp ($l1);
		if ($temp[0] == $chr && $temp[1] == $pos && $temp[3] == $ref && $temp[4] == $alt)
		{
				#print "$l1\n+++$lastLine\n"; 
		}
		else
		{
			print "$l1\tUNIUQE\n";
			$chr = $temp[0];
			$pos = $temp[1];
			$ref = $temp[3];
			$alt = $temp[4];
		}		
		$lastLine = $l1; 
	} 
	
}

