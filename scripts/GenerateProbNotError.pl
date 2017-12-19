#!/usr/bin/perl 


use strict;
use Math::BigFloat ':constant'


sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	chomp($string);
	return $string;
}


open (Fastq , $ARGV[0]) || die "ERROR could not open sam file";
my $l1; 
my @temp;
$l1 = <Fastq>;
$l1 = <Fastq>;
$l1 = <Fastq>;
$l1 = <Fastq>;
while ($l1 = <Fastq>)
{
	
			
		@temp = split("\t", $l1);
		my $num = 0;
		my $denom=0; 
		for (my $i =3; $i<length(@temp); $i++)
		{
			$num+=$temp[$i];
		}
		
		my $value = $num/($num+$temp[1]); 
		print "$temp[0]\t$num\t$temp[$0]\t$value\n";

					
	 
	
}
#Call should be GFFfile, FastaReff, SNPFilePath

