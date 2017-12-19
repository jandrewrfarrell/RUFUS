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
my @temp;
my @scores; 
for (my $i = 0; $i < 30; $i ++)
{
 	$scores[$i] = 0;	
}
my $counter = 0; 
while ($l1 = <Fastq>)
{
	
	if (substr ( $l1, 1, 1) eq "#")
	{}
	else
	{
			
		@temp = split("\t", $l1);
		$scores[$temp[5]]++; #  = $scores[$temp[6]] +1;
	}
			
	 
	
}
for (my $i = 0; $i < 30; $i ++)
{
	print "$i; ";
        for (my $j = 0; $j < $scores[$i]; $j++)
	{
		print "+"; 
	}
	print "; $scores[$i] \n";
}
print "~~~~~~~~~\n";
for (my $i = 30; $i < 10000; $i ++)
{
	if ( $scores[$i] > 0)
        {
		print "$i; ";
        	for (my $j = 0; $j < $scores[$i]; $j++)
        	{
        	        print "+";
        	}
        	print "; $scores[$i] \n";
	}
}
#Call should be GFFfile, FastaReff, SNPFilePath

