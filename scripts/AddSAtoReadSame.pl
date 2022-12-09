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



my @temp;
my $l1; 
my $counter = 0; 
while ($l1 = <>)
{
	my $SA=0; 
	@temp = split(/\t/, $l1);	
	for (my $i = 10; $i < @temp; $i++)
	{
		if (substr($temp[$i], 0, 4) eq "SA:Z")
		{
			my @temp2 = split(/\;/, $temp[$i]);
			$SA=@temp2;
		}
	}
	print "$temp[0]:SA=$SA";
	
	for (my $i = 1; $i < @temp; $i++)
	{
		print "\t$temp[$i]";
	}

}
#Call should be GFFfile, FastaReff, SNPFilePath

