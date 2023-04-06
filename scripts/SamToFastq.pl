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
	@temp = split(/\t/, $l1);	
	$counter = $counter+1;
	if (length($temp[9] > 25))
	{
		print "\@$temp[0]\n";
		print "$temp[9]\n";
		print "+\n";
		print "$temp[10]\n";
	}
	
}
#Call should be GFFfile, FastaReff, SNPFilePath

