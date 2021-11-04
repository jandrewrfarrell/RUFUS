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


my $l; 
my @temp;
my $last = "first"; 
while ($l = <STDIN>)
{
	chomp($l); 
	@temp = split("\t", $l);
	if ($temp[0] =~ /^#/)
	{
		print "$l\n";
	}
	else
	{
		if ($last eq "first")
		{
			print "$l\n"; 
		}
		else
		{
			my @lastTemp = split("\t", $last); 
			if ($lastTemp[0] eq $temp[0] && $lastTemp[1] eq $temp[1] && $lastTemp[2] eq $temp[2] && $lastTemp[3] eq $temp[3] && $lastTemp[4] eq $temp[4])
			{
				####I should do something better here, merge the calls or use the one with the higher value
			}
			else
			{
				print "$l\n";
			}
		}
		$last = $l;
	}
}



