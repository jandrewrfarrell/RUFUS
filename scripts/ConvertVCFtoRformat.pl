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
my @samples; 
my $lines = 0; 
my @fields; 
while ($l = <STDIN>)
{
	chomp($l); 
	@temp = split("\t", $l);
	if ($temp[0] =~ /^#/)
	{
		if ($temp[0] eq "#CHROM")
		{
			for (my $i = 9; $i < @temp; $i++)
			{
				push @samples, $temp[$i]; 
			}
		}
	}
	else
	{
		if ($lines ==0)
		{
			@fields = split(":", $temp[8]); 
			print	"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"; 
			for (my $i = 0; $i < @samples; $i++)
			{
				for (my $j = 0; $j < @fields; $j++)
				{
					print "	$samples[$i]-$fields[$j]";
				}
			}
			print "\n"; 
			$lines = 1; 
		}
		for (my $i = 0; $i <= 7; $i++)
		{
			print "$temp[$i]	"
		}
		print "$temp[8]"; 
		for (my $i = 9; $i < @temp; $i++)
		{
			@fields = split(":", $temp[$i]);
			for (my $j = 0; $j < @fields; $j++)
			{
				print "	$fields[$j]";
			}
		}
		print "\n"; 
	}
}

