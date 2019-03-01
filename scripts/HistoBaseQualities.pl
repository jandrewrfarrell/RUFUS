#!/usr/bin/perl

use strict;
my $line; 
while ($line = <>)
{
	my @fields = split(" ", $line); 
	my @char = split("", $fields[10]);
	my $b =0; 
	for ( $b = 0;  $b < scalar @char; $b++)
	{
		my $ord = ord($char[$b])-33; 
		print "$ord\n";
	}
}


