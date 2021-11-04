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
my $counter = 0; 
while ($l1 = <Fastq>)
{
	
	if (substr ( $l1, 1, 1) eq "#")
	{}
	else
	{
			
		my @temp = split("\t", $l1);
		my @info = split(";", $temp[7]); 
		$temp[0] =~ s/chr//;
		if ($temp[4] =~/\</)
		{
			my $end = -1; 
			for (my $i = 0; $i < @info; $i++)
			{
				my @f = split("=", $info[$i]);
				if ("$f[0]" eq "END")
				{
					$end = $f[1]; 
				}
			}
			my $start = $temp[1] - 1;
			print "$temp[0]	$start	$end	$temp[2]-$temp[5]\n";
		}
		elsif (abs(length($temp[3]) - length($temp[4])) > 50)
		{
			my $start = $temp[1] - 1;
			if ( $temp[2] =~ /Y/ )
			{
				$start = $start - length($temp[4]);
			}

			my $end = $temp[1] - 1 + length($temp[3]);
			print "$temp[0]	$start	$end	$temp[2]-$temp[5]\n"; 
		}
	}
			
	 
	
}
#Call should be GFFfile, FastaReff, SNPFilePath

