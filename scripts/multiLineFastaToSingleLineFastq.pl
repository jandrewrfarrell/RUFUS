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

my $l;
my $seq = "";
$l = <Fastq>;
my @a = split(" ", $l);
my @len = split("=", $a[1]); 
my @reads = split("=", $a[2]); 

print "$a[0]_L$len[1]_D$reads[1]:5:5\n";
while ($l = <Fastq>)
{
	chomp($l);
	my $FC = substr($l , 0, 1);
	if ($FC eq ">")
	{
	#	print "yay header\n";
		print "$seq\n";
		print "+\n";
		print "$seq\n";
		

		my @a = split(" ", $l);
		my @len = split("=", $a[1]);
		my @reads = split("=", $a[2]);
		print "$a[0]_L$len[1]_D$reads[1]:5:5\n";	
		
		#print "$l\n";
		$seq = ""; 
	}
	else
	{
		$seq=$seq . $l;	
	} 
	
}
print "$seq\n";
print "+\n";
print "$seq\n";

