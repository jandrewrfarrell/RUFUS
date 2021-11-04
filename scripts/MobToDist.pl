#!/usr/bin/perl

use strict;

open (Fastq , $ARGV[0]) || die "ERROR could not open sam file";

my $l;
$l = <Fastq>;
print "$l"; 
$l = <Fastq>;
print "$l";
$l = <Fastq>;
print "$l";
$l = <Fastq>;
print "$l";
my @fa; 
$l = <Fastq>;
print "$l";
chomp($l);
@fa = split(" ", $l); 


my @counts;
my $lines = 0; 
while ($l = <Fastq>)
{
	$lines++; 
	chomp($l);
	my @a = split(" ", $l);
	push @counts, @a; 
	
}
my $total; 
for (my $i = 1; $i < $lines; $i++)
{
	print " i = $i \n"; 
	print " adding $i $counts[3][3]"; 
	#$total += $counts[$i][1]; 
}
print "total = $total"; 
