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
my $l2 = "";
my $l3 = "";
my $l4 = "";
my @temp;

my $counter = 0; 
while ($l1 = <Fastq>)
{
	
	$counter = $counter+1;
	$l2 = <Fastq>;
	$l3 = <Fastq>;
	$l4 = <Fastq>;
	chomp $l1; 
	chomp $l2; 
	chomp $l4;
	#my $name = substr($l1, 1); 
	my @t = split ' ', $l1;
	print "$t[0]	0	*	0	*	*	*	0	0	$l2	$l4	\n";
}
#Call should be GFFfile, FastaReff, SNPFilePath

