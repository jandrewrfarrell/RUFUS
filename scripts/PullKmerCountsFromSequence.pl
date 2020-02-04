#!/usr/bin/perl



use strict;
my $jellyfishPath="~/bin/RUFUS//bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
my $string = $ARGV[0];
my $HS = $ARGV[1];
print "sequence\t$ARGV[2]\t"; 
for (my $j = 3; $j < scalar(@ARGV); $j++)
	        {print "$ARGV[$j]\t";}
		print "\n"; 
for ( my $i = 0; $i < length($string) - $HS; $i++)
{

	my $hash =  substr($string, $i, $HS); 
		
		my $first = `$jellyfishPath  query $ARGV[2] $hash`;
                chomp $first; 
                my @temp1 = split / /, $first;
                print "$temp1[0]\t$temp1[1]\t";

	for (my $j = 3; $j < scalar(@ARGV); $j++)
	{
		my $first = `$jellyfishPath  query $ARGV[$j] $hash`;
		chomp $first; 
		my @temp1 = split / /, $first;
		print "$temp1[1]\t";
	}
	print "\n"; 
}
print "\n";
