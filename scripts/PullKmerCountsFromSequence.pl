#!/usr/bin/perl


use strict;
my $kgJhash = "/scratch/ucgd/lustre/work/u0991464/RUFUS.1000g.reference/1000G.RUFUSreference.min45.Jhash";
my $jellyfishPath="~/bin/RUFUS//bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
my $string = $ARGV[0];
my $HS = $ARGV[1];
print "sequence\t$ARGV[2]\t"; 
for (my $j = 3; $j < scalar(@ARGV); $j++)
	        {print "$ARGV[$j]\t";}
		print "1kg\n";

for ( my $i = 0; $i < length($string) - $HS; $i++)
{

	my $hash =  substr($string, $i, $HS); 
	print "$hash\t";	
	my $first = `$jellyfishPath  query $ARGV[2] $hash`;
        chomp $first; 
        my @temp1 = split / /, $first;
        print "$temp1[1]\t";

	for (my $j = 3; $j < scalar(@ARGV); $j++)
	{
		my $first = `$jellyfishPath  query $ARGV[$j] $hash`;
		chomp $first; 
		my @temp1 = split / /, $first;
		print "$temp1[1]\t";
	}

	my $first = `$jellyfishPath  query $kgJhash $hash`;
	chomp $first;
	my @temp1 = split / /, $first;
	if ($temp1[1] eq 0){
		my $revcomp = reverse $hash;
		$revcomp =~ tr/ATGCatgc/TACGtacg/;
		$first = `$jellyfishPath  query $kgJhash $revcomp`;
		chomp $first;
		@temp1 = split / /, $first;
	}
	print "$temp1[1]\t";
	print "\n"; 
}
print "\n";
