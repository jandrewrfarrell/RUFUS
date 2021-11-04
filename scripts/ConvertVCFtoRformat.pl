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
my %SVids; 
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
			print	"CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	SIZE	TYPE	COMPLEX	FORMAT"; 
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
		#####size stuff goes here ########
		my $size = 0; 
		my $type = "none";
		my $complex = "no" ; 
		if ( index($temp[4], "<DEL>") != -1)
		{
			my $pull1 = substr($temp[7], index($temp[7], "SVLEN=") +6 ); 
			$size = substr($pull1, 0, index($pull1, ";"));
			$type = "del";
			if (length $temp[4] > 5)
			{
				$complex = "yes"; 
			}

		}
		elsif ( index($temp[4], "<INS>") != -1)
		{
			my $pull1 = substr($temp[7], index($temp[7], "SVLEN=") +6 );
                        $size = substr($pull1, 0, index($pull1, ";"));
			$type = "INS";
                        if (length $temp[4] > 5)
                        {
                                $complex = "yes";
                        }		
		}
		elsif ( index($temp[4], "<DUP>") != -1)
                {
                        my $pull1 = substr($temp[7], index($temp[7], "SVLEN=") +6 );
                        $size = substr($pull1, 0, index($pull1, ";"));
                	$type = "dup";
                        if (length $temp[4] > 5)
                        {
                                $complex = "yes";
                        }
		}
		elsif ( index($temp[4], "<INV>") != -1)
                {
                        my $pull1 = substr($temp[7], index($temp[7], "END=") +4 );
                        $size = substr($pull1, 0, index($pull1, ";"));
			$size = $size - $temp[1]; 
                        $type = "INV";
			my $des = substr($temp[2], 0, index($temp[2], "-"));
                        if (length $temp[4] > 5 || index($des, "Y") != -1 || index($des, "D") != -1 || index($des, "I") != -1 || index($des, "bnd") != -1)
                        {
                                $complex = "yes";
                        }
                }
		elsif ( index($temp[4], "<INS:ME:MOB>") !=-1)
		{
			$size = 3000000001; 
			$type = "MOB"; 
			$complex = "no"; 	
		}
		elsif ( index($temp[2], "bnd_") != -1 )
		{

			my $pull1 = substr($temp[7], index($temp[7], "SVID=") +5 );
			#my $svid = substr($pull1, 0, index($pull1, ";"));
			my $svid = substr($pull1, 0, 1);
			####fix cuase I forgot the ; in some bnd svid fields
			     
			if ( exists $SVids{$svid} ) 
			{
				$size ="dup";
			}
			elsif (index($temp[2], "OrphanBND") > -1)
			{
				$size = "orphan-$temp[2]"; 
			}
			else
			{
				$size = 3000000002;
				$SVids{$svid} = 1; 
			}
			if (index($temp[2], "_bnd)_") != -1)
			{
				$complex = "yes";
			}
			$type = "trans";

		}
		else
		{
			if ( length($temp[3]) == 1 && length($temp[4]) == 1)
			{
				$size = 0; 
				$type = "snv"; 
				$complex = "no"; 
			}
			elsif(length($temp[3]) == 1)
			{
#ins stuff goes here 
				$size = length($temp[4])-1; 
				my $des = substr($temp[2], 0, index($temp[2], "-"));
				if (index($des, "Y") != -1 )
				{
					$type = "dup";
				}
				elsif (index($des, "I") != -1 )
				{
					$type = "ins"; 
				}
				else
                                {
                                        $type = "ERRORns"; 
                                }

				$complex = "no";
			
                        }
			elsif(length($temp[4]) == 1)
			{
				$size = length($temp[3])-1;
				$size = $size *-1;
                                my $des = substr($temp[2], 0, index($temp[2], "-"));
                                if (index($des, "D") != -1) 
                                { 
                                        $type = "del";   
                                }
				else
				{
					$type = "ERRORdel"; 
				}
                                $complex = "no";
#del stuff goes here 
			}
			elsif(length($temp[4]) ==length($temp[3]))
			{
				$size = length($temp[4])-1;
				my $des = substr($temp[2], 0, index($temp[2], "-"));
				if (index($des, "X") != -1 && index($des, "D") == -1 && index($des, "I") == -1 && index($des, "Y") == -1)
				{
					$size = 0;
					$type = "multisnp"; 
					$complex = "no"; 
				}
				else
				{
					$size = length($3); 
					$type = "complex"; 
					$complex = "yes"; 
				}

			}
			else
			{

###make this beter, should cehck for del dup ins tandem 
				 $size = length($temp[3]);
				 $type = "multi"; 
				 $complex = "yes"; 
				 if (length($4) > $size)
				 {
					 $type = "multi"; 
					 $size = -1*length($4); 
				 }
			}

		}
		print "$size	$type	$complex	"; 
		print "$temp[8]"; 
		for (my $i = 9; $i < @temp; $i++)
		{
			@fields = split(":", $temp[$i]);
			for (my $j = 0; $j < @fields; $j++)
			{
				if ($fields[$j] eq "./.")
				{
					print " 0/0";
				}
				elsif($fields[$j] eq ".")
				{
					print "	0";
				}
				else
				{
					print "	$fields[$j]";
				}
			}
		}
		print "\n"; 
	}
}

