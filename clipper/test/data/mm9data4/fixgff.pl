use warnings;
use strict;

open(FI, "<".$ARGV[0]);
while (<FI>)
{
   chomp;	
   my $line = $_;
   if (	$line =~ /mRNA/)
   {
      $line =~ s/;/;Parent=/g;
   }
   print $line."\n";
}
