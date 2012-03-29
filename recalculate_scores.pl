#!/usr/bin/perl
use warnings;
use strict;
my $minlen = 10;
my $species = $ARGV[1];

open(FI, "<$ARGV[0]");

while (<FI>)
{
   chomp;
   my ($chr, $left, $right, $name, $score, $strand) = split (/\t/, $_);
   $chr =~ s/chr//g;
   my $HEIGHTref = &get_phast($species, $chr, $left, $right);
   my %HEIGHT = %$HEIGHTref;

   my @scores;
   for (my $l =$left; $l<$right;$l++) #extract phastcons within region
   {
      if ($HEIGHT{$l})
      {
	 push @scores, $HEIGHT{$l};
      }else
      {
	 push @scores, 0;
      }
   }
   $DB::single=1;
   my $len = $right - $left + 1;
   my $d = scalar @scores;
   if ($len < $minlen || $d < 3)
   {
#      print STDERR "region too short ($len nt) or there aren't enough data points to consider ($d):\nchr".$chr."\t".$left."\t".$right."\t".$gid."_".$region_n.."\n\n";
      next;
   }
   my ($n,$median,$mean,$average_deviation, $standard_deviation, $variance, $skew, $kurtosis) = &moments(@scores);
   $DB::single=1 if ($n eq "break");
   die "left > right!\n" if ($left > $right);
   print "chr".$chr."\t".$left."\t".$right."\t".$name."\t".$mean."\t".$strand."\n";
}


sub get_phast{
   my ($species, $chr, $start, $stop) = @_;
   my %HASH;
   my $index;
   if ($species eq "hg18")
   {
      $index = "/nas3/yeolab/Conservation/phastCons/hg18_28way/index";
   }elsif ($species eq "hg19")
   {
      $index = "/nas/nas0/yeolab/Conservation/phastCons/hg19_46way/vertebrate/reformat/index";
   } elsif ($species eq "mm9")
   {
      $index = "/nas3/yeolab/Conservation/phastCons/mm9_30way/vertebrate/index";
   }else
   {
      die "I don't know what index to use for $species";
   }
   
#   open(PH, "perl /nas3/yeolab/Conservation/scripts/get_phastcons_fromindex.pl --index $index --chr $chr --from $start --to $stop 2>/dev/null |" )||die;
   open(PH, "perl /nas3/yeolab/Conservation/scripts/get_phastcons_fromindex.pl --index $index --chr $chr --from $start --to $stop |" )||die;
   while (<PH>)
   {
      chomp;
      my @a = split (/\t/, $_);
      $HASH{$a[0]} = $a[1];
   }
   return (\%HASH);
}




sub moments {
   my @nums = @_;
   my $sum = 0;
   my $n = scalar(@nums);
   if ($n < 3)
   {
      return "break";
   }
   foreach my $nn (@nums) {
      $sum += $nn;
   }


   my $mean = $sum/$n;
   my $average_deviation = 0;
   my $standard_deviation = 0;
   my $variance = 0;
   my $skew = 0;
   my $kurtosis = 0;
   if ($n > 3)
   {
      

   foreach (@nums) {
      my $deviation = $_ - $mean;
      $average_deviation += abs($deviation);
      $variance += $deviation**2;
      $skew += $deviation**3;
      $kurtosis += $deviation**4;
   }
   $average_deviation /= $n;
   $variance /= ($n - 1);
   
   $standard_deviation = sqrt($variance);
   
   if ($variance) {
      $skew /= ($n * $variance * $standard_deviation);
      $kurtosis = $kurtosis/($n * $variance * $variance) - 3.0;
   }
   }
   
   @nums = sort { $a <=> $b } @nums;
   my $mid = int($n/2);
   my $median = ($n % 2) ? $nums[$mid] : ($nums[$mid] + $nums[$mid-1])/2;
   
   #printf("n:                  %d\n", $n);
   #printf("median:             %f\n", $median);
   #printf("mean:               %f\n", $mean);
   #printf("average_deviation:  %f\n", $average_deviation);
   #printf("standard_deviation: %f\n", $standard_deviation);
   #printf("variance:           %f\n", $variance);
   #printf("skew:               %f\n", $skew);
   #printf("kurtosis:           %f\n", $kurtosis);
   
   return($n,$median,sprintf('%.2f',$mean),$average_deviation, $standard_deviation, $variance, $skew, $kurtosis);

}

