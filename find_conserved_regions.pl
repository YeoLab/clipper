#!/nas3/yeolab/Software/perl-5.12.0/bin/perl
use warnings;
use strict;
use Getopt::Long;
use IO::Handle;
use Fcntl ":seek";
use Fcntl ":flock";
use lib '/nas3/yeolab/Software/EXTRA_perlmodules/Parallel-ForkManager-0.7.5/blib/lib';
use Parallel::ForkManager;
use Benchmark;
use File::Temp qw/ tempfile tempdir/;
my $t_start = Benchmark->new();
use Cwd;
my $cwd = getcwd;
$cwd = $cwd."/";
my $species = '';		#species (string)
my $outdir = $cwd;		#prefix (string)
my $AS_ID = "";	#AS_ID gid1,gid2...	#only do peakfinding on a particular AS.STRUCTURE gene, comma delimited.
my $usechr = ""; #chr (chr1,chr2...)	#only use a particular chromosome.
my $children = 6; #childs (int > 1)	#number of parallel child processes allowed to run simultaneously (genes analyzed at once)
my $serial=0;
my $margin =1000;
my $minlen = 20;
my $dontdelete = 0;
my $usage = 0;
my $floor = .9;
my $ceiling = 1.0;
my $tmpdir;
my $index;

GetOptions(
	   'species=s' => \$species,
	   'prefix=s' => \$outdir,
	   'chr=s' => \$usechr,
	   'childs=s' => \$children,
	   'serial' => \$serial,
	   'ceiling=f' => \$ceiling,
	   'floor=f' => \$floor,
	   'dontdelete|keep-tmp' => \$dontdelete,
	   'usage' =>\$usage,
	   'margin=s' => \$margin,
	   'minlen=s' => \$minlen,
	   'tmpdir=s' => \$tmpdir,
	   'index=s' => \$index
	  );

if ($species eq "hg18")
{
   $index = "/nas3/yeolab/Conservation/phastCons/hg18_28way/index" unless $index;
}elsif ($species eq "hg19")
{
   $index = "/nas3/yeolab/Conservation/phastCons/hg19_46way/vertebrate/reformat/index" unless $index;
} elsif ($species eq "mm9")
{
   $index = "/nas3/yeolab/Conservation/phastCons/mm9_30way/vertebrate/index" unless $index;
}else
{
   die "I don't know what index to use for $species" unless $index;
}

print STDERR "Using index: $index\n";

#my $tmpdir = $outdir."tmp/";
my $cleanup = 1;
if ($dontdelete || $tmpdir)
{
   $cleanup = 0;
}
unless (-e "$outdir")
{
   mkdir("$outdir");
}
my %OPTS;
$OPTS{"DIR"} = $outdir;
$OPTS{"CLEANUP"}= $cleanup;
unless ($tmpdir)
{
   $tmpdir = File::Temp->newdir('tmp_XXXXXXX', %OPTS);
}
unless ($tmpdir =~ /\/$/)
{
   $tmpdir .= "/";
}


print "Using a temporary directory: $tmpdir\n";
#print "starting";
if (!$AS_ID && $dontdelete)
{
   print STDERR "Are you sure you want to use the --dontdelete option when AS_ID is not used?  This may create too many files.\nSay \"yes\" if you want to do it anyways...\n";
   my $response = <STDIN>;
   chomp $response;
   unless ($response eq "yes")
   {
      print STDERR "not using --dontdelete\n";
      $dontdelete=0;
   }
}
package other;
use POSIX ":sys_wait_h";	#imports WNOHANG

my $pm = new Parallel::ForkManager($children);
my @AS_IDS = split (/\,/, $AS_ID);
if ($usage)
{
   print STDERR "must specify ingenes, species\n";
   print "Must be run with options, the options are:
#deprecated	species (string)	#required. supported:hg18,mm8,ce6
	prefix (string)		#output location, default:cwd
#deprecated	AS_ID gid1,gid2,...	#only do peakfinding on a particular list of AS.STRUCTURE genes, comma delimited. default:OFF (all gids in ingenes considered)
	chr (chr1,chr2...)	#only use a particular chromosome. default:OFF
        childs (int > 1)	#number of parallel child processes allowed to run simultaneously (genes analyzed at once), default:6
        serial (flag)           #do not run child processes, run genes one-after-another. default:OFF
        margin (int)            #number of bases with floor<phastcons < cieling allowed before starting a new regions. default:1000
        ceiling (float)		#upper limit for \" conservation\"
        floor (float)		#lower limit for \" conservation\"
        index (string) 		#index to use. there is a default for mm9,hg18,hg19,ce6 but you'll need to set this manally for other species
";
	exit(0);
}
my @chr;
if ($usechr)
{
   $usechr =~ s/chr//g;
   @chr= split (/\,/, $usechr);
} else
{
   if ($species eq "hg18" || $species eq "hg19")
   {
      @chr = (1..22, "X", "Y");
   }
   if ($species eq "hg18CMV")
   {
      die "no phastcons for $species\n";
   @chr = (1..22, "X", "Y", "CMV");
   }
   if ($species eq "mm8" ||$species eq "mm9")
   {
      @chr = (1..19, "X", "Y");
   }
   if ($species eq "ce6")
   {
      die "no phastcons for $species\n";
      @chr = ("I", "II", "III", "IV", "V", "X");
   }
}
#my ($aliref, $giref) = &build_alignable_mRNA($species, \@chr);
#my %ALI = %$aliref;
#my %GENEINFO = %$giref;

#my $tmpdir = $outdir."tmp/";
unless (-e $tmpdir)
{
   mkdir($tmpdir);
}
my %kids;
#my @GIDS =  keys %GENEINFO;
#my $total = scalar @GIDS;
#my $go = 0;
my $iteration =0;
my %checked; #for some reason, the queue insists on submitting some genes multiple times, restrict that. fixed, Proc::Queue does not work.  use Parallel::Forkmanager instead.
my @GIDS;


open(GENEBED, "</nas3/lovci/projects/ucscBED/".$species."/".$species.".AS.STRUCTURE_genes.BED");

while (<GENEBED>)
{
   chomp;
   my ($chr, $g_start, $g_stop, $gid, $score, $signstrand) = split (/\t/, $_);

   $chr =~ s/chr//g;
   my $CONT = 0;
   if ($AS_ID)
   {
      foreach my $AS (@AS_IDS)
      {
	 if ($gid eq $AS)
	 {
	    $CONT=1;
	 }
      }
   } else
   {
      $CONT=1;
   }
   next if $CONT eq 0;
   $CONT=0;
   if ($usechr)
   {
      foreach my $tchr (@chr)
      {
	 if ($tchr eq $chr)
	 {
	    $CONT=1;
	 }
      }
   } else
   {
      $CONT=1;
   }
   next if $CONT eq 0;


#   $iteration++;
#   if ($iteration %100 eq 0)
#   {
#      my $percent_complete = sprintf('%6.3f',100*$iteration/$total);
#      print STDERR $iteration." iterations out of ".$total.", ".$percent_complete."% complete\n";
#   }

   push @GIDS, $gid;
   unless ($serial)
   {
      my $pid = $pm->start and next; #here is where a child process starts.  The parent will continue to spit out children and keep $children many alive at once. so sad.
      die "$0: fork: $!" unless defined $pid; #clearly lifted      
   }
   my $tmpfile = $tmpdir.$gid.".tmp";
#   next unless ($gid eq "ENSG00000184221");
   open(REGIONS, ">".$tmpfile.".regions") ||die;

   my $HEIGHTref = &get_phast($index, $species, $chr, $g_start, $g_stop);
   my %HEIGHT = %$HEIGHTref;
   my $region_n = 0;
   my $loc = $g_start;
   while ($loc <= $g_stop)
   {
#      print $loc."\n";
      if ($HEIGHT{$loc} && $HEIGHT{$loc} > $floor && $HEIGHT{$loc} <= $ceiling)
      {
#	 open(TMP, ">".$tmpfile.".dat") ||die;

	 my $gap = 0;
	 my $left = $loc;
#	 until ($gap == $margin || $region_n == 1) # make a gap $margin wide around reads and fit only chunks of the data.
#	 {
#	    if ($HEIGHT{$left} && $HEIGHT{$left} > $floor && $HEIGHT{$loc} <= $ceiling)	    {
#	       $gap=0;
#	    }
#	    elsif ($left eq $g_start)
#	    {
#	       last;
#	    }else
#	    {
#	       $gap++;
#	    }
#	    $left--;
#	 }
#	 $left = $left +$gap;
	 $gap = 0;
	 until ($gap == $margin)
	 {
	    if ($HEIGHT{$loc} && $HEIGHT{$loc} > $floor && $HEIGHT{$loc} <= $ceiling)
	    {
	       $loc++;
	       $gap=0;
	    }
	    elsif ($loc eq $g_stop)
	    {
	       last;
	    }
	    else
	    {
	       $loc++;
	       $gap++;
	    }
	 }
	 my $right = $loc - $gap;
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
	 my $len = $right - $left + 1;
	 my $d = scalar @scores;
	 if ($len < $minlen || $d < 3)
	 {
	    print STDERR "region too short ($len nt) or there aren't enough data points to consider ($d):\nchr".$chr."\t".$left."\t".$right."\t".$gid."_".$region_n.."\n\n";
	    next;
	 }


	 my ($n,$median,$mean,$average_deviation, $standard_deviation, $variance, $skew, $kurtosis) = &moments(@scores);
	 $DB::single=1 if ($n eq "break");
	 die "left > right!\n" if ($left > $right);
	 $region_n++;
	 print REGIONS "chr".$chr."\t".$left."\t".$right."\t".$gid."_".$region_n."\t".$mean."\t".$signstrand."\n";
#	 print "chr".$chr."\t".$left."\t".$right."\t".$gid."_".$region_n."\t".$mean."\t".$signstrand."\n";

      }
   }continue
   {
      $loc++;
   }
   unless ($serial)
   {
      exit(0);
      $pm->finish;      
   }
}
unless ($serial)
{
   $pm->wait_all_children;
}
my $allfile;
$ceiling = sprintf("%.0f", $ceiling*10);
$floor = sprintf("%.0f", $floor*10);
if ($usechr)
{
   if ($AS_ID)
   {
      open ($allfile, ">".$outdir."conserved_regions.PARTIAL.".$species.".chr".$usechr.".margin_".$margin.".".$floor."-".$ceiling.".BED") ||die;
   }else
   {
      open ($allfile, ">".$outdir."conserved_regions.".$species.".chr".$usechr.".margin_".$margin.".".$floor."-".$ceiling.".BED") ||die; 
   }
}else
{
   if ($AS_ID)
   {
      open ($allfile, ">".$outdir."conserved_regions.PARTIAL.".$species.".margin_".$margin.".".$floor."-".$ceiling.".BED") ||die;
   }else
   {
      open ($allfile, ">".$outdir."conserved_regions.".$species.".margin_".$margin.".".$floor."-".$ceiling.".BED") ||die;
   }
}

printf $allfile "track name=".$species.".conserved description=\"highly conserved regions in $species genome\" itemRgb=On\n";

GID: foreach my $gid (@GIDS)
{
   my $last_loc = 0;
   my $finish =0;
   my $left_side = 0;
   my @all_clusters;
   my $tmpfile = $tmpdir.$gid.".tmp.regions";
   open(EDGES, "<".$tmpfile) || next;
   unless ($dontdelete)
   {
      system("rm -rf $tmpfile");
   }
   my @edges = <EDGES>;
   system("rm -rf ".$tmpfile.".regions") unless ($dontdelete);
   close(EDGES);
   for (my $i = 0; $i< scalar(@edges); $i++)
   {
      chomp $edges[$i];
      print $allfile $edges[$i]."\n";
   }
}

my $t_finish = Benchmark->new;
my $timetotal = Benchmark::timediff($t_finish,$t_start);



print STDERR "Finished with find_conserved_regions on $usechr, it took ".Benchmark::timestr($timetotal)."\n";

sub get_phast{
   my ($index, $species, $chr, $start, $stop) = @_;
   my %HASH;
   $DB::single=1;
#   open(PH, "perl /nas3/yeolab/Conservation/scripts/get_phastcons_fromindex.pl --index $index --chr $chr --from $start --to $stop |" )||die;
   open(PH, "perl /nas3/yeolab/Conservation/scripts/get_phastcons_fromindex.pl --index $index --chr $chr --from $start --to $stop 2>/dev/null |" )||die;

   while (<PH>)
   {
      chomp;
      my @a = split (/\t/, $_);
      $HASH{$a[0]} = $a[1];
   }
   return (\%HASH);
}






sub fac {
   my $n = shift;
   my $fact = $n;
   if ($fact ==1)
   {
      return(1);
   }
   if ($n ==0)
   {
      print "you suck, you can't get 0!";
      die;
   }
   for (my $i = $n-1; $i>0; $i--)
   {
      $fact *=$i;
   }
   return($fact);
}


sub log10 {
   my $n = shift;
   return log($n)/log(10);
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

