#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Long;
use Cwd;


use File::Spec;



my $job_type;
my $hostname = `hostname`;
our $base = "belongtous";
if ($hostname=~ /optiputer/ || $hostname =~ /compute/)
{
   $base = "nas/nas0";
   $job_type = "SGE";


}elsif ($hostname =~ /tcc/ || $hostname =~ /triton/)
{
   $base = "projects";
   $job_type= "PBS";

}

require "/".$base."/yeolab/Software/perl-5.12.0/lib/File/Temp.pm";
use File::Temp qw/ tempfile tempdir/;
my $this_script = File::Spec->rel2abs($0);

#extract conserved kmers
#input: BED FILE, species, index, nspecies
my $bed;
my $species;
my $index;
my $pivot;
my $s_list;
my $k_size = 6;
my $cutoff = 10;
my $job_name = "XXXX";
my $range;
my $tmpdir;
my $min_length=10; #minimum region length for kmers within that region to be considered
my $rangesize = 100;
my $chunk = 100;
my $start;
my $resources;

my $finish;
my $outfile ="kmers.out";
GetOptions("bed=s" => \$bed,
	   "species=s" => \$species,
	   "index=s" => \$index,
	   "s_list=s" => \$s_list,
	   "pivot=s" => \$pivot,
	   "k=s" => \$k_size,
	   "cutoff=s" => \$cutoff,
	   "range=s" => \$range,
	   "start" => \$start,
	   "chunk=s" => \$chunk,
	   "tmpdir=s" => \$tmpdir,
	   "rangesize=s" => \$rangesize,
	   "outfile|out=s" => \$outfile,
	   "resources=s" =>\$resources,
	   "finish" => \$finish
);

if ($resources)
{
   $resources = join(" ", $resources, "-j", "n");   
}else
{
   $resources = "-j n";
}
if ($species eq "mm9")
{
   $pivot = "mm9_30" unless ($pivot);
   $index = "/".$base."/yeolab/Conservation/mm9_30way" unless ($index);
   $s_list = "mm9,rn4,cavPor2,oryCun1,hg18,panTro2,ponAbe2,rheMac2,calJac1,otoGar1,tupBel1,sorAra1,eriEur1,canFam2,felCat3,equCab1,bosTau3,dasNov1,loxAfr1,echTel1,monDom4,ornAna1,galGal3,anoCar1,xenTro2,tetNig1,fr2,gasAcu1,oryLat1,danRer5" unless ($s_list);
}
if ($species eq "hg18")
{
   $pivot = "hg18_28" unless ($pivot);
   $index = "/".$base."/yeolab/Conservation/hg18_28way" unless ($index);
   $s_list = "hg18-panTro2-rheMac2-otoGar1-tupBel1-mm8-rn4-cavPor2-oryCun1-sorAra1-eriEur1-canFam2-felCat3-equCab1-bosTau3-dasNov1-loxAfr1-echTel1-monDom4-ornAna1-anoCar1-galGal3-xenTro2-fr2-tetNig1-gasAcu1-oryLat1-danRer4" unless ($s_list);


}

if ($species eq "hg19")
{
   $pivot = "hg19_46" unless ($pivot);
   $index = "/".$base."/yeolab/Conservation/hg19_46way" unless ($index);
   $s_list = "hg19-panTro2-gorGor1-ponAbe2-rheMac2-papHam1-calJac1-tarSyr1-micMur1-otoGar1-tupBel1-mm9-rn4-dipOrd1-cavPor3-speTri1-oryCun2-ochPri2-vicPac1-turTru1-bosTau4-equCab2-felCat3-canFam2-myoLuc1-pteVam1-eriEur1-sorAra1-loxAfr3-proCap1-echTel1-dasNov2-choHof1-macEug1-monDom5-ornAna1-galGal3-taeGut1-anoCar1-xenTro2-tetNig2-fr2-gasAcu1-oryLat2-danRer6-petMar1" unless ($s_list);


}

$s_list =~ s/\,/-/g; #comma or dash are acceptable separators
my @species_list = split (/\-/, $s_list); #not used
push @species_list, "++";
push @species_list, "CS";
#hg18 28-way, all species =  hg18,panTro2,rheMac2,otoGar1,tupBel1,mm8,rn4,cavPor2,oryCun1,sorAra1,eriEur1,canFam2,felCat3,equCab1,bosTau3,dasNov1,loxAfr1,echTel1,monDom4,ornAna1,anoCar1,galGal3,xenTro2,fr2,tetNig1,gasAcu1,oryLat1,danRer4

my $ns = (scalar @species_list);
$ns-=3;  #remove self, ++ and CS
my $kmerlistfile = $outfile.".k".$k_size.".c".$cutoff.".ns".$ns;
if ($finish)
{
   open(OUT, ">".$kmerlistfile);
   open(OUTMA, ">".$outfile.".multiz");
   open(OUTFA, ">".$outfile.".fa");
   open(OUTCONS, ">".$outfile.".cons.fa");
   open(OUTNONCONS, ">".$outfile.".noncons.fa");
   my $list = $tmpdir."scripts_out/list";
   open(LIST, "<".$list) || die "counldn't read $list because $!";
   my %K;
   while (<LIST>)
   {
      chomp;
      my $r = $_;
      my $file = $tmpdir."trash/file.".$r;
      my $fileMA = $tmpdir."trash/multiz.".$r;
      my $fileFA = $tmpdir."trash/fa.".$r;
      my $filecons = $tmpdir."trash/cons.".$r;
      my $filenoncons = $tmpdir."trash/noncons.".$r;


      if (-e $filecons)
      {	
	 open(FI, "<".$filecons);
	 while (<FI>)
	 {
	    print OUTCONS $_;
	 }
	 close(FI); 
      }

      if (-e $filenoncons)
      {	
	 open(FI, "<".$filenoncons);
	 while (<FI>)
	 {
	    print OUTNONCONS $_;
	 }
	 close(FI); 
      }

      if (-e $fileMA)
      {
	 open(FI, "<".$fileMA);
	 while (<FI>)
	 {
	    print OUTMA $_;
	 }
	 close(FI);
      }


      if (-e $fileFA)
      {
	 open(FI, "<".$fileFA);
	 while (<FI>)
	 {
	    print OUTFA $_;
	 }
	 close(FI);
      }

      if (-e $file)
      {
	 open(FI, "<".$file);
	 while (<FI>)
	 {
	    chomp;
	    my ($k, $c, $nc) = split (/\t/, $_);
	    $K{$k}{conserved} +=$c;
	    $K{$k}{notconserved} +=$nc;
	 }
#	 system("rm -rf $file");
      }
   }
   foreach my $k (keys %K)
   {
      print OUT $k."\t".$K{$k}{conserved}."\t".$K{$k}{notconserved}."\n";
   }
#   system("rm -rf $tmpdir");
   exit;
}

if ($start)
{
   $bed = File::Spec->rel2abs($bed);
   my $cleanup = 0;
   my $outdir = File::Spec->rel2abs(getcwd);
   my %OPTS;
   $OPTS{"DIR"} = $outdir;
   $OPTS{"CLEANUP"}= $cleanup;
   if ($bed =~ /gz$/)
   {
      open(BED, "gunzip -c $bed |") ||die;
   }else
   {
      open(BED, "<$bed") ||die;
   }
   my @a = <BED>;
   my $nlines = scalar @a;
   unless ($tmpdir)
   {
      $tmpdir = File::Temp->newdir('tmp_kmXXXX', %OPTS);
      $tmpdir = File::Spec->rel2abs($tmpdir);
      $tmpdir .= "/";
   }
   unless (-e $tmpdir)
   {
      mkdir($tmpdir);
   }
   if ($job_name eq "XXXX")
   {
      if ($tmpdir =~ m/tmp_(\w+)/)
      {
	 $job_name = $1;
      }
   }
   die if ($job_name =~ /XXXX/);
  print STDERR "Using a temporary directory: $tmpdir\n";
   mkdir($tmpdir."cluster_scripts/");
   mkdir($tmpdir."scripts_out/");
   mkdir($tmpdir."trash/");
   my $argfile = $tmpdir."cluster_scripts/run_args.txt";
   open(ARGS, ">".$argfile) || die "counldn't write to ".$tmpdir."cluster_scripts/run_args.txt because $!";
   my $list = $tmpdir."scripts_out/list";
   open(LIST, ">".$list) || die "counldn't write to $list because $!";
   my $nmany=0;
   my $run_args = "--species $species --bed $bed --index $index --s_list $s_list --k $k_size --cutoff $cutoff --tmpdir $tmpdir";
   if ($resources)
   {
      $run_args .= " --resources \"".$resources."\"";
   }

   my $r=0;
   until ($r >= $nlines)
   {
      $nmany++;
      $r +=$rangesize;
      my $b = $r - $rangesize+1;
      if ($r > $nlines)
      {
	 $r = $nlines;
      }

      print ARGS $this_script." ".$run_args." --range ".$b."-".$r."\n";
      my $range_local = $b."-".$r;
      print LIST $range_local."\n";
   }

   my $runscript = $tmpdir."cluster_scripts/run_kmers.sh";
   my $stdout = $tmpdir."scripts_out/$job_name.out";
   my $stderr = $tmpdir."scripts_out/$job_name.err";
   my $rsc;
   if ($resources)
   {
      $rsc = join(" ", $resources, "-o", $stdout, "-e", $stderr);
   }else
   {
      $rsc = join(" ", "-o", $stdout, "-e", $stderr);
   }


   my $first_job_id = &submit_job($argfile, $runscript, $job_name, $job_type, $chunk."|".$nmany, $rsc);

   my $runscript2 = $tmpdir."cluster_scripts/consolidate.sh";
   my $run_args2 = "--outfile $outfile --tmpdir $tmpdir --finish --bed $bed --s_list $s_list --cutoff $cutoff --k $k_size";
   my $stdout2 = $tmpdir."scripts_out/$job_name.out2";
   my $stderr2 = $tmpdir."scripts_out/$job_name.err2";
   #wait for jobs to finish then submit the consolidate script
   if ($resources)
   {
      $rsc = join(" ", $resources, "-o", $stdout2, "-e", $stderr2);
   }else
   {
      $rsc = join(" ", "-o", $stdout2, "-e", $stderr2);
   }
   if ($job_type eq "PBS")
   {
      my @x = split (/\./, $first_job_id);
      my $job_int_id = $x[0];
      $rsc = join(" ", $rsc, "-W", "depend=afterokarray:".$job_int_id."[]"); #-W..blah... means wait for first_job_id to finish
   }elsif ($job_type eq "SGE")
   {
      $rsc = join(" ", $rsc, "-hold_jid", $job_name);
   }
   my @cmd;
   push @cmd, $this_script." ".$run_args2;
   push @cmd, "perl /".$base."/lovci/projects/conservation/compute_sig.pl $kmerlistfile";
   &submit_job(\@cmd, $runscript2, $job_name, $job_type, "0", $rsc);
   exit;
}


#script starts here 
my %C_OR_NO;


my ($range_start, $range_stop) =split (/-/, $range);
my $head = $range_stop;
my $tail = $range_stop -$range_start + 1;

open(MA, ">".$tmpdir."trash/multiz.".$range); # store mutiz alignments here
open(CONS, ">".$tmpdir."trash/cons.".$range); # store mutiz alignments here
open(FA, ">".$tmpdir."trash/fa.".$range); # store mutiz alignments here
open(NONCONS, ">".$tmpdir."trash/noncons.".$range); # store mutiz alignments here
my $n = 0;

if ($bed =~ /gz$/)
{
   open(BED, "gunzip -c $bed | head -n $head | tail -n $tail |") ||die;
}else
{
   open(BED, "head -n $head $bed | tail -n $tail |") ||die;
}

while (<BED>)
{
   chomp;
   my %MA;
   next unless ($_ =~ /^chr/);
   my ($chr, $start, $stop, $name, $score, $signstrand, $thick_start, $thick_stop, $color) = split (/\t/, $_);
   my $region_length = $stop -$start+1;
   next unless ($region_length > $min_length);
   my $strand = 1;
   if ($signstrand eq "-")
   {
      $strand = -1;
   }
   my @stats = &getstats_maf($chr, $start, $stop, $strand, $species, $pivot, $s_list);
   foreach (@stats)
   {
      my ($s, $seq) = split (/\t/, $_);
      $MA{$s} = $seq;
   }
   my ($cts, $c_seq, $nc_seq) = &count_cons_non(\%MA, $s_list, $cutoff, $k_size);
   my %KMERS = %$cts;
   foreach my $kmer (keys %{$KMERS{conserved}})
   {
      $C_OR_NO{$kmer}{conserved}++;
   }
   foreach my $kmer (keys %{$KMERS{notconserved}})
   {
      $C_OR_NO{$kmer}{notconserved}++;
   }
   $n++;
   if ($n == 100 )
   {
#      last;
   }
   my $id = ">".$chr.":".$start."-".$stop."|".$signstrand."|".$name;


   print CONS $id."\n".$c_seq."\n";
   print NONCONS $id."\n".$nc_seq."\n";

   print MA $id."\n";
   print FA ">".$id."\n".$MA{$species}."\n";
   foreach my $spe (@species_list)
   {
      print MA "=".$spe."\t".$MA{$spe}."\n";
   }
}
close(MA);
open(TRASH, ">".$tmpdir."trash/file.".$range);

foreach my $kmer (keys %C_OR_NO)
{
   if (!$C_OR_NO{$kmer}{conserved})
   {
      $C_OR_NO{$kmer}{conserved}=0;
   }
   if (!$C_OR_NO{$kmer}{notconserved})
   {
      $C_OR_NO{$kmer}{notconserved}=0;
   }
   print TRASH $kmer."\t".$C_OR_NO{$kmer}{conserved}."\t".$C_OR_NO{$kmer}{notconserved}."\n";
}
close(TRASH);
exit;

sub count_cons_non{
   my $href = shift;
   my $list = shift;
   my $cutoff = shift;
   my $k = shift;
   my @list = split (/\-/, $list);
   my $s = shift @list; ##reference species
   my %MULTIZ = %$href;
   my %KMERS;
   $KMERS{'conserved'} = ();
   $KMERS{'notconserved'} = ();
   my @kmers;
   my $size = length($MULTIZ{$s});

   my @cons_seq;
   my @nc_seq;

WORD:   for (my $i = 0; $i<($size-$k); $i++)
   {
      my $kmer = substr($MULTIZ{$s}, $i, $k);
      if ($kmer =~ /^\-/)
      {
	 next WORD; #methods to handle gaps are below, but this won't work if it's the first letter
      }
      my $nuc = substr($MULTIZ{$s}, $i, 1);
      if ($kmer =~ /\#/)#skip gaps in the alignment... proabably should be improved.
      {
	 next;
      }
      my $CS = substr($MULTIZ{"++"}, $i, $k);

      unless ($CS =~ /\./) #if it's a perfect alignment, count it and move on.
      {
	 $KMERS{'conserved'}{$kmer}++;

	 push @cons_seq, $nuc;
	 push @nc_seq, "N";

	 next;
      }

      while ($kmer =~ /\-/) #remove gaps in the alignment
      {
	 my $n =()= $kmer =~ /\-/gi;
	 if ($i+$k+$n >= $size)
	 {
	    last WORD;#could not get a k-mer at offset i
	 }
	 $kmer = substr($MULTIZ{$s}, $i, $k+$n);
	 my $newn =()= $kmer =~ /\-/gi;
	 if ($newn == $n){
	    $kmer =~ s/\-//g;
	    last;
	 }
      }
      my $n_same = 0;
      my $nuc_same =0;
      SUBWORD:foreach my $ali (@list)
      {
	 my $t_kmer = substr($MULTIZ{$ali}, $i, $k);
	 my $t_nuc = substr($MULTIZ{$ali}, $i, 1);

	 while ($t_kmer =~ /\-/) #remove gaps in the alignment
	 {
	    my $n =()= $t_kmer =~ /\-/gi;
	    if ($i+$k+$n >= $size)
	    {
	       next SUBWORD;#could not get a k-mer at offset i
	    }
	    $t_kmer = substr($MULTIZ{$ali}, $i, $k+$n);
	    my $newn =()= $t_kmer =~ /\-/gi;
	    if ($newn == $n){
	       $kmer =~ s/\-//g;
	       last;
	    }
	 }
	 my $same = "-";
	 if ($t_kmer eq $kmer)
	 {
	    $same="+";
	    $n_same++;
	 }
	 if ($t_nuc eq $nuc)
	 {
	    $nuc_same++;
	 }

#	 print $same."\t";
#	 foreach (0..$i)
##	 {
#	    print " ";
#	 }
#	 print $kmer."\n";

	 if ($n_same >= $cutoff)
	 {
	    last SUBWORD;
	 }
      }
      if ($n_same >= $cutoff)
      {
	 $KMERS{'conserved'}{$kmer}++
      }else
      {
	 $KMERS{'notconserved'}{$kmer}++
      }
      if ($nuc_same >= $cutoff)
      {
 	 push @cons_seq, $nuc;
 	 push @nc_seq, "N";
      }else
      {
 	 push @cons_seq, "N";
 	 push @nc_seq, $nuc;
      }
   }

 LETTER: for (my $i=($size-$k); $i<$size; $i++) #do the rest
    {
       my $nuc = substr($MULTIZ{$s}, $i, 1);
       if ($nuc =~ /\#|\-/) #skip gaps
       {
 	 next;
       }
       my $n_same = 0;
       if ($nuc eq "N")
       {
 	 push @cons_seq, "N";
 	 push @nc_seq, "N";
 	 next LETTER;
       }
     SUBWORD: foreach my $ali (@list)
       {
 	 my $subnuc = substr($MULTIZ{$ali}, $i, 1);
 	 if ($subnuc eq $nuc)
 	 {
 	    $n_same++;
 	 }
       }
       if ($n_same >= $cutoff)
       {
 	 push @cons_seq, $nuc;
 	 push @nc_seq, "N";
 	 next LETTER;
       }else
       {
 	 push @cons_seq, "N";
 	 push @nc_seq, $nuc;
       }
    }
   my $c_seq = join("", @cons_seq);
   my $nc_seq = join("", @nc_seq);
   return(\%KMERS, $c_seq, $nc_seq);
}

sub count_cons_non_old{
   my $href = shift;
   my $list = shift;
   my $cutoff = shift;
   my $k = shift;
   my @list = split (/\-/, $list);
   my $s = shift @list; ##reference species
   my %MULTIZ = %$href;
   my %KMERS;
   $KMERS{'conserved'} = ();
   $KMERS{'notconserved'} = ();
   my @kmers;
   my $size = length($MULTIZ{$s});
#   print $s."\t".$MULTIZ{$s}."\n";
#   foreach my $ali (@list)
#   {
#      print $ali."\t".$MULTIZ{$ali}."\n";
#   }
#   print "++"."\t".$MULTIZ{"++"}."\n";
WORD:   for (my $i = 0; $i<($size-$k); $i++)
   {
      my $kmer = substr($MULTIZ{$s}, $i, $k);
      die "fix bug: kmers starting with gaps should be handled better";
      if ($kmer =~ /\#/)#skip gaps in the alignment... proabably should be improved.
      {
	 next;
      }
      my $CS = substr($MULTIZ{"++"}, $i, $k);

      unless ($CS =~ /\./) #if it's a perfect alignment, count it and move on.
      {
	 $KMERS{'conserved'}{$kmer}++;
	 next;
      }

      while ($kmer =~ /\-/) #remove gaps in the alignment
      {
	 my $n =()= $kmer =~ /\-/gi;
	 if ($i+$k+$n >= $size)
	 {
	    last WORD;#could not get a k-mer at offset i
	 }
	 $kmer = substr($MULTIZ{$s}, $i, $k+$n);
	 my $newn =()= $kmer =~ /\-/gi;
	 if ($newn == $n){
	    $kmer =~ s/\-//g;
	    last;
	 }
      }
      my $n_same = 0;
      SUBWORD:foreach my $ali (@list)
      {
	 my $t_kmer = substr($MULTIZ{$ali}, $i, $k);

	 while ($t_kmer =~ /\-/) #remove gaps in the alignment
	 {
	    my $n =()= $t_kmer =~ /\-/gi;
	    if ($i+$k+$n >= $size)
	    {
	       next SUBWORD;#could not get a k-mer at offset i
	    }
	    $t_kmer = substr($MULTIZ{$ali}, $i, $k+$n);
	    my $newn =()= $t_kmer =~ /\-/gi;
	    if ($newn == $n){
	       $kmer =~ s/\-//g;
	       last;
	    }
	 }
	 my $same = "-";
	 if ($t_kmer eq $kmer)
	 {
	    $same="+";
	    $n_same++;
	 }

#	 print $same."\t";
#	 foreach (0..$i)
##	 {
#	    print " ";
#	 }
#	 print $kmer."\n";

	 if ($n_same >= $cutoff)
	 {
	    last SUBWORD;
	 }
      }
      if ($n_same >= $cutoff)
      {
	 $KMERS{'conserved'}{$kmer}++
      }else
      {
	 $KMERS{'notconserved'}{$kmer}++
      }
   }
   return(\%KMERS);
}

sub getstats_maf {
   my $chr = shift;
      $chr =~ s/chr//g;
   my $x = shift;
   my $y = shift;
   my $strand = shift;
   my $species = shift;
   my $pivot = shift;
   my $s_list = shift;

   my $input = "$species $pivot $s_list $chr $x $y $strand -1";
   print $input."\n";
   exit;
#   print $input;
#   die;
#   print $input."\n";
   open(CMD,"perl /".$base."/yeolab/Conservation/scripts/consensus_norealign2.pl $input 2>/dev/null |");



#   open(CMD,"perl /nas3/yeolab/Conservation/scripts/consensus_norealign2.pl $input |");
   my @tmp; 
   while (<CMD>)
   {
      chomp;
      push(@tmp,$_);
   }
   close(CMD);

   return @tmp;
}

sub submit_job{
    my $command = shift; # an array ref of  commands to run or if $parallel is not zero then this is a file with commands to run listed on new lines
    my $job_file = shift; # a file to put the generated script in
    my $job_name = shift; # a name for submitted job
    my $job_type = shift; # pbs or sge
    my $parallel = shift; # 0 or chunk|nmany
    my $resources = shift; # additional resources for job submission
    my @resources = split (/\s/, $resources);
    my $length = scalar @resources;
    my $njobs = 1;
#    print STDERR "Attempting to submit $job_file as $job_name ....";
    unless ($length % 2 == 0){
       print $resources."\n";
	die "uneven number of resources, specify <-option>, <resource>, <-option>, <resource>";
    }
    open (OUT, ">".$job_file)||die;

    my $marker;
    my $job_id;
    my $cwd;
    if ($job_type eq "SGE")
    {
       $marker = '$';
       $job_id= '$SGE_TASK_ID';
       $cwd = '#$ -cwd';
    }elsif ($job_type eq "PBS")
    {
       $marker = 'PBS';
       $job_id = '$PBS_ARRAYID';
       $cwd = 'cd $PBS_O_WORKDIR';
    }else
    {
       die;
    }
    print OUT "#!/bin/sh 
#".$marker." -N $job_name
#".$marker." -V
#".$marker." -S /bin/sh\n";

    for (my $i =0; $i < scalar @resources; $i+=2)
    {
       my $resource = $resources[$i];
       my $value = $resources[$i+1];
       print OUT "#".$marker." ".$resource." ".$value."\n";
    }
    if ($parallel ne 0)
    {
       my $argfile = $command;
       my ($chunk, $nmany) = split (/\|/, $parallel);
       if (!$nmany)
       {
	  die;
       }
       my $chunk2=$chunk-1;

       if ($nmany%$chunk == 0)
       {
	  $njobs = $nmany/$chunk
       } else
       {
	  $njobs = int(1+($nmany/$chunk));
       }
       print OUT "#".$marker." -t 1-$njobs
$cwd
#echo \"hostname is:\"
#hostname

index=0
lastindex=0
let \"index=$job_id*$chunk\"
let \"lastindex=\$index-$chunk2\"\n";

       if ($chunk2 > 0)
       {
	  print OUT "for i in `seq \$lastindex \$index`\n";
       } else {
	  print OUT "for i in \$index\n"; # if user chooses 1 for chunk size
       }
       print OUT "do
ARGFILE=\"$argfile\"
line=\$(cat \$ARGFILE | head -n \$i | tail -n 1)
eval \$line
done
";
    }else
    {
print OUT "$cwd
echo \"hostname is:\"
hostname\n\n";
       my @commands =@$command;

       foreach my $cmd (@commands){
	  print OUT $cmd."\n";
       }
    }
    close(OUT);
    my $jobid = `qsub $job_file`;
#    if ($job_type eq "PBS")
#    {
#       my @job_id_list = ();
#       my ($id, $triton, $sdsc, $edu) = split (/\./, $jobid);
#       foreach my $j (1..$njobs)
#       {
#	  push @job_id_list, $id."-".$j
#       }
#       $jobid = join(",", @job_id_list);
#    }
    print $jobid;
    return($jobid);
#    print STDERR "Successfully submitted $job_file as $job_name\n";
}

