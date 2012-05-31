use warnings;
use strict;

my $n = 0;
my ($d,$target,$lengthtarget,$miRNA,$lengthmiR,$mfe,$pvalue,$position, $matchmiRNA, $matchTarget);
while (<STDIN>)
{
   chomp;
   if (/^target:/)
   {
      ($d,$target) = split(/\:\s/,$_); 
      ($d,$lengthtarget) = split(/\:\s/,<STDIN>);chomp $lengthtarget;
      ($d,$miRNA) = split(/\:\s/,<STDIN>); chomp $miRNA;
      ($d,$lengthmiR) = split(/\:\s/,<STDIN>);<STDIN>;chomp $lengthmiR;
      ($d,$mfe) = split(/\s/,<STDIN>);
      ($d,$pvalue) = split(/\:\s/,<STDIN>);chomp $pvalue; <STDIN>;
      ($d,$position) = split(/\s+/,<STDIN>);chomp $position;

      my $target_unmatched = <STDIN>; chomp $target_unmatched;
      $target_unmatched =~ s/(target 5\' )|( 3')//g;
      $target_unmatched = lc($target_unmatched);
      my $target_matched = <STDIN>; chomp $target_matched;
      $target_matched =~ s/(^          )|(   $)//g;

      my $targetseq = "";
      $target_matched =~ s/\s/-/g;
      $target_unmatched =~ s/\s/-/g;
      my @matched = split(//, $target_matched);
      
      my @unmatched = split(//, $target_unmatched);

      for (my $i =0; $i< scalar(@unmatched); $i++)
      {
	 if ($matched[$i] eq "-")
	 {
	    $matched[$i] = $unmatched[$i];
	 }
      }
      $targetseq = join("", @matched);

      my $miRNA_matched = <STDIN>; chomp $miRNA_matched;
      $miRNA_matched =~ s/(^          )|(   $)//g;
      my $miRNA_unmatched = <STDIN>; chomp $miRNA_unmatched;
      $miRNA_unmatched =~ s/(miRNA  3\' )|( 5')//g;
      $miRNA_unmatched = lc($miRNA_unmatched);
      $miRNA_matched =~ s/\s/-/g;
      $miRNA_unmatched =~ s/\s/-/g;
      my $miRNAseq = "";

      @matched = split(//, $miRNA_matched);
      @unmatched = split(//, $miRNA_unmatched);

      for (my $i =0; $i< scalar(@unmatched); $i++)
      {
	 if ($matched[$i] eq "-")
	 {
	    $matched[$i] = $unmatched[$i];
	 }
      }
      $miRNAseq = reverse(join("", @matched));
      print join("\t", $target, $lengthtarget, $position, $miRNA, $lengthmiR, $mfe, $pvalue, $targetseq, $miRNAseq)."\n";
   }
}
