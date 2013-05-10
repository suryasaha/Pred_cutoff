#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha  
# http://citrusgreening.org/pred_cutoff

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
# Development of this software was done with support from the Citrus Research 
# and Development Foundation and it cannot be copyrighted. This software is freely 
# available to the public for use without restriction on its use or reproduction.
#
# Although all reasonable efforts have been taken to ensure the accuracy
# and reliability of the software and data, the Citrus Research and Development 
# Foundation does not and cannot warrant the performance or results that may
# be obtained by using this software or data. The Citrus Research and Development
# Foundation disclaims all warranties, express or implied, including warranties
# of performance, merchantability or fitness for any particular purpose.
#
# Please cite the author in any work or product based on this material.
#
# ===========================================================================

use strict;
use warnings;
use Getopt::Long;
use POSIX;

if(!($^O=~ /linux/)){ print STDERR "Script has not been tested on non-linux operating systems. Exiting.."; exit 1;}


=head1 NAME

Build_model.pl - Build a model (HMM or PWM) using a ClustalW alignment of DNA motifs

=head1 SYNOPSIS

  % Build_model.pl --file seq.fna --alignment motifs.aln --motif_length length --engine HMMER/patser 
  
=head1 DESCRIPTION

This script uses a alignment produced in Clustal format to create a model to be used by Search_pred.pl and 
Pred_cutoff.pl scripts. It uses hmmbuild from the HMMER suite for a hidden markov model and make-matrix from
the Consensus suite for a position weight matrix model. 

=head2 PATSER REQUIREMENTS

Comment out non-sequence lines with '#' in the alignment file other make-matrix will produce a unerecognized
character error. Please make sure that you do not have gaps in your alignment if using patser. If you do have 
gaps, we recommend that you use multiple position weight matrices, one for each ungapped alignment 
and remove redundant predictions at the end. 

=head2 NOTES

This script has been tested on Linux and requires make-matrix (ftp://www.genetics.wustl.edu/pub/stormo/Consensus)
and HMMER 2.3.2 (ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz) to be installed and 
accessible. We do not recommend using HMMER 3.0 since it is not optimized for DNA/DNA comparisons. 

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. All options
are mandatory (see below).  

   --file         <file>      Fasta file whole genome DNA sequence (required)
   --alignment    <ALN>       A ClustalW formatted alignment file (required)
   --motif_length <length>    Length of target DNA motif (required)
   --engine       <engine>    Type of engine to be used to create the model (HMMER or patser) (required)
   --verbose                  Print progress messages
   
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut


my ($i,$file,$aln,$engine,$mlen,$verbose);

GetOptions (
	'file=s' => \$file,
	'alignment=s' => \$aln,
	'motif_length=i' => \$mlen, 
	'engine=s' => \$engine,
	verbose => \$verbose) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($file) or (system('pod2text',$0), exit 1);
if (!(-e $file)){print STDERR "$file not found: $!\n"; exit 1;}
defined($aln) or (system('pod2text',$0), exit 1);
if (!(-e $aln)){print STDERR "$aln not found: $!\n"; exit 1;}
defined($mlen) or (system('pod2text',$0), exit 1);
defined($engine) or (system('pod2text',$0), exit 1);
if(($engine ne 'HMMER') && ($engine ne 'patser')){print STDERR "Incorrect engine name supplied. $engine should be HMMER or patser.\n"; exit 1;}
if($verbose){print STDERR "Setting engine type to $engine\n";}
print STDERR "\n";

# Prep work
my ($rec,$ctr,$j,$k,@temp,$seq,%counts);
unless(open(IN,$file)){print "not able to open $file\n\n";exit 1;}
$seq='';
while ($rec=<IN>){
	if ($rec=~ /^>/){ next;}
	else{ chomp $rec; $seq=$seq.$rec;}
}
close(IN);

#create model for engine
if($engine eq 'patser'){
	#create alphabet file
	if($verbose){print STDERR "Creating alphabet file for patser..\n";}
	$i= length $seq;
	$counts{'A'}=($seq =~ tr/A//);
	if($counts{'A'}==0){ $counts{'A'}=($seq =~ tr/a//);}
	$counts{'T'}=($seq =~ tr/T/X/);
	if($counts{'T'}==0){ $counts{'T'}=($seq =~ tr/t//);}
	$counts{'G'}=($seq =~ tr/G/Y/);
	if($counts{'G'}==0){ $counts{'G'}=($seq =~ tr/g//);}
	$counts{'C'}=($seq =~ tr/C/Z/);
	if($counts{'C'}==0){ $counts{'C'}=($seq =~ tr/c//);}
	if (-e "$file.alphabet"){ unlink "$file.alphabet" or warn "cannot delete old alphabet file: $!\n";}
	unless(open(OUT,">$file.alphabet")){print "not able to open $file.alphabet\n";exit 1;}
	print OUT "\#\n\# Alphabet file for patser\n\#\n\n";
	$j=$counts{'A'}+$counts{'T'}; $j=sprintf("%.1f", $j/$i); print OUT "a:t $j\n";
	$j=$counts{'G'}+$counts{'C'}; $j=sprintf("%.1f", $j/$i); print OUT "g:c $j\n";
	close (OUT);
	%counts=();
	
	if($verbose){print STDERR "Creating PWM for patser using $aln ..";}
	$k=system "cat $aln | make-matrix-v2 -a $file.alphabet > $aln.pwm";
	if($k!=0){ print STDERR "make-matrix execution failed: $?\n"; exit 1;}
	unlink "$file.alphabet";
}
elsif($engine eq 'HMMER'){
	#create null file
	if($verbose){print STDERR "Creating null file for HMMER..\n";}
	$i= length $seq;
	$counts{'A'}=($seq =~ tr/A//);
	if($counts{'A'}==0){ $counts{'A'}=($seq =~ tr/a//);}
	$counts{'T'}=($seq =~ tr/T/X/);
	if($counts{'T'}==0){ $counts{'T'}=($seq =~ tr/t//);}
	$counts{'G'}=($seq =~ tr/G/Y/);
	if($counts{'G'}==0){ $counts{'G'}=($seq =~ tr/g//);}
	$counts{'C'}=($seq =~ tr/C/Z/);
	if($counts{'C'}==0){ $counts{'C'}=($seq =~ tr/c//);}
	if (-e "$file.null"){ unlink "$file.null" or warn "cannot delete old null file: $!\n";}
	unless(open(OUT,">$file.null")){print "not able to open $file.null\n";exit 1;}
	print OUT "\#\n\# Null model for HMMER\n\#\n\nNucleic\n\n";
	$j=sprintf("%.6f", $counts{'A'}/$i); print OUT "$j \# A\n";
	$j=sprintf("%.6f", $counts{'C'}/$i); print OUT "$j \# C\n";
	$j=sprintf("%.6f", $counts{'G'}/$i); print OUT "$j \# G\n";
	$j=sprintf("%.6f", $counts{'T'}/$i); print OUT "$j \# T\n\n";
	$j=sprintf("%.6f", $mlen/($mlen+1)); print OUT "$j \# p1 = $mlen/($mlen+1)\n";
	close (OUT);
	%counts=();
	
	#run hmmbuild and calibrate
	if($verbose){print STDERR "Creating HMM for HMMER using $aln ..";}
	$k=system "hmmbuild --archpri 0.95 --nucleic -F --null $file.null $file.$aln.hmm $aln 1>$file.hmmbuild.out";
	if($k!=0){ print STDERR "hmmbuild execution failed: $?\nSee $file.hmmbuild.out for errors\n"; exit 1;}
	unlink "$file.null"; unlink "$file.hmmbuild.out";
	$k=system "hmmcalibrate --num 50000 $file.$aln.hmm 1>$file.hmmcalib.out";
	if($k!=0){ print STDERR "hmmcalibrate execution failed: $?\nSee $file.hmmcalib.out for errors\n"; exit 1;}
	unlink "$file.hmmcalib.out";
}

if($verbose){
	my($user_t,$system_t,$cuser_t,$csystem_t);	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print STDERR "\n\nSystem time for process: $system_t\n"; print STDERR "User time for process: $user_t\n";
}