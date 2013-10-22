#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha  
# https://github.com/suryasaha/Pred_cutoff

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
if(!($^O=~ /linux/)){ print STDERR "Script has not been tested on non-linux operating systems. Exiting.."; exit 1;}

=head1 NAME

Create_art_genomes.pl - Create optimally ordered artificial genomes using the seqpp toolkit 

=head1 SYNOPSIS

  % Create_art_genomes.pl --file seq.fna --copies copies
  
=head1 DESCRIPTION

This script reads in a Fasta formatted DNA sequence representing the whole genome sequence of 
a organism. It requires the seqpp toolkit (http://stat.genopole.cnrs.fr/seqpp/, ver 4.2.0+) to 
be installed. It uses the estim_m tool to compute the optimal order for the markov chain to 
simulate the genome. The simul_m tool is then used to create artificial genomes that are placed 
inside the genomes directory. 

=head2 NOTES

Please make sure that the estim_m and simul_m tools are accessible from the command-line 
in the directory where you are running this script. The Fasta file should have a single DNA 
sequence. The artificial genomes created will be stored in genomes directory. This script has 
been tested on Linux.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some
options are mandatory (see below).  

   --file     <file>      Fasta file whole genome DNA sequence (required)
   --copies   <copies>    Number of artificial genomes to create. Default is 600
   --verbose              Print progress messages
   
=head1 AUTHOR

Surya Saha, ss2489 near cornell.edu , \@SahaSurya

=cut

my ($i,$file,$copies,$verbose);
GetOptions ('file=s' => \$file, 
	'copies:s' => \$copies,
	verbose => \$verbose) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($file) or (system('pod2text',$0), exit 1);
if (!(-e $file)){print STDERR "$file not found: $!\n"; exit 1;}
if(!defined($copies)){print STDERR "Creating 600 artificial genomes..\n";}
$copies ||= 600;
$i=localtime();
if (-e 'genomes'){rename 'genomes',"genomes_Create_art_genomes_$i"; unlink glob "genomes/* genomes/.*";rmdir ("genomes");}
mkdir ("genomes", 0755) or warn "Cannot make genomes directory: $!\n";

my ($rec,@bic,$ctr,$order,$j,$k,$seq,@temp);

# find optimal order
if (-e 'temp.bic'){ unlink 'temp.bic' or warn "cannot delete old temporary files: $!\n"; exit 1;}
for $i (1..10){
	if($verbose){print STDERR "Computing BIC for degree $i..\n";}
	else{print STDERR '.';}
	$j=`estim_m $file -d $i --Bic`;
	@temp = split(' ',$j);
	chomp $temp[$#temp]; $bic[$i-1]=$temp[$#temp];
	@temp=();
}

$j=$bic[0];
foreach $i (1..9){
	if ($bic[$i]>$bic[$i-1]){ $order=$i; last;}
}

# create model for optimal order
if (-e 'temp.model'){ unlink 'temp.model' or warn "cannot delete old temporary files: $!\n";}
if($verbose){print STDERR "\nCreating model file of order $order..\n\n";}
else{print STDERR '.';}
$k=system "estim_m $file -d $order -o temp.model";
if($k!=0){ print STDERR "estim_m execution failed: $?\n"; exit 1;}

# create artificial genomes
unless(open(INFNA,$file)){print "not able to open $file\n\n";exit 1;}
$seq='';
while ($rec=<INFNA>){
	if ($rec=~ /^>/){ next;}
	else{ chomp $rec; $seq=$seq.$rec;}
}
close(INFNA);

$i= length $seq;
for $j (1..$copies){
	if($verbose){print STDERR "Creating artificial genome $j..\n";}
	else{print STDERR '.';}
	$k=system "simul_m -l $i -m temp.model -o genomes/$j.$file";
	if($k!=0){ print STDERR "simul_m execution failed: $?\n"; exit 1;}
}
unlink 'temp.model' or warn "cannot delete temporary files: $!\n";
if(!$verbose){print STDERR "\n";}

if($verbose){
	my($user_t,$system_t,$cuser_t,$csystem_t);	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print STDERR "\n\nSystem time for process: $system_t\n"; print STDERR "User time for process: $user_t\n";
}