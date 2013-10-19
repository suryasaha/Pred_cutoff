#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha  
# https://github.com/suryasaha/Pred_cutoff

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
# This software is freely available to the public for use and the 
# authors have not placed any restriction on its use or reproduction.
# Although all reasonable efforts have been taken to ensure the accuracy
# and reliability of the software and data, the authors do not and 
# cannot warrant the performance or results that may be obtained by
# using this software or data. The authors disclaim all warranties, 
# express or implied, including warranties of performance, merchantability 
# or fitness for any particular purpose.
#
# Please cite the author in any work or product based on this material.
#
# ===========================================================================

use strict;
use warnings;
use Getopt::Long;
use POSIX;
eval {require Bio::SearchIO::hmmer};
if ($@){print STDERR "Cannot find Bio::SearchIO::hmmer\n"; exit 1;}
use Bio::SearchIO;

if(!($^O=~ /linux/)){ print STDERR "Script has not been tested on non-linux operating systems. Exiting.."; exit 1;}

=head1 NAME

Search_pred.pl - Run a model (HMM or PWM) on whole genome DNA sequence using a pattern matching engine and 
                 report the predictions in GFF format. 

=head1 SYNOPSIS

  % Search_pred.pl --file seq.fna --engine HMMER/patser --model HMM/PWM
  
=head1 DESCRIPTION

This script scans a whole genome DNA sequence for a model. It uses HMMER for a hidden markov model and patser 
for a position weight matrix model. A hidden markov model can be created using the hmmbuild tool in the HMMER 
2.3.2 suite. A position weight matrix can be built using a number of tools such as the make_matrix tool. The 
predictions are reported in GFF format. This GFF file is required by the Pred_cutoff.pl script to add in E-value 
estimation.

=head2 NOTES

This script has been tested on Linux and requires patser (http://ural.wustl.edu/patser-v3e.1.tar.gz) and
HMMER 2.3.2 (ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz) to be installed and 
accessible. We do not recommend using HMMER 3.0 since it is not optimized for DNA/DNA comparisons. You 
need to have the BioPerl library installed and accessible (http://www.bioperl.org/wiki/Installing_BioPerl).

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options. e.g. -f instead of --file. Some options
are mandatory (see below).  

   --file         <file>      Fasta file whole genome DNA sequence (required)
   --engine       <engine>    Tool used to create the report (HMMER or patser) (required)
   --model        <HMM/PWM>   Model created from motifs. Hidden Markov Model (.hmm) for
                              HMMER or Position Weight Matrix(.pwm) for patser (required) 
   --domE         <x>         Domain E-value cutoff for hmmsearch on artificial genomes. Default 
                              is 0.3. All predictions with lower E-values will be recorded.
   --score        <x>         Score cutoff for patser on artificial genomes. Default is 1. All 
                              predictions with higher scores will be recorded.
   --verbose                  Print progress messages
   
=head1 AUTHOR

Surya Saha, ss2489 near cornell.edu , \@SahaSurya

=cut

my ($i,$file,$engine,$model,$verbose,$domE,$pscore);

GetOptions (
	'file=s' => \$file, 
	'engine=s' => \$engine,
	'model=s' => \$model,
	'domE:f' => \$domE,
	'score:i' => \$pscore,
	verbose => \$verbose) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($file) or (system('pod2text',$0), exit 1);
if (!(-e $file)){print STDERR "$file not found: $!\n"; exit 1;}
defined($engine) or (system('pod2text',$0), exit 1);
if(($engine ne 'HMMER') && ($engine ne 'patser')){print STDERR "Incorrect engine name supplied. $engine should be HMMER or patser.\n"; exit 1;}
if($verbose){print STDERR "Setting engine type to $engine\n";}
defined($model) or (system('pod2text',$0), exit 1);
if (!(-e $model)){print STDERR "$model not found: $!\nPlease supply Hidden Markov Model (.hmm) for HMMER or Position Weight Matrix(.pwm) for patser on the command line\n"; exit 1;}
if($engine eq 'HMMER'){ defined($domE) or ($domE=0.3); if($verbose){print STDERR "Setting hmmsearch domain E value cutoff to $domE\n\n";}}
elsif($engine eq 'patser'){ defined($pscore) or ($pscore=1); if($verbose){print STDERR "Setting patser score cutoff to $pscore\n";}}
if($verbose){print STDERR "Setting model name to $model\n";}
print STDERR "\n";

# Supporting functions
# get the complement
sub comp{
	my $DNA;
	$DNA=$_[0];	$DNA=~ s/\s*//g; $DNA=~ tr/ACGTacgt/TGCAtgca/;
	$DNA=reverse($DNA);	return $DNA;
}

# Prep work
my ($rec,$ctr,$j,$k,$l,@temp,$seq,%counts,$name);
unless(open(IN,$file)){print "not able to open $file\n\n";exit 1;}
$seq='';
while ($rec=<IN>){
	if ($rec=~ /^>/){ 
		chomp $rec; $rec=~ s/^>//; 
		@temp=split(' ',$rec);
		$name=$temp[0]; @temp=();
	}
	else{ chomp $rec;	$seq=$seq.$rec;}
}
close(IN);

#run engine on genome
if($engine eq 'patser'){
	my ($mlen);
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
	$j=$counts{'A'}+$counts{'T'}; $j=sprintf("%.1f", $j/$i); print OUT "a:t $j\n";
	$j=$counts{'G'}+$counts{'C'}; $j=sprintf("%.1f", $j/$i); print OUT "g:c $j\n";
	close (OUT);
	%counts=();
	
	#create patser compatible artificial genome files, run patser on them and get counts from reports
	unless(open(IN,$file)){print "not able to open $file\n\n";exit 1;}
	unless(open(OUT,">$file.temp")){print "not able to open $file.temp\n\n";exit 1;}
	while($rec=<IN>){
		if($rec =~ /^>/){ print OUT "$name \\\n";}
		else {print OUT $rec;}
	}
	print OUT "\\";
	close (OUT); close (IN);
	if($verbose){print STDERR "Running patser on $file ..";}
	$k=system "patser-v3e -m $model -a $file.alphabet -f $file.temp -c -ls $pscore > $file.patser.out";
	if($k!=0){ print STDERR "patser execution failed: $?\n"; exit 1;}
	unlink "$file.temp"; unlink "$file.alphabet";
	if($verbose){print STDERR " getting results..";}
	unless(open(IN,"$file.patser.out")){print "not able to open $file.patser.out\n\n";exit 1;}
	unless(open(OUT,">$file.patser.gff")){print "not able to open $file.patser.gff\n\n";exit 1;}
	print OUT "##gff-version 3\n"; $rec=<IN>;
	while (!($rec =~ /^average score above/)){
		$rec=<IN>;
		if ($rec =~ /^width of the alignment matrix/){
				$rec =~ s/^width of the alignment matrix: //;
				chomp $rec;	$mlen=$rec;
		}
	}$rec=<IN>;
	while($rec=<IN>){
		$rec=~ s/^ *//;
		if($rec =~ /^$name  /){ 
			chomp $rec;	@temp = split(' ',$rec);
			print OUT "$temp[0]\tpatser\tprediction\t";
			if($temp[2] =~ /C$/){ $temp[2] =~ s/C$//; print OUT "$temp[2]\t",$temp[2]+$mlen,"\t$temp[4]\t-\t.\t";}
			else{ print OUT "$temp[2]\t",$temp[2]+$mlen,"\t$temp[4]\t+\t.\t";}
			print OUT "note \"From PWM $model\"\; ln(p-value) $temp[6]\; score_cutoff $pscore\; evidence not_experimental\;\n";
			@temp=();
		}
	}
	close (IN); close (OUT);
	unlink "$file.patser.out";
} 
elsif($engine eq 'HMMER'){
	#run hmmsearch on genome
	my($result,$hit,$hsp);
	if($verbose){print STDERR "Running hmmsearch on $file ..";}
	#pos strand
	$k=system "hmmsearch --domE $domE $model $file > $file.hmmsearch.pos.out";
	if($k!=0){ print STDERR "hmmsearch execution failed: $?\n"; exit 1;}
	#comp strand
	unless(open(OUT,">comp.$file")){print "not able to open comp.$file\n";exit 1;}
	print OUT ">comp_seq\n"; print OUT &comp($seq);	close(OUT);
	$k=system "hmmsearch --domE $domE $model comp.$file > $file.hmmsearch.comp.out";
	if($k!=0){ print STDERR "hmmsearch execution failed: $?\n"; exit 1;}
	unlink "comp.$file";
	if($verbose){print STDERR " getting results..";}
	unless(open(OUT,">$file.hmmsearch.gff")){print "not able to open $file.hmmsearch.gff\n\n";exit 1;}
	print OUT "##gff-version 3\n";
	
	#read in hmmsearch report for pos
	$j=Bio::SearchIO->new(-format => 'hmmer',-file => "$file.hmmsearch.pos.out");
	while($result=$j->next_result()){# class of $result: Bio::Search::Result::HMMERResult
		while ($hit = $result->next_hit()){# class of $hit: Bio::Search::Hit::HMMERHit
			while ($hsp = $hit->next_hsp()){# class of $hsp: Bio::Search::HSP::HMMERHSP
				print OUT "$name\thmmer\tprediction\t",$hsp->start('subject'),"\t",$hsp->end('subject'),"\t";
				#print OUT $hsp->score(),"\t+\t.\tnote \"From model ",$result->hmm_name,"\"\; E-value ",$hsp->evalue(),"\;";
				print OUT $hsp->score(),"\t+\t.\tnote \"From model ",$result->hmm_name,"\"\;";
				print OUT " evidence not experimental\;\n";
	   		}
	   	}
	}
	#read in hmmsearch report for comp
	$j=Bio::SearchIO->new(-format => 'hmmer',-file => "$file.hmmsearch.comp.out");
	while($result=$j->next_result()){# class of $result: Bio::Search::Result::HMMERResult
		while ($hit = $result->next_hit()){# class of $hit: Bio::Search::Hit::HMMERHit
			while ($hsp = $hit->next_hsp()){# class of $hsp: Bio::Search::HSP::HMMERHSP
				#converting coordinates to pos strand space
				print OUT "$name\thmmer\tprediction\t",length($seq) - $hsp->end('subject') + 1,"\t",
				 length($seq) - $hsp->start('subject') + 1,"\t";
				print OUT $hsp->score(),"\t-\t.\tnote \"From model ",$result->hmm_name,"\"\;";
				print OUT " evidence not experimental\;\n";
	   		}
	   	}
	}
	close(OUT);
	unlink "$file.hmmsearch.pos.out"; unlink "$file.hmmsearch.comp.out";
}

if($verbose){
	my($user_t,$system_t,$cuser_t,$csystem_t);	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print STDERR "\n\nSystem time for process: $system_t\n"; print STDERR "User time for process: $user_t\n";
}
