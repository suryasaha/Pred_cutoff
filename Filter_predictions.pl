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
use POSIX;

if(!($^O=~ /linux/)){ print STDERR "Script has not been tested on non-linux operating systems. Exiting.."; exit 1;}

=head1 NAME

Filter_predictions.pl - Filter out the predictions that are located within a user defined range 
                        of genomic features such as genes.

=head1 SYNOPSIS

  % Filter_predictions.pl --predictions pred.gff --features features.gff --range range
  
=head1 DESCRIPTION

This script scans the list of predictions and simply pulls out those within a user defined range. The 
input files need to be in GFF3 format such as those produced by Pred_cutoff.pl.

Please cite the following publication if you use this pipeline
Bound to succeed: Transcription factor binding site prediction and its contribution to understanding virulence and environmental adaptation in bacterial plant pathogens
Saha, Surya and Lindeberg, Magdalen (2013) Molecular Plant-Microbe Interactions. doi: 10.1094/MPMI-04-13-0090-CR

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options. e.g. -f instead of --file. All options
are mandatory (see below).  

   --predictions  <GFF>       Predictions in GFF3 format (required)
   --features     <GFF>       Genomic features such as genes in GFF3 format (required)
   --range        <range>     Window for qualifying predictions (required)
   --verbose                  Print progress messages
   
=head1 AUTHOR

Surya Saha, ss2489 near cornell.edu , \@SahaSurya

=cut

my ($i,$preds,$features,$verbose,$range);

GetOptions (
	'predictions=s' => \$preds, 
	'features=s' => \$features,
	'range=i' => \$range,
	verbose => \$verbose) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($preds) or (system('pod2text',$0), exit 1);
if (!(-e $preds)){print STDERR "$preds not found: $!\n"; exit 1;}
if($verbose){print STDERR "Setting prediction file to $preds\n";}
defined($features) or (system('pod2text',$0), exit 1);
if (!(-e $features)){print STDERR "$features not found: $!\n"; exit 1;}
if($verbose){print STDERR "Setting feature file to $features\n";}
defined($range) or (system('pod2text',$0), exit 1);
if($verbose){print STDERR "Setting range to $range\n";}
print STDERR "\n";

my ($rec,$ctr,@table,@sel,$j,$k,$l,@temp);
unless(open(INPREDS,$preds)){print "not able to open $preds\n\n";exit 1;}
unless(open(INFTRS,$features)){print "not able to open $features\n\n";exit 1;}
unless(open(OUT,">filtered.$preds")){print "not able to open filtered.$preds\n";exit 1;}

# read in files
while($rec=<INPREDS>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec); $temp[2]='prediction';
	push @table, [@temp];
}
close (INPREDS); if($verbose){print STDERR "Finished reading in predictions\n";}
while($rec=<INFTRS>){
	if($rec =~ /#/){next;}
	@temp=split("\t",$rec); $temp[2]='feature';
	push @table, [@temp];
}
close (INFTRS); if($verbose){print STDERR "Finished reading in features\n";}

# write out filtered predictions
# sorting the records
@temp = sort {$a->[3] <=> $b->[3]} @table; @table=@temp;
# iterating over the combined table to find qualifying promoters
@temp=();
for ($ctr=0;$ctr<@table;$ctr++){
	if ($table[$ctr][2] eq 'prediction' && $table[$ctr][6] eq '+'){
		# if pos strand, look at genes on right
		for($i=$ctr+1;$i<@table;$i++){# if gene, same strand and within range
			if (($table[$i][2] eq 'feature') && ($table[$i][6] eq '+') && ($table[$i][3]<=($table[$ctr][4]+$range)) && ($table[$i][3]>$table[$ctr][4])){#no overlap
				$temp[0]=$ctr; $temp[1]=$table[$i][8]; $temp[2]=$table[$i][3]; $temp[3]=$table[$i][4];
				push @sel,[@temp]; @temp=(); last;
			}
			elsif(($table[$i][3]>($table[$ctr][4]+$range))){ last;}
		}
	}
	elsif ($table[$ctr][2] eq 'prediction' && $table[$ctr][6] eq '-'){
		# if neg strand, look at genes on left
		for($i=$ctr-1;$i>=0;$i--){# if gene, same strand and within range
			if (($table[$i][2] eq 'feature') && ($table[$i][6] eq '-') && ($table[$i][4]>=($table[$ctr][3]-$range)) && ($table[$i][4]<$table[$ctr][3])){# no overlap
				$temp[0]=$ctr; $temp[1]=$table[$i][8]; $temp[2]=$table[$i][3]; $temp[3]=$table[$i][4];
				push @sel,[@temp]; @temp=(); last;
			}
			elsif(($table[$i][4]<($table[$ctr][3]-$range))){ last;}
		}
	}
}

# print out .gff file if qualified promoters with related gene info
if($verbose){print STDERR "Printing out predictions within $range of features\n";}
foreach $i (@sel){
	print OUT "$table[$i->[0]][0]\t$table[$i->[0]][1]\t$table[$i->[0]][2]\t$table[$i->[0]][3]\t$table[$i->[0]][4]\t";
	print OUT "$table[$i->[0]][5]\t$table[$i->[0]][6]\t$table[$i->[0]][7]\t";
	chomp($table[$i->[0]][8]); chomp($i->[1]);
	print OUT $table[$i->[0]][8],"$i->[1]\; matchstart=$i->[2]\; matchend=$i->[3]\;\n";
}
close (OUT);

if($verbose){
	my($user_t,$system_t,$cuser_t,$csystem_t);	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print STDERR "\nSystem time for process: $system_t\n"; print STDERR "User time for process: $user_t\n";
}
