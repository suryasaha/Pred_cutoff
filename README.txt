
				PRED_CUTOFF README
			http://citrusgreening.org/pred_cutoff




STEP 1
Build_model.pl - Build a model (HMM or PWM) using a ClustalW alignment of DNA motifs
  Build_model.pl --file seq.fna --alignment motifs.aln --motif_length length --engine HMMER/patser

STEP 2
Search_pred.pl - Run a model (HMM or PWM) on whole genome DNA sequence using a pattern matching engine and report the predictions in GFF format. 
  Search_pred.pl --file seq.fna --engine HMMER/patser --model HMM/PWM

STEP 3
Create_art_genomes.pl - Create optimally ordered artificial genomes using the seqpp toolkit 
  Create_art_genomes.pl --file seq.fna --copies copies

STEP 4
Pred_cutoff.eval.pl - Add E values to HMMER/patser report (GFF) and prints out a new GFF file. 
  Pred_cutoff.eval.pl --file seq.fna --copies copies --engine HMMER/patser --report report.gff

STEP 5 (optional)
Filter_predictions.pl - Filter out the predictions that are located within a user defined range of genomic features such as genes.
  Filter_predictions.pl --predictions pred.gff --features features.gff --range range
  




BUILD_MODEL.PL

Build_model.pl - Build a model (HMM or PWM) using a ClustalW alignment of DNA motifs

SYNOPSIS
Build_model.pl --file seq.fna --alignment motifs.aln --motif_length length --engine HMMER/patser 
  
DESCRIPTION
This script uses a alignment produced in Clustal format to create a model to be used by Search_pred.pl and Pred_cutoff.pl scripts. It uses hmmbuild from the HMMER suite for a hidden markov model and make-matrix from the Consensus suite for a position weight matrix model. 

PATSER REQUIREMENTS
Comment out non-sequence lines with '#' in the alignment file other make-matrix will produce a unerecognized character error. Please make sure that you do not have gaps in your alignment if using patser. If you do have gaps, we recommend that you use multiple different position weight matrices, one for each ungapped alignment and remove redundant predictions at the end. 

NOTES
This script has been tested on Linux and requires make-matrix (ftp://www.genetics.wustl.edu/pub/stormo/Consensus) and HMMER 2.3.2 (ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz) to be installed and accessible. We do not recommend using HMMER 3.0 since it is not optimized for DNA/DNA comparisons. 

COMMAND-LINE OPTIONS
Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. All options are mandatory (see below).  
   --file         <file>      Fasta file whole genome DNA sequence (required)
   --alignment    <ALN>       A ClustalW formatted alignment file (required)
   --motif_length <length>    Length of target DNA motif (required)
   --engine       <engine>    Type of engine to be used to create the model (HMMER or patser) (required)
   --verbose                  Print progress messages



SEARCH_PRED.PL

Search_pred.pl - Run a model (HMM or PWM) on whole genome DNA sequence using a pattern matching engine and 
                 report the predictions in GFF format. 
SYNOPSIS
Search_pred.pl --file seq.fna --engine HMMER/patser --model HMM/PWM
  
DESCRIPTION
This script scans a whole genome DNA sequence for a model. It uses HMMER for a hidden markov model and patser for a position weight matrix model. A hidden markov model can be created using the hmmbuild tool in the HMMER 2.3.2 suite. A position weight matrix can be built using a number of tools such as the make_matrix tool. The predictions are reported in GFF format. This GFF file is required by the Pred_cutoff.pl script to add in error rates for E-value estimation. 

NOTES
This script has been tested on Linux and requires patser (http://ural.wustl.edu/patser-v3e.1.tar.gz) and HMMER 2.3.2 (ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz) to be installed and accessible. We do not recommend using HMMER 3.0 since it is not optimized for DNA/DNA comparisons. You need to have the BioPerl library installed and accessible (http://www.bioperl.org/wiki/Installing_BioPerl).

COMMAND-LINE OPTIONS
Command-line options can be abbreviated to single-letter options. e.g. -f instead of --file. Some options are mandatory (see below).  
   --file         <file>      Fasta file whole genome DNA sequence (required)
   --engine       <engine>    Tool used to create the report (HMMER or patser) (required)
   --model        <HMM/PWM>   Model created from motifs. Hidden Markov Model (.hmm) for
                              HMMER or Position Weight Matrix(.pwm) for patser (required) 
   --domE         <x>         Domain E-value cutoff for hmmsearch on artificial genomes. Default 
                              is 0.3. All predictions with lower E-values will be recorded.
   --score        <x>         Score cutoff for patser on artificial genomes. Default is 1. All 
                              predictions with higher scores will be recorded.
   --verbose                  Print progress messages



CREATE_ART_GENOMES.PL

Create_art_genomes.pl - Create optimally ordered artificial genomes using the seqpp toolkit 

SYNOPSIS
Create_art_genomes.pl --file seq.fna --copies copies
  
DESCRIPTION
This script reads in a Fasta formatted DNA sequence representing the whole genome sequence of a organism. It requires the seqpp toolkit (http://stat.genopole.cnrs.fr/seqpp/) to be installed. It uses the estim_m tool to compute the optimal order for the markov chain to simulate the genome. The simul_m tool is then used to create artificial genomes that are placed inside the genomes directory. 

NOTES
Please make sure that the estim_m and simul_m tools are accessible from the command-line in the directory where you are running this script. The Fasta file should have a single DNA sequence. The artificial genomes created will be stored in genomes directory. This script has been tested on Linux.

COMMAND-LINE OPTIONS
Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options are mandatory (see below).  
   --file     <file>      Fasta file whole genome DNA sequence (mandatory)
   --copies   <copies>    Number of artificial genomes to create. Default is 600
   --verbose              Print progress messages



PRED_CUTOFF.PL

Pred_cutoff.pl - Add E value estimates to HMMER/patser report (GFF) and prints out 
                 a new GFF file. 

SYNOPSIS
  % Pred_cutoff.pl --file seq.fna --copies copies --engine HMMER/patser --model <HMM/PWM> --domE <x> --score <x> --report report.gff
  
DESCRIPTION
This script reads in a HMMER (for HMM search) or patser (for PWM search) report file in GFF format. The 
appropriate search engine (HMMER or patser) is then run against artificial genomes created using 
Create_art_genomes.pl script. E value estimates are computed for each prediction and the predictions are 
reported in GFF format.

NOTES
Please make sure this script is run in the directory containing the genomes directory created by
Create_art_genomes.pl. This script requires patser (http://ural.wustl.edu/patser-v3e.1.tar.gz) and
HMMER 2.3.2 (ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz) to be installed and 
in your path. We do not recommend using HMMER 3.0 since it is not optimized for DNA/DNA comparisons. 
You need to have the BioPerl library installed and accessible (http://www.bioperl.org/wiki/Installing_
BioPerl). This script has been tested on Linux.

COMMAND-LINE OPTIONS
Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --file         <file>      Fasta file whole genome DNA sequence (required)
   --copies       <copies>    Number of artificial genomes created using Create_art_genomes (required)
   --engine       <engine>    Tool used to create the report (HMMER or patser) (required)
   --model        <HMM/PWM>   Model created from motifs. Hidden Markov Model (.hmm) for
                              HMMER or Position Weight Matrix(.pwm) for patser  
   --report       <GFF>       File name of GFF file generated by Search_pred.pl on whole 
                              genome DNA sequence (required)
   --domE         <x>         Domain E-value cutoff for hmmsearch on artificial genomes. Default 
                              is 0.3. All predictions with lower E-values will be recorded.
   --score        <x>         Score cutoff for patser on artificial genomes. Default is 1. All 
                              predictions with higher scores will be recorded.
   --verbose                  Print progress messages
   

FILTER_PREDICTIONS.PL

Filter_predictions.pl - Filter out the predictions that are located within a user defined range 
                        of genomic features such as genes.

SYNOPSIS
Filter_predictions.pl --predictions pred.gff --features features.gff --range range
  
DESCRIPTION
This script scans the list of predictions and simply pulls out those within a user defined range. The input files need to be in GFF3 format such as those produced by Pred_cutoff.pl.

COMMAND-LINE OPTIONS
Command-line options can be abbreviated to single-letter options. e.g. -f instead of --file. All options are mandatory (see below).  
   --predictions  <GFF>       Predictions in GFF3 format (required)
   --features     <GFF>       Genomic features such as genes in GFF3 format (required)
   --range        <range>     Window for qualifying predictions (required)
   --verbose                  Print progress messages

