# Call:  R --slave --no-save --no-restore --no-environ --silent --args 100 < gpd.R 
# count of art genomes is 100


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

rm(list=ls())

#data
counts = read.table("hitcounts.txt")
x = counts[,1]
y = counts[,2]
x2 = x^2
count = as.numeric(commandArgs(TRUE));

#lm() fit
lmfit_x2=lm(y~x+x2)
coefs_x2=coefficients (lmfit_x2)
score_cutoff_x2 = (-sqrt(-4*coefs_x2[1]*coefs_x2[3] + (coefs_x2[2]^2) + 4*coefs_x2[3]*count) - coefs_x2[2]) / (2*coefs_x2[3])

hits=scan("allhits.txt", quiet = TRUE)
#GPD analysis using evir
#install evir if required
if(!require("evir", quietly=TRUE)) { 
 install.packages("evir", repos = "http://cran.r-project.org")
}

library(evir)
gpd=gpd(hits,threshold=score_cutoff_x2)

#print Evalues to file
cutoff=score_cutoff_x2
#cat(file="new.Evals.out","Score   Eval\n")
for (score in seq(min(x),max(x),0.1)){
  Eval=((1+((gpd$par.ests[1]*(score-cutoff))/gpd$par.ests[2]))^-(1/gpd$par.ests[1]))
  cat(file="Evals.out",append=TRUE, score," ",Eval,"\n")
}

