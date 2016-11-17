#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2013 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

#combn_sub is a function to take a combinatorial sampling of x choose m with a maximum nset of samples
#the code was found on the internet and was said to be a modification of the standard combn code - seems to work`
combn_sub <- function (x, m, nset = 1000, seed=123, simplify = TRUE, ...) {
#    print("Entering combn_sub")
#    print(x)
#    print(m)
#    print(nset)
#    print(seed)
#    print(simplify)
    stopifnot(length(m) == 1L)
    if (m < 0) 
        stop("m < 0", domain = NA)
    if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == 
        x) 
        x <- seq_len(x)
#    print(x)
    n <- length(x)
#    print(n)
    if (n < m) 
        stop("n < m", domain = NA)
    m <- as.integer(m)
    e <- 0
    h <- m
    a <- seq_len(m)
#    print(a)
    len.r <- length(r <-  x[a] )
#    print(len.r)
#    print(r)
    count <- as.integer(round(choose(n, m)))
    if (is.na(count)){
       count <- .Machine$integer.max
    }
#    print(count)
    if( count < nset ) nset <- count
    dim.use <- c(m, nset) 
#    print(dim.use)      

    ##-----MOD 1: Change the output matrix size--------------
    out <- matrix(r, nrow = len.r, ncol = nset, byrow = FALSE) 
#    print(out)

    if (m > 0) {
        if (count < (10 * nset)) {
#	    print("small")
            i <- 1L
            nmmp1 <- n - m + 1L

            ##----MOD 2: Select a subset of indices
            set.seed(seed)
            samp <- sort(sample( 1:count, nset ))  

            ##----MOD 3: Start a counter.
            counter <- 1L    

            while (a[1L] != nmmp1 ) {

                #-----MOD 4: Whenever the counter matches an index in samp, 
                #a combination of row indices is produced and stored in the matrix `out`
                if(samp[i] == counter){ 
                    out[, i] <- x[a]
                    if( i == nset ) break
                    i <- i + 1L
                }
		#iterate through all combinations
                if (e < n - h) {
                    h <- 1L
                    e <- a[m]
                    j <- 1L
                } else {
                    e <- a[m - h]
                    h <- h + 1L
                    j <- 1L:h
                }
                a[m - h + j] <- e + j
                #-----Increase the counter by 1 for each iteration of the while-loop
                counter <- counter + 1L
            }
        } else {
#	    print("large")
	    ##-----MOD 5: For large count just sample twice as many as needed and remove duplicates
	    i <- 0L
	    extra <- 2L * nset
	    outr <- matrix(r, nrow = extra, ncol = len.r, byrow = TRUE) 
#	    print(outr)
	    notDone = TRUE
	    while (notDone) {
	        while (i < extra) {
		    i <- i + 1L
            	    outr[i,] <- sort(sample(x, m))
		}
		outr <- unique(outr)
		i <- nrow(outr)
		if (i >= nset) {
		    notDone = FALSE
		    outr <- outr[1:nset,]
		}
	    }
	    out <- t(outr)
	}
    }
    array(out, dim.use)
}

library("compiler")
combn_sub <- cmpfun(combn_sub)

library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character",
	"perc_core", "p", 2, "character", #percenage of genomes to be considered core
	"perc_novel", "q", 2, "character", #percenage of genomes to be considered novel
	"num_samples", "s", 1, "character" #maximum number of combinatorial samples to run
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited 0/1 pan-genome cluster table with column and row headers>\n",
	"	-o <output file name for combinatorial counts>\n",
	"	-p <percentage of genomes to be considered core - default 100>\n",
	"	-q <percentage of genomes to be considered novel - default 0 but clearly at least one genome>\n",
	"	-s <maximum number of combinatorial samples to generate>\n",
	"\n",	
	"This script will generate a table of combinatorial counts for subsets of the pan-genome.\n",
	"\n",
	"\n");

if((!length(opt$input_file)) | (!length(opt$output_file))){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n\n");

OutputFileName=opt$output_file;

cat("\n")
cat("Output File Name: ", OutputFileName, "\n\n");

if (is.null(opt$perc_core)) {
    PercCore = 100;
} else {
    PercCore=as.numeric(opt$perc_core);
}

cat("\n")
cat("Percentage of genomes to be considered core: ", PercCore, "\n\n");

if (is.null(opt$perc_novel)) {
    PercNovel = 0;
} else {
    PercNovel=as.numeric(opt$perc_novel);
}

cat("\n")
cat("Percentage of genomes to be considered novel: ", PercNovel, "\n\n");

NumSamples=as.integer(opt$num_samples);

cat("\n")
cat("Maximum Number of Combinatorial Samples: ", NumSamples, "\n\n");

###############################################################################
###############################################################################

# Load data - assumes all values in table are integer to save memory on input colClasses="integer"
# header=FALSE assumes there is no header line with column names
#row.names=1 speicfies that the first column is a name for the cluster/row
z <- read.delim(InputFileName, sep="\t", header=FALSE, colClasses="integer", check.names=FALSE, comment.char="", quote="", row.names=1)
z[z>0] <- 1 #set all paralog values greater than 1 to 1 to indicate presence or absence making later calculations easier
#get the number of columns/genomes not including the cluster identifiers do not hard code this as done before
n <- ncol(z)
r <- nrow(z)
start <- 2
stop <- n
write(c("Genome", "Core", "Dispose", "Novel", "Unique", "Pan_genome"), file=OutputFileName, append = FALSE, sep = "\t", ncolumns=6)
for (i in start:stop) {
    num_core <- trunc((i * PercCore) / 100)
    if ((i < 10) && (num_core < i)) num_core <- num_core + 1 # correction doesn't really work for small numbers of genomes so adjust
    if (num_core < 1) num_core <- 1
    num_novel <- ceiling((i * PercNovel) / 100)
    if (num_novel <= 1) num_novel <- 1
    print(i) # print out the number of genomes in the permutation.
    print(num_core)
    print(num_novel)
    ic <- combn_sub(n, i, nset=NumSamples) # get all (maximum NumSamples) i combinations of n genomes
#    print(ic)
    num_combs <- ncol(ic)
    for (j in 1:num_combs) {
        if ((((i == 2) || (i == n)) && (num_combs < NumSamples)) || ((i == (n - 1)) && (((n - 1) * n) < NumSamples))) { # for small sample sizes we need to use all genomes as the last added not just a random one
	    # this will generate redundant values for core, dis, spec, and pan-genome but will not affect the later median calculation
	    start_genome <- 1
	    last_genome <- i
	} else {
  	    last_genome <- sample(i, 1) #designate an arbitrary one of the i genomes as the last one added
	    start_genome <- last_genome
	}
	for (k in start_genome:last_genome) {
    	    next_comb <- ic[,j]
	    sub_z <- z[,next_comb] #a subset of columns/genomes from z
	    if (i == 1) { # for only one genome rowSums does not work - silly
	        cluster_sums <- sub_z #number of genomes that have a gene in this cluster
	        last_z <- sub_z
	    } else {
	        cluster_sums <- rowSums(sub_z) #number of genomes that have a gene in this cluster
	        last_z <- sub_z[,k] == 1
	    }
	    core <- length(cluster_sums[cluster_sums >= num_core])
	    spec <- length(cluster_sums[(cluster_sums <= num_novel) & (cluster_sums > 0)])
	    spec_v <- cluster_sums == 1
	    not <- length(cluster_sums[cluster_sums == 0])
	    dis <- length(cluster_sums[(cluster_sums < num_core) & (cluster_sums > num_novel)])
	    last_v <- last_z & spec_v
	    last <- length(last_v[last_v == TRUE])
	    if (i == 1) { # for only one genome pan_genome just equals core
	        pan_genome <- core
	    } else {
	        pan_genome <- core + dis + spec
	    }
	    if ((pan_genome + not) != r) {
	        print("Categories of genes do not sum up to total number of genes")
	   	print(next_comb)
	   	print(core)
	   	print(dis)
	   	print(spec)
	   	print(not)
	   	print(r)
	   	print(last)
	   	print(k)
	    }
	    write(c(i, core, dis, last, spec, pan_genome), file=OutputFileName, append = TRUE, sep = "\t", ncolumns=6)
	}
    }
}
