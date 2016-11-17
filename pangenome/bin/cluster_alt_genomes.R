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


library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character",
	"max_size", "m", 1, "character", #maximum size of a cluster
	"num_clusters", "n", 1, "character" #optimal number of clusters to generate
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited vectors with row and column headers>\n",
	"	-o <output file name for pdf plots, cluster and txt files>\n",
	"	-m <maximum size of a cluster>\n",
	"	-n <optimal number of clusters to generate>\n",
	"\n",	
	"This script will generate PDF and txt files for clusters based on hclust given the input vectors.\n",
	"\n",
	"\n");

if((!length(opt$input_file)) | (!length(opt$output_file))){
	cat(usage);
	q(status=-1);
}

###############################################################################

OutputFileName=opt$output_file;

OutputFileNamePDF <- paste(OutputFileName, ".pdf", sep="")
OutputFileNameTXT <- paste(OutputFileName, ".txt", sep="")
OutputFileNameCLU <- paste(OutputFileName, ".cluster", sep="")

options(width=160)
sink(OutputFileNameTXT)

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n\n");

cat("\n")
cat("Output File Name: ", OutputFileName, "\n\n");

NumClusters <- as.integer(opt$num_clusters);

cat("\n")
cat("Optimal number of clusters to generate: ", NumClusters, "\n\n");

MaxSize <- as.integer(opt$max_size)

cat("\n")
cat("Maximum size of clusters to generate: ", MaxSize, "\n\n");

MaxSize <- MaxSize / 2

###############################################################################
###############################################################################

# Load data - assumes all values in table are integer to save memory on input colClasses="integer"
# header=TRUE assumes there is a header line with column names
#row.names=1 specifies that the first column is a name for the vector/row
vecs <- read.delim(InputFileName, sep="\t", header=TRUE, colClasses=c("character", "integer"), check.names=FALSE, comment.char="", quote="", row.names=1)
#print(vecs)
num_genomes <- nrow(vecs)
width_pdf <- ceiling(num_genomes / 3)
if (width_pdf < 7) {
    width_pdf <- 7
}
pdf(OutputFileNamePDF, width=width_pdf)
dist_mat <- dist(vecs)
print(dist_mat)
clus <- hclust(dist_mat)
plot(clus)
plot(as.dendrogram(clus))
print(clus)
print(clus$merge)
print(clus$height)
print(clus$labels)
num_merges <- nrow(clus$merge)
cluster_size <- numeric(num_merges)
merge_labels <- character(num_merges)
level_batch <- rep(1, num_merges)
max_level <- numeric(num_merges)
cat(file=OutputFileNameCLU, append=FALSE, "")
for (i in 1:num_merges) {
    if (clus$merge[i,1] < 0) {
        size1 <- 1
        label1 <- clus$labels[-1 * clus$merge[i,1]]
	level1 <- 1
    } else {
        size1 <- cluster_size[clus$merge[i,1]]
	label1 <- merge_labels[clus$merge[i,1]]
	level1 <- max_level[clus$merge[i,1]]
    }
    if (clus$merge[i,2] < 0) {
        size2 <- 1
        label2 <- clus$labels[-1 * clus$merge[i,2]]
	level2 <- 1
    } else {
        size2 <- cluster_size[clus$merge[i,2]]
	label2 <- merge_labels[clus$merge[i,2]]
	level2 <- max_level[clus$merge[i,2]]
    }
    if (label1 == "") {
        next_label <- label2
    } else {
        if (label2 == "") {
	    next_label <- label1
	} else {
	    next_label <- paste(label1, label2, sep=",")
	}
    }
    cluster_size[i] <- size1 + size2
    merge_labels[i] <- next_label
    max_level[i] <- max(level1, level2)
    if ((cluster_size[i] > MaxSize) || (i == num_merges)) {
	cat(file=OutputFileNameCLU, append=TRUE, "L", max_level[i], "B", level_batch[max_level[i]], "(", merge_labels[i], ")\n", sep="")
        merge_labels[i] <- paste("L", max_level[i], "B", level_batch[max_level[i]], sep="")
	level_batch[max_level[i]] <- level_batch[max_level[i]] + 1
        max_level[i] <- max_level[i] + 1
	cluster_size[i] <- 1
    }
}
#print (level_batch)
#print (max_level)
#print (merge_labels)
#print (cluster_size)
