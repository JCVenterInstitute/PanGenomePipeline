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


library('getopt')

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 1, "character",
	"hclust_method", "m", 1, "character",
	"dist_convert", "d", 2, "character", #optional conversion of similarity to distance matrix using offset provided
	"threshold_clusters", "h", 2, "character", #optional threshold of clusters to generate
	"num_clusters", "n", 2, "character" #optional number of clusters to generate
)

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE)

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2]

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited distance matrix with row and column headers>\n",
	"	-o <output file name for pdf plots, cluster and txt files>\n",
	"	-n <optional number of medoids to generate either n or h must be specified but not both>\n",
	"	-h <optional threshold for medoids to generate either n or h must be specified but not both>\n",
	"	-d <optional offset to convert a similarity matrix to a distance matrix>\n",
	"	-m <method to use for hclust such as complete, single or average>\n",
	"\n",	
	"This script will generate genome medoids for the given number or threshold specified.\n",
	"\n",
	"\n")

if((!length(opt$input_file)) | (!length(opt$output_file))){
	cat(usage)
	q(status=-1)
}

###############################################################################

OutputFileName <- opt$output_file

#OutputFileNamePDF <- paste(OutputFileName, ".pdf", sep="")
#OutputFileNameTXT <- paste(OutputFileName, ".txt", sep="")
OutputFileNameMED <- paste(OutputFileName, ".medoids", sep="")
OutputFileNameNWK <- paste(OutputFileName, ".newick", sep="")

options(width=160)
#sink(OutputFileNameTXT)

InputFileName <- opt$input_file

cat("\nInput File Name: ", InputFileName, "\n\n")
cat("\nOutput File Name: ", OutputFileName, "\n\n")

if (is.null(opt$num_clusters)) {
    if (is.null(opt$threshold_clusters)) {
        cat("\nMust specify either -n or -h\n")
	cat(usage)
	q(status=-1)
    } else {
        ThreshMedoids <- as.numeric(opt$threshold_clusters)
    	cat("\nThreshold of medoids to generate: ", ThreshMedoids, "\n\n")
    }
} else {
    if (is.null(opt$threshold_clusters)) {
        NumMedoids <- as.integer(opt$num_clusters)
    	cat("\nNumber of medoids to generate: ", NumMedoids, "\n\n")
    } else {
        cat("\nCannot specify both -n and -h\n")
	cat(usage)
	q(status=-1)
    }
}

HclustMethod <- opt$hclust_method

cat("\nHclust Method: ", HclustMethod, "\n\n")

###############################################################################
###############################################################################

# Load data - assumes all values in table are numeric to save memory on input colClasses="numeric"
# header=TRUE assumes there is a header line with column genome names
#row.names=1 specifies that the first column is a genome name for the row
distance_dataframe <- read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1, strip.white=TRUE, fill=TRUE)
#print(distance_dataframe)
num_genomes <- nrow(distance_dataframe)
width_pdf <- ceiling(num_genomes / 3)
if (width_pdf < 7) {
    width_pdf <- 7
}
#pdf(OutputFileNamePDF, width=width_pdf)
distance_matrix <- as.matrix(distance_dataframe)
#print(distance_matrix)
if (!is.null(opt$dist_convert)) {
    cat("\nConverting similarity matrix to distance matrix using threshold: ", opt$dist_convert, "\n\n")
    distance_matrix <- as.numeric(opt$dist_convert) - distance_matrix
}
if (!is.na(distance_matrix[nrow(distance_matrix),1]))
{
	if (is.na(distance_matrix[nrow(distance_matrix), nrow(distance_matrix)]))
	{
		distance_matrix= rbind(rep(NA, nrow(distance_matrix)+1), cbind( distance_matrix, rep(NA, nrow(distance_matrix)),))
		colnames(distance_matrix)[ncol(distance_matrix)] = rownames(distance_matrix)[ncol(distance_matrix)]
		rownames(distance_matrix)[1] = colnames(distance_matrix)[1]
		num_genomes = num_genomes + 1; # Adding one to the genome count
	
	}
	dist_mat <- as.dist(distance_matrix)
	distance_matrix = as.matrix(dist_mat)

}else
{
	#Adding blank rows to the boundaries to correct for missing identity
	distance_matrix= rbind(cbind(rep(NA, nrow(distance_matrix)), distance_matrix), rep(NA, nrow(distance_matrix)+1))
	
	colnames(distance_matrix)[1] = rownames(distance_matrix)[1]
	rownames(distance_matrix)[nrow(distance_matrix)] = colnames(distance_matrix)[nrow(distance_matrix)]
	num_genomes = num_genomes + 1; # Adding one to the genome count
	dist_mat <- as.dist(t(distance_matrix))
	distance_matrix = as.matrix(dist_mat)
}
if (HclustMethod == "nj") {
    library('ape')
    nj_tree <- nj(dist_mat)
    write.tree(nj_tree,file=OutputFileNameNWK)
    q()
} else {
    clus <- hclust(dist_mat, method = HclustMethod)
}
# This snippet of R code should take a complete linkage hclust object and generate a Newick formatted string (label:branch_length,label:branch_length)
# where label is either a leaf node label or recursively replaced by a subtree to generate a binary tree
# hclust is generating distances at each merge step
# assumptions are that this distance should be equally divided between the two clusters (subtrees) being merged,
# each subtree already has a partial distance equal to its height and leaf nodes have zero height
# an initial join of two leaf nodes results in a subtree with two equal length branches of 1/2 the merge distance
# new branch lengths equal 1/2 the merge distance - the subtree height
num_labels <- nrow(clus$labels)
num_merges <- nrow(clus$merge)
subtree_height <- numeric(num_merges)
merge_labels <- character(num_merges)
for (i in 1:num_merges) {
    height <- clus$height[i] / 2.0
    subtree_height[i] <- height
    if (clus$merge[i,1] < 0) {
        label1 <- clus$labels[-1 * clus$merge[i,1]]
	height1 <- height
    } else {
	label1 <- merge_labels[clus$merge[i,1]]
	height1 <- height - subtree_height[clus$merge[i,1]]
    }
    if (clus$merge[i,2] < 0) {
        label2 <- clus$labels[-1 * clus$merge[i,2]]
	height2 <- height
    } else {
	label2 <- merge_labels[clus$merge[i,2]]
	height2 <- height - subtree_height[clus$merge[i,2]]
    }
    merge_labels[i] <- paste("(", label1, ":", height1, ",", label2, ":", height2, ")", sep="")
}
write(merge_labels[num_merges],file=OutputFileNameNWK)
#plot(clus)
#plot(as.dendrogram(clus))
#print(clus)
#print(clus$merge)
#print(clus$height)
#print(clus$labels)
if (is.null(opt$num_clusters)) {
    clusters <- cutree(clus, h=ThreshMedoids)
} else {
    clusters <- cutree(clus, k=NumMedoids)
}
#print(clusters)
sum_check <- 0
NumMedoids <- max(clusters)
size_cluster <- numeric(NumMedoids)
for (i in 1:NumMedoids) {
    size_cluster[i] <- sum(clusters == i)
    sum_check <- sum_check + size_cluster[i]
}
#print(size_cluster)
if (sum_check != num_genomes) {
    print("Sum of cluster sizes not equal to number of genomes")
    print(sum_check)
    print(num_genomes)
}
cat(file=OutputFileNameMED, append=FALSE, "")

clust.medoid = function(i, pdistmat, pclusters) {
    ind = (pclusters == i)
    if (sum(ind) <= 1) {
        return (rownames(distance_matrix)[ind]) ## medoid of a single object is the object
    } else {
        return(names(which.min(rowSums( pdistmat[ind, ind] ))))
    }
}

medoids <- sapply(unique(clusters), clust.medoid, as.matrix(dist_mat), clusters)
#print(medoids)

for (i in 1:NumMedoids) {
    cat(medoids[i], size_cluster[i], file=OutputFileNameMED, append=TRUE, sep="\t")
    cat("\t", file=OutputFileNameMED, append=TRUE)
    cat(clus$labels[clusters == i], file=OutputFileNameMED, append=TRUE, sep=",")
    cat("\n", file=OutputFileNameMED, append=TRUE)
}

