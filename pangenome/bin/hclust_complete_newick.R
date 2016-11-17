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

# This snippet of R code should take a complete linkage hclust object and generate a Newick formatted string (label:branch_length,label:branch_length)
# where label is either a leaf node label or recursively replaced by a subtree to generate a binary tree
# hclust is generating distances at each merge step
# assumptions are that this distance should be equally divided between the two clusters (subtrees) being merged,
# each subtree already has a partial distance equal to its height and leaf nodes have zero height
# an initial join of two leaf nodes results in a subtree with two equal length branches of 1/2 the merge distance
# new branch lengths equal 1/2 the merge distance - the subtree height
clus <- hclust(dist_mat)
print(clus$merge)
print(clus$height)
print(clus$labels)
num_labels <- nrow(clus$labels)
num_merges <- nrow(clus$merge)
subtree_height <- numeric(num_merges)
merge_labels <- character(num_merges)
print(subtree_height)
print(merge_labels)
for (i in 1:num_merges) {
    #print(i)
    #print(clus$merge[i,1])
    height <- clus$height[i] / 2.0
    subtree_height[i] <- height
    if (clus$merge[i,1] < 0) {
        label1 <- clus$labels[-1 * clus$merge[i,1]]
	height1 <- height
    } else {
	#print(subtree_height[clus$merge[i,1]])
	#print(merge_labels[clus$merge[i,1]])
	label1 <- merge_labels[clus$merge[i,1]]
	height1 <- height - subtree_height[clus$merge[i,1]]
    }
    #print(clus$merge[i,2])
    if (clus$merge[i,2] < 0) {
        label2 <- clus$labels[-1 * clus$merge[i,2]]
	height2 <- height
    } else {
	#print(subtree_height[clus$merge[i,2]])
	#print(merge_labels[clus$merge[i,2]])
	label2 <- merge_labels[clus$merge[i,2]]
	height2 <- height - subtree_height[clus$merge[i,2]]
    }
    merge_labels[i] <- paste("(", label1, ":", height1, ",", label2, ":", height2, ")", sep="")
}
print(merge_labels[num_merges])
