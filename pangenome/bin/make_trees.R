#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (C) 2016-2017 J. Craig Venter Institute (JCVI).             #
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

argList = commandArgs(trailingOnly=TRUE);

if(length(argList)!= 2)
{
	cat("Need filename & type")
	stop(0);
}
file.name <- argList[1];
tree.type <- tolower(argList[2]);

if (!file.exists(file.name))
{
    cat(paste("\nFile does not ", file.name," exist...\n\nQuiting..\n"));
    q()
}
 
distance_dataframe <- read.delim(file.name, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1, strip.white=TRUE)
distance_matrix <- as.matrix(distance_dataframe)

dist_mat <- as.dist(distance_matrix);

accept_hclust <- c("complete", "average", "single", "nj", "upgma", "neighbor")
require(ape, quietly=TRUE);

if (tree.type %in% accept_hclust)
{
	if (tree.type == "nj" | tree.type == "neighbor")
      {
		
			
			phy <- nj(dist_mat)
		
      }
	  else
      {
		if (tree.type == "upgma")
		{
			tree.type = "average"
		}
        hc_phy <- hclust(dist_mat, method=tree.type)
        phy <- as.phylo(hc_phy)
      }
      
    }else
    {
      cat("\nSuggested unacceptible phylogenetic method...\n\nQuiting..\n");
      q()
    }
  
 write.tree(phy, file="outtree")
 
