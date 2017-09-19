#!/usr/bin/env Rscript

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
 