#!/usr/bin/env Rscript

library(ape)

dist.dna.alt <- function(fa, as.matrix=F)
{
	fa.mat <- matrix(toupper(unlist(as.character(fa))), nrow = length(fa), byrow = T)
	out.mat <- matrix(NA, nrow(fa.mat), nrow(fa.mat))
	rownames(out.mat) = names(fa);
	colnames(out.mat) = names(fa);	
	for (i in 1:(nrow(fa.mat)-1))
	{
		for (j in (i+1):(nrow(fa.mat)))
		{
			if (ncol(fa.mat) > 1)
			{
				fa.trim = fa.mat[c(i,j), colSums(fa.mat[c(i,j),] == "-" | fa.mat[c(i,j),] == "N") < 2];
			}
			else
			{
				fa.trim = fa.mat[c(i,j),sum(fa.mat[c(i,j),] == "-" | fa.mat[c(i,j),] == "N") < 2];
			}
			if (length(fa.trim) > 2)
			{
				out.mat[i,j] = mean(fa.trim[1,] != fa.trim[2,]);
			}
			else
			{
				out.mat[i,j] = sum(fa.trim[1] != fa.trim[2]);
			}
			out.mat[j,i] = out.mat[i,j];
		}
	}
	if (as.matrix == F)
	{	return(as.dist(out.mat)); }
	return(out.mat);
}

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  cat("\nAt least one argument must be supplied.\nThe first argument (required) is the input multi-fasta file.\nThe second field is optional, and is the name of a genome you wish to get statistics for individually.\n
       e.g. summarize_alignment.R /path/to/multifasta target_genome \n\n")
  quit()()
}


input_file = args[1]
target_genome = args[2]



in_aligned_fasta <- read.FASTA(input_file)
distance_matrix <- dist.dna.alt(in_aligned_fasta, as.matrix = T)
#Adding in NAs to remove the diagonal
distance_matrix_nodiag = distance_matrix;
for (i in 1:nrow(distance_matrix_nodiag)) { distance_matrix_nodiag[i,i] = NA; }




find_unique_allels <- function(fa, max.gap = 0.5)
{
	fa.mat <- matrix(toupper(unlist(as.character(fa))), nrow = length(fa), byrow = T)
	fa.mat.trim = fa.mat[,colSums(fa.mat == "-")/nrow(fa.mat) <= max.gap];
	fa.unique.row.cnt = rep(0, nrow(fa.mat));
	for (bse in c("A", "C", "T", "G"))
	{
		#cat(paste(bse, "\n"));
		if (sum(colSums(fa.mat == bse)==1)> 1)
		{
			fa.unique.row.cnt <- fa.unique.row.cnt + (rowSums(fa.mat[,(1:ncol(fa.mat))[colSums(fa.mat == bse)==1]]==bse));
		}
		if (sum(colSums(fa.mat == bse)==1) == 1)
		{
			fa.unique.row.cnt <- fa.unique.row.cnt + (0+(fa.mat[,(1:ncol(fa.mat))[colSums(fa.mat == bse)==1]]==bse))
		}
	}
	names(fa.unique.row.cnt) = names(fa);
	return(fa.unique.row.cnt);
}

get_stats <- function(data){
  header <- c("Mean", "Median", "Min", "Max", "SD", "0%", "5%", "10%", "15%", "20%", "25%",
              "30%", "35%", "40%", "45%", "50%", "55%", "60%", "65%", "70%", "75%",
              "80%", "85%", "90%", "95%", "100%")
  
  data_list <- c()
  #print(data)
  data_list <- append(data_list, c(mean(data, na.rm=T), median(data, na.rm=T), min(data, na.rm=T), max(data, na.rm=T), sd(data, na.rm=T)))
  quant_data <- unname(quantile(data, c(0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1), na.rm=T))
  data_list <- append(data_list, quant_data)
  out.df <- t(data.frame(data_list))
  colnames(out.df) <- header
  return(out.df)
}

check_target_genome <- function(df, target_genome){
  if (!(target_genome %in% row.names(df))){
    cat("Not a valid target genome. Check how the first five genomes are formatted and then adjust your target genome to reflect this format:\n")
    cat(head(row.names(df)),"\n")
    quit()
  }
}

if(is.na(target_genome)){
  out_df <- get_stats(distance_matrix_nodiag)
  row.names(out_df) <- "All"
  out_df <- round(out_df, 5)
  #Making the matrix of values:
  add.df <- get_stats(find_unique_allels(in_aligned_fasta));
  row.names(add.df) <- "UniqueAlleleCount";
  out_df <- rbind(out_df, round(add.df, 5))
  write.table(out_df, file = "", quote = F, col.names = NA, sep = "\t")
} else{
  check_target_genome(distance_matrix, target_genome)
  out_df <- get_stats(distance_matrix_nodiag[target_genome,])
  col_names = colnames(out_df)
  row.names(out_df) <- target_genome
  out_df <- cbind(round(out_df, 5), find_unique_allels(in_aligned_fasta)[target_genome])
  colnames(out_df) = c(col_names, "UniqueAlleleCount")
  write.table(out_df, file = "", quote = F, col.names = NA, sep = "\t")
}

