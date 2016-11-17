#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2015 J. Craig Venter Institute (JCVI).                  #
#       All rights reserved.                                                  #
#       Written by Derrick E. Fouts, Ph.D.                                    #
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
	"input_file", "i", 1, "character"
	);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited core_clusters_histogram.txt generated from overview_stats.txt>\n",
        "\n",	
	"This script will generate a core_cluster_histogram.pdf file, showing the PanOCT cluster size distribution from singletons (min) to core genes (max).\n",
	"\n",
	"\n");

if((!length(opt$input_file))){
	cat(usage);
	q(status=-1);
}

###############################################################################

InputFileName=opt$input_file;

cat("\n")
cat("Input File Name: ", InputFileName, "\n\n");

OutputFileName="core_cluster_histogram";

OutputFileNamePDF <- paste(OutputFileName, ".pdf", sep="")

cat("\n")
cat("Output File Name: ", OutputFileNamePDF, "\n\n");



###############################################################################
###############################################################################

pdf(OutputFileNamePDF)

# read in datafile
x <- as.matrix(read.delim(InputFileName, sep="\t", header=FALSE, colClasses=c("numeric","numeric"), check.names=FALSE))

# generate histogram

#plot(x[,1], x[,2], type = "h", ylab="# of clusters", xlab="Cluster Size",main = "PanOCT Cluster Size Distribution")
plot(x[,1], x[,2], type = "h", ylab="# of clusters", xlab="Cluster Size",main = "PanOCT Cluster Size Distribution", bty='L', xlim = range(pretty(c(x[1,1], x[,1]))), ylim = range(pretty(c(0, x[,2]))), lwd = 1, lend = "square")