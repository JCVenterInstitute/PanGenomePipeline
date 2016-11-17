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
	"zero_flag", "z", 0, "logical", # if specified eliminate rows that are all 0/blank/null or sum to zero
	"col_only_flag", "C", 0, "logical", # if specified only remove columns not rows
	"input_file", "i", 1, "character", # file containing matrix/data frame/table with row and column headers
	"output_file", "o", 1, "character", # output file name to write subset of input matrix
	"column_labels", "c", 1, "character", # column names to keep
	"row_labels", "r", 2, "character" # row names to keep (optional if not specified then same as column names
)

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE)

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2]

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited matrix with row and column headers>\n",
	"	-o <output file name for subset of matrix>\n",
	"	-c <labels of columns to keep>\n",
	"	-r <labels of rows to keep (optional) use column labels if not specified>\n",
	"	-C <flag to indicate that rows are not being deleted by label>\n",
	"	-z <flag to indicate that rows/cols whose rows/cols are all 0/blank/null should be deleted>\n",
	"\n",	
	"This script will generate the specified subset of the inputted matrix.\n",
	"\n",
	"\n")

if((!length(opt$input_file)) | (!length(opt$output_file))){
	cat(usage)
	q(status=-1)
}

###############################################################################

OutputFileName <- opt$output_file
#OutputFileNameTXT <- paste(OutputFileName, ".txt", sep="")

options(width=160)
#sink(OutputFileNameTXT)

InputFileName <- opt$input_file

cat("\nInput File Name: ", InputFileName, "\n\n")
cat("\nOutput File Name: ", OutputFileName, "\n\n")

ColsToKeep <- scan(file=opt$column_labels, what="character")

cat("\nColumns to keep: ", "\n")
print(ColsToKeep)

if (is.null(opt$row_labels)) {
    RowsToKeep <- ColsToKeep
} else {
    RowsToKeep <- scan(file=opt$row_labels, what="character")
}
if (is.null(opt$cols_only_flag)) {
    cat("\nRows to keep: ", "\n")
    print(RowsToKeep)
}

###############################################################################
###############################################################################

# Load data - assumes all values in table are numeric to save memory on input colClasses="numeric"
# header=TRUE assumes there is a header line with column genome names
# row.names=1 specifies that the first column is a genome name for the row
input_dataframe <- read.delim(InputFileName, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1, strip.white=TRUE)
if (is.null(opt$cols_only_flag)) {
    sub_dataframe <- input_dataframe[RowsToKeep,ColsToKeep]
} else {
    sub_dataframe <- input_dataframe[ ,ColsToKeep]
}
write.table(sub_dataframe, file=OutputFileName, sep="\t", quote=FALSE)
