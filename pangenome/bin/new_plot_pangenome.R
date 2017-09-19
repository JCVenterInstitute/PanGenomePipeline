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
	"verb_on", "v", 0, "binary"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited combinatorial pan-genome counts>\n",
	"	-o <output file name for pdf plots and txt files>\n",
	"	-v <turns on verbose setting, i.e. printing summaries to the standard output>\n", 
	"\n",	
	"This script will generate PDF and txt files for exponential and power law models of Core, Novel, and Pan-genome gene sizes.\n",
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

OutputFileNamePDF <- paste(OutputFileName, ".pdf", sep="")
OutputFileNameTXT <- paste(OutputFileName, ".txt", sep="")

###############################################################################
###############################################################################

pdf(OutputFileNamePDF)
sink(OutputFileNameTXT)
# attempt to get more tick marks
par (lab = c(10, 10, 7))
###Novel:Power:mean###
# only plot medians - comment out means
z <- read.delim(InputFileName, sep="\t", header=TRUE, colClasses="integer", check.names=FALSE, comment.char="", quote="")
znox<-z[,"Genome"]
min_genomes <- min(znox)
max_genomes <- max(znox)
znoy<-z[,"Novel"]
min_novel <- min(znoy)
max_novel <- max(znoy)
#y <- vector(mode="numeric", length=0)
#x <- vector(mode="integer", length=0)
ymed <- vector(mode="numeric", length=0)
xmed <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
#    y <- c(y, mean(z[,"Novel"][z[,"Genome"]==i]))
#    x <- c(x, i)
    ymed <- c(ymed, median(z[,"Novel"][z[,"Genome"]==i]))
    xmed <- c(xmed, i)
}
#print("Novel Genes : Power Model: Means")
#nov_pwm_nlsfit<-nls(y~k*x^(-a), start=list(k=min_genomes*max(y),a=1),trace=TRUE)
#summary(nov_pwm_nlsfit)
#rsq.nls<-1-(sum((residuals(nov_pwm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Novel Genes : Power Model: Medians")
nov_pwmed_nlsfit<-nls(ymed~k*xmed^(-a), start=list(k=min_genomes*max(ymed),a=1),trace=TRUE,control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(nov_pwmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(nov_pwmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
maxXlim <- 4 * max_genomes
minYlim <- min_novel
if (minYlim < 1) { #for log-log plots cannot have this less than 1
    minYlim <- 1
}
maxYlim <- max_novel
plot (znoy ~ znox, xlim=c(2,maxXlim), ylim=c(minYlim,maxYlim), log="xy", col="gray", lend="square", main="Novel Genes Medians: Power Red; Exponential Blue", xlab="# of genomes", ylab="# of new genes", pch=1, cex=0.25)
points(xmed, ymed, col="purple", pch=1, cex=0.5)
#points(x, y, col="purple", pch=5, cex=0.5)
xf<-seq(from=1, to=maxXlim, by=1)
xf
#dxf <- data.frame(x=xf)
dxfmed <- data.frame(xmed=xf)
#lines(xf, predict(nov_pwm_nlsfit, dxf), col="red", lwd=1, lty="dashed")
lines(xf, predict(nov_pwmed_nlsfit, dxfmed), col="red", lwd=1)

###Novel:Exponential:mean###
#print("Novel Genes: Exponential Model: Means")
#nov_expm_nlsfit<-nls(y~a*exp(-x/b)+c, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
#summary(nov_expm_nlsfit)
#rsq.nls<-1-(sum((residuals(nov_expm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Novel Genes: Exponential Model: Medians")
nov_expmed_nlsfit<-nls(ymed~a*exp(-xmed/b)+c, start=list(a=exp(1)*(max(ymed)-min(ymed)),b=min_genomes,c=min(ymed)),trace=TRUE, )
if (!is.null(opt$verb_on)){
	summary(nov_expmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(nov_expmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
#lines(xf, predict(nov_expm_nlsfit, dxf), col="blue", lwd=1, lty="dashed")
lines(xf, predict(nov_expmed_nlsfit, dxfmed), col="blue", lwd=1)

###Core:Power:mean###
znoy<-z[,"Core"]
min_core <- min(znoy)
max_core <- max(znoy)
#y <- vector(mode="numeric", length=0)
#x <- vector(mode="integer", length=0)
ymed <- vector(mode="numeric", length=0)
xmed <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
#    y <- c(y, mean(z[,"Core"][z[,"Genome"]==i]))
#    x <- c(x, i)
    ymed <- c(ymed, median(z[,"Core"][z[,"Genome"]==i]))
    xmed <- c(xmed, i)
}
#print("Core Genes: Power Model: Means")
#cor_pwm_nlsfit<-nls(y~k*x^(-a), start=list(k=max(y),a=2/(min_genomes+max_genomes)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
#summary(cor_pwm_nlsfit)
#rsq.nls<-1-(sum((residuals(cor_pwm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Core Genes: Power Model: Medians")
cor_pwmed_nlsfit<-nls(ymed~k*xmed^(-a), start=list(k=max(ymed),a=2/(min_genomes+max_genomes)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(cor_pwmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(cor_pwmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
maxYlim <- max_core
minYlim <- min_core
if (minYlim < 1) { #for log-log plots cannot have this less than 1
    minYlim <- 1
}
plot (znoy ~ znox, xlim=c(2,maxXlim), ylim=c(minYlim, maxYlim), log="xy", col="gray", lend="square", main="Core Genes Medians: Power Red; Exponential Blue", xlab="# of genomes", ylab="# of core genes", pch=1, cex=0.25)
points(xmed, ymed, col="purple", pch=1, cex=0.5)
#points(x, y, col="purple", pch=5, cex=0.5)
#lines(xf,predict(cor_pwm_nlsfit, dxf),col="red", lwd=1, lty="dashed")
lines(xf,predict(cor_pwmed_nlsfit, dxfmed),col="red", lwd=1)

###Core:Exponential:mean###
#print("Core Genes: Exponential Model: Means")
#cor_expm_nlsfit<-nls(y~a*exp(-x/b)+c, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
#summary(cor_expm_nlsfit)
#rsq.nls<-1-(sum((residuals(cor_expm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Core Genes: Exponential Model: Medians")
cor_expmed_nlsfit<-nls(ymed~a*exp(-xmed/b)+c, start=list(a=exp(1)*(max(ymed)-min(ymed)),b=min_genomes,c=min(ymed)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(cor_expmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(cor_expmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
#lines(xf,predict(cor_expm_nlsfit, dxf), col="blue", lwd=1, lty="dashed")
lines(xf,predict(cor_expmed_nlsfit, dxfmed), col="blue", lwd=1)

###Pan_genome:Power:mean###
znoy<-z[,"Pan_genome"]
min_pan <- min(znoy)
max_pan <- max(znoy)
#y <- vector(mode="numeric", length=0)
#x <- vector(mode="integer", length=0)
ymed <- vector(mode="numeric", length=0)
xmed <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
#    y <- c(y, mean(z[,"Pan_genome"][z[,"Genome"]==i]))
#    x <- c(x, i)
    ymed <- c(ymed, median(z[,"Pan_genome"][z[,"Genome"]==i]))
    xmed <- c(xmed, i)
}
#print("Pan-Genome: Power Model: Means")
#pan_pwm_nlsfit<-nls(y~k*x^(a), start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
#summary(pan_pwm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Pan-Genome: Power Model: Medians")
pan_pwmed_nlsfit<-nls(ymed~k*xmed^(a), start=list(k=max(ymed)/min_genomes,a=1),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{	
	summary(pan_pwmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_pwmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
maxYlim <- 1.5 * max_pan
minYlim <- min_pan
if (minYlim < 1) { #for log-log plots cannot have this less than 1
    minYlim <- 1
}
plot (znoy ~ znox, xlim=c(2,maxXlim), ylim=c(minYlim, maxYlim), log="xy", col="gray", lend="square", main="Pan-Genome Medians: Power Red; Exponential Blue", xlab="# of genomes", ylab="Pan_genome size", pch=1, cex=0.25)
points(xmed, ymed, col="purple", pch=1, cex=0.5)
#points(x, y, col="purple", pch=5, cex=0.5)
#lines(xf,predict(pan_pwm_nlsfit, dxf), col="red", lwd=1, lty="dashed")
lines(xf,predict(pan_pwmed_nlsfit, dxfmed), col="red", lwd=1)

###Pan_genome:Exponential:mean###
#print("Pan-Genome: Exponential Model: Means")
#pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
#summary(pan_expm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Pan-Genome: Exponential Model: Medians")
pan_expmed_nlsfit<-nls(ymed~-a*exp(-xmed/b)+c, start=list(a=max(ymed),b=(min_genomes+max_genomes)/2,c=max(ymed)+min(ymed)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(pan_expmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_expmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
#lines(xf,predict(pan_expm_nlsfit, dxf),col="blue", lwd=1, lty="dashed")
lines(xf,predict(pan_expmed_nlsfit, dxfmed),col="blue", lwd=1)

###Unique:Power:mean###
znoy<-z[,"Unique"]
min_uniq <- min(znoy)
max_uniq <- max(znoy)
#y <- vector(mode="numeric", length=0)
#x <- vector(mode="integer", length=0)
ymed <- vector(mode="numeric", length=0)
xmed <- vector(mode="integer", length=0)
if (min_genomes < 2) {
    min_genomes <- 2 # unique genes don't make much sense for only 1 genome so messes up model fitting
}
for (i in min_genomes:max_genomes) {
#    y <- c(y, mean(z[,"Unique"][z[,"Genome"]==i]))
#    x <- c(x, i)
    ymed <- c(ymed, median(z[,"Unique"][z[,"Genome"]==i]))
    xmed <- c(xmed, i)
}
#print("Unique: Power Model: Means")
#pan_pwm_nlsfit<-nls(y~k*x^(a), start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
#summary(pan_pwm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Unique: Power Model: Medians")
pan_pwmed_nlsfit<-nls(ymed~k*xmed^(a), start=list(k=max(ymed)/min_genomes,a=1),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(pan_pwmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_pwmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
maxYlim <- 1.5 * max_uniq
minYlim <- min_uniq
if (minYlim < 1) { #for log-log plots cannot have this less than 1
    minYlim <- 1
}
plot (znoy ~ znox, xlim=c(2,maxXlim), ylim=c(minYlim, maxYlim), log="xy", col="gray", lend="square", main="Unique Medians: Power Red; Exponential Blue", xlab="# of genomes", ylab="Unique size", pch=1, cex=0.25)
points(xmed, ymed, col="purple", pch=1, cex=0.5)
#points(x, y, col="purple", pch=5, cex=0.5)
xf<-seq(from=2, to=maxXlim, by=1)
xf
#dxf <- data.frame(x=xf)
dxfmed <- data.frame(xmed=xf)
#lines(xf,predict(pan_pwm_nlsfit, dxf),col="red", lwd=1, lty="dashed")
lines(xf,predict(pan_pwmed_nlsfit, dxfmed),col="red", lwd=1)

###Unique:Exponential:mean###
#print("Unique: Exponential Model: Means")
#pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
#summary(pan_expm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Unique: Exponential Model: Medians")
pan_expmed_nlsfit<-nls(ymed~-a*exp(-xmed/b)+c, start=list(a=max(ymed),b=(min_genomes+max_genomes)/2,c=max(ymed)+min(ymed)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(pan_expmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_expmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
#lines(xf,predict(pan_expm_nlsfit, dxf),col="blue", lwd=1, lty="dashed")
lines(xf,predict(pan_expmed_nlsfit, dxfmed),col="blue", lwd=1)

###Dispose:Power:mean###
znoy<-z[,"Dispose"]
min_disp <- min(znoy)
max_disp <- max(znoy)
#y <- vector(mode="numeric", length=0)
#x <- vector(mode="integer", length=0)
ymed <- vector(mode="numeric", length=0)
xmed <- vector(mode="integer", length=0)
if (min_genomes < 3) {
    min_genomes <- 3 # dispensible genes are always equal to 0 bydefinition for 1 or 2 genomes so messes up model fitting
}
for (i in min_genomes:max_genomes) {
#    y <- c(y, mean(z[,"Dispose"][z[,"Genome"]==i]))
#    x <- c(x, i)
    ymed <- c(ymed, median(z[,"Dispose"][z[,"Genome"]==i]))
    xmed <- c(xmed, i)
}
#print("Dispensable: Power Model: Means")
#pan_pwm_nlsfit<-nls(y~k*x^(a), start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
#summary(pan_pwm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Dispensable: Power Model: Medians")
pan_pwmed_nlsfit<-nls(ymed~k*xmed^(a), start=list(k=max(ymed)/min_genomes,a=1),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(pan_pwmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_pwmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
maxYlim <- 1.5 * max_disp
minYlim <- min_disp
if (minYlim < 1) { #for log-log plots cannot have this less than 1
    minYlim <- 1
}
plot (znoy ~ znox, xlim=c(3,maxXlim), ylim=c(minYlim, maxYlim), log="xy", col="gray", lend="square", main="Dispensable Medians: Power Red; Exponential Blue", xlab="# of genomes", ylab="Dispensable size", pch=1, cex=0.25)
points(xmed, ymed, col="purple", pch=1, cex=0.5)
#points(x, y, col="purple", pch=5, cex=0.5)
xf<-seq(from=3, to=maxXlim, by=1)
xf
#dxf <- data.frame(x=xf)
dxfmed <- data.frame(xmed=xf)
#lines(xf,predict(pan_pwm_nlsfit, dxf),col="red", lwd=1, lty="dashed")
lines(xf,predict(pan_pwmed_nlsfit, dxfmed),col="red", lwd=1)

###Dispose:Exponential:mean###
#print("Dispensable: Exponential Model: Means")
#pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
#summary(pan_expm_nlsfit)
#rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
#rsq.nls
print("Dispensable: Exponential Model: Medians")
pan_expmed_nlsfit<-nls(ymed~-a*exp(-xmed/b)+c, start=list(a=max(ymed),b=(min_genomes+max_genomes)/2,c=max(ymed)+min(ymed)),trace=TRUE, control = nls.control(warnOnly=T, minFactor =0))
if (!is.null(opt$verb_on))
{
	summary(pan_expmed_nlsfit)
}
rsq.nls<-1-(sum((residuals(pan_expmed_nlsfit))^2)/sum((ymed-mean(ymed))^2))
rsq.nls
#lines(xf,predict(pan_expm_nlsfit, dxf),col="blue", lwd=1, lty="dashed")
lines(xf,predict(pan_expmed_nlsfit, dxfmed),col="blue", lwd=1)
