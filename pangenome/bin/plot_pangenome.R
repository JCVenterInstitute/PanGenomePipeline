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
	"output_file", "o", 1, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];

usage = paste (
	"\nUsage:\n\n", script_name,
	"\n",
	"	-i <input file name for tab delimited combinatorial pan-genome counts>\n",
	"	-o <output file name for pdf plots and txt files>\n",
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
###Novel:Power:mean###
print("Novel Genes: Power Model: Means")
z <- read.delim(InputFileName, sep="\t", header=TRUE, colClasses="integer", check.names=FALSE, comment.char="", quote="")
znox<-z[,"Genome"]
min_genomes <- min(znox)
max_genomes <- max(znox)
znoy<-z[,"Novel"]
zno<-data.frame(Genome=znox,Novel=znoy)
y <- vector(mode="numeric", length=0)
x <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
    y <- c(y, mean(z[,"Novel"][z[,"Genome"]==i]))
    x <- c(x, i)
}
nov_pwm_nlsfit<-nls(y~k*x^(-a), data=z, start=list(k=min_genomes*max(y),a=1),trace=TRUE)
summary(nov_pwm_nlsfit)
rsq.nls<-1-(sum((residuals(nov_pwm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(1,maxYlim), col="blue", lend="square", main="Novel Genes: Power Model: Means", xlab="# of genomes", ylab="# of new genes", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(nov_pwm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Novel:Exponential:mean###
print("Novel Genes: Exponential Model: Means")
nov_expm_nlsfit<-nls(y~a*exp(-x/b)+c, data=z, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
summary(nov_expm_nlsfit)
rsq.nls<-1-(sum((residuals(nov_expm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(1,maxYlim), col="blue", lend="square", main="Novel Genes: Exponential Model: Means", xlab="# of genomes", ylab="# of new genes", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(nov_expm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Novel:Power###
print("Novel Genes: Power Model: All")
no_nlsfit<-nls(znoy~k*znox^(-a), data=z, start=list(k=min_genomes*max(y),a=1),trace=TRUE)
summary(no_nlsfit)
rsq.nls<-sum((predict(no_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(1,maxYlim), col="blue", lend="square", main="Novel Genes: Power Model: All", xlab="# of genomes", ylab="# of new genes", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(no_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Novel:Exponential###
print("Novel Genes: Exponential Model: All")
no_nlsfit<-nls(znoy~a*exp(-znox/b)+c, data=z, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
summary(no_nlsfit)
rsq.nls<-sum((predict(no_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(1,maxYlim), col="blue", lend="square", main="Novel Genes: Exponential Model: All", xlab="# of genomes", ylab="# of new genes", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(no_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Core:Power:mean###
print("Core Genes: Power Model: Means")
znoy<-z[,"Core"]
zno<-data.frame(Genome=znox,Core=znoy)
y <- vector(mode="numeric", length=0)
x <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
    y <- c(y, mean(z[,"Core"][z[,"Genome"]==i]))
    x <- c(x, i)
}
cor_pwm_nlsfit<-nls(y~k*x^(-a), data=z, start=list(k=max(y),a=2/(min_genomes+max_genomes)),trace=TRUE)
summary(cor_pwm_nlsfit)
rsq.nls<-1-(sum((residuals(cor_pwm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Core Genes: Power Model: Means", xlab="# of genomes", ylab="# of core genes", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(cor_pwm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Core:Exponential:mean###
print("Core Genes: Exponential Model: Means")
cor_expm_nlsfit<-nls(y~a*exp(-x/b)+c, data=z, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
summary(cor_expm_nlsfit)
rsq.nls<-1-(sum((residuals(cor_expm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Core Genes: Exponential Model: Means", xlab="# of genomes", ylab="# of core genes", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(cor_expm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Core:Power###
print("Core Genes: Power Model: All")
cor_pw_nlsfit<-nls(znoy~k*znox^(-a), data=z, start=list(k=max(y),a=2/(min_genomes+max_genomes)),trace=TRUE)
summary(cor_pw_nlsfit)
rsq.nls<-sum((predict(cor_pw_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Core Genes: Power Model: All", xlab="# of genomes", ylab="# of core genes", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(cor_pw_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Core:Exponential###
print("Core Genes: Exponential Model: All")
cor_exp_nlsfit<-nls(znoy~a*exp(-znox/b)+c, data=z, start=list(a=exp(1)*(max(y)-min(y)),b=min_genomes,c=min(y)),trace=TRUE)
summary(cor_exp_nlsfit)
rsq.nls<-sum((predict(cor_exp_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Core Genes: Exponential Model: All", xlab="# of genomes", ylab="# of core genes", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(cor_exp_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Pan_genome:Power:mean###
print("Pan-Genome: Power Model: Means")
znoy<-z[,"Pan_genome"]
zno<-data.frame(Genome=znox,Pan_genome=znoy)
y <- vector(mode="numeric", length=0)
x <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
    y <- c(y, mean(z[,"Pan_genome"][z[,"Genome"]==i]))
    x <- c(x, i)
}
pan_pwm_nlsfit<-nls(y~k*x^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pwm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Pan-Genome: Power Model: Means", xlab="# of genomes", ylab="Pan_genome size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_pwm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Pan_genome:Exponential:mean###
print("Pan-Genome: Exponential Model: Means")
pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_expm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Pan-Genome: Exponential Model: Means", xlab="# of genomes", ylab="Pan_genome size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_expm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Pan_genome:Power###
print("Pan-Genome: Power Model: All")
pan_pw_nlsfit<-nls(znoy~k*znox^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pw_nlsfit)
rsq.nls<-sum((predict(pan_pw_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Pan-Genome: Power Model: All", xlab="# of genomes", ylab="Pan_genome size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_pw_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Pan_genome:Exponential###
print("Pan-Genome: Exponential Model: All")
pan_exp_nlsfit<-nls(znoy~-a*exp(-znox/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_exp_nlsfit)
rsq.nls<-sum((predict(pan_exp_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Pan-Genome: Exponential Model: All", xlab="# of genomes", ylab="Pan_genome size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_exp_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Dispose:Power:mean###
print("Dispensable: Power Model: Means")
znoy<-z[,"Dispose"]
zno<-data.frame(Genome=znox,Dispose=znoy)
y <- vector(mode="numeric", length=0)
x <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
    y <- c(y, mean(z[,"Dispose"][z[,"Genome"]==i]))
    x <- c(x, i)
}
pan_pwm_nlsfit<-nls(y~k*x^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pwm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Dispensable: Power Model: Means", xlab="# of genomes", ylab="Dispensable size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_pwm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Dispose:Exponential:mean###
print("Dispensable: Exponential Model: Means")
pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_expm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Dispensable: Exponential Model: Means", xlab="# of genomes", ylab="Dispensable size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_expm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Dispose:Power###
print("Dispensable: Power Model: All")
pan_pw_nlsfit<-nls(znoy~k*znox^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pw_nlsfit)
rsq.nls<-sum((predict(pan_pw_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Dispensable: Power Model: All", xlab="# of genomes", ylab="Dispensable size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_pw_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Dispose:Exponential###
print("Dispensable: Exponential Model: All")
pan_exp_nlsfit<-nls(znoy~-a*exp(-znox/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_exp_nlsfit)
rsq.nls<-sum((predict(pan_exp_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Dispensable: Exponential Model: All", xlab="# of genomes", ylab="Dispensable size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_exp_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Unique:Power:mean###
print("Unique: Power Model: Means")
znoy<-z[,"Unique"]
zno<-data.frame(Genome=znox,Unique=znoy)
y <- vector(mode="numeric", length=0)
x <- vector(mode="integer", length=0)
for (i in min_genomes:max_genomes) {
    y <- c(y, mean(z[,"Unique"][z[,"Genome"]==i]))
    x <- c(x, i)
}
pan_pwm_nlsfit<-nls(y~k*x^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pwm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_pwm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Unique: Power Model: Means", xlab="# of genomes", ylab="Unique size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_pwm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

###Unique:Exponential:mean###
print("Unique: Exponential Model: Means")
pan_expm_nlsfit<-nls(y~-a*exp(-x/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_expm_nlsfit)
rsq.nls<-1-(sum((residuals(pan_expm_nlsfit))^2)/sum((y-mean(y))^2))
rsq.nls
maxXlim <- 1.5 * max_genomes
maxYlim <- 2 * max(y)
minYlim <- 0.5 * min(y)
plot (y ~ x, xlim=c(1,maxXlim), ylim=c(minYlim, maxYlim), col="blue", lend="square", main="Unique: Exponential Model: Means", xlab="# of genomes", ylab="Unique size", pch=3, cex=0.5)
xf<-seq(0,maxXlim,length=maxXlim)
lines(xf,predict(pan_expm_nlsfit,list(x=xf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Unique:Power###
print("Unique: Power Model: All")
pan_pw_nlsfit<-nls(znoy~k*znox^(a), data=z, start=list(k=max(y)/min_genomes,a=1),trace=TRUE)
summary(pan_pw_nlsfit)
rsq.nls<-sum((predict(pan_pw_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Unique: Power Model: All", xlab="# of genomes", ylab="Unique size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_pw_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)

##Unique:Exponential###
print("Unique: Exponential Model: All")
pan_exp_nlsfit<-nls(znoy~-a*exp(-znox/b)+c, data=z, start=list(a=max(y),b=(min_genomes+max_genomes)/2,c=max(y)+min(y)),trace=TRUE)
summary(pan_exp_nlsfit)
rsq.nls<-sum((predict(pan_exp_nlsfit)-mean(znoy))^2)/sum((znoy-mean(znoy))^2) 
rsq.nls
plot (znoy ~ znox, xlim=c(1,maxXlim), ylim=c(minYlim,maxYlim), col="blue", lend="square", main="Unique: Exponential Model: All", xlab="# of genomes", ylab="Unique size", pch=3, cex=0.5)
znoxf<-seq(0,maxXlim,length=maxXlim)
lines(znoxf,predict(pan_exp_nlsfit,list(znox=znoxf)),col="red", lwd=1)
#points(zno, cex=0.3)
