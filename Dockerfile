FROM ubuntu:14.04

MAINTAINER Alban Gaignard <alban.gaignard@univ-nantes.fr>

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-get update

RUN apt-get install -y git curl wget bzip2 libcurl4-openssl-dev libxml2-dev

RUN apt-get install -y --force-yes r-base-dev

## install R packages
RUN Rscript -e 'install.packages("devtools",dependencies=TRUE, repos="https://cloud.r-project.org/")' \
	&& Rscript -e 'install.packages("optparse", dependencies=TRUE, repos="https://cloud.r-project.org/")' \
	&& Rscript -e 'install.packages("xml", dependencies=TRUE, repos="https://cloud.r-project.org/")'

## install additional R packages using R
RUN > rscript.R \
	&& echo 'source("https://bioconductor.org/biocLite.R")' >> rscript.R \
	&& echo 'biocLite(ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("edgeR",ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("MASS",ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("DESeq2",ask=FALSE)' >> rscript.R \
	&& echo 'biocLite("preprocessCore",ask=FALSE)' >> rscript.R \
	&& Rscript rscript.R

# Cleanup
RUN rm rscript.R

WORKDIR /home

SHELL ["/bin/bash", "-c"]

COPY Readme.md .
COPY src .

#RUN /home/testMicroArraySimul.sh
#RUN /home/testRnaSeqSimul.sh
#RUN /home/testNormRnaseq.sh
#RUN /home/testFusionSimul.sh

ENTRYPOINT [ "bash" ]
