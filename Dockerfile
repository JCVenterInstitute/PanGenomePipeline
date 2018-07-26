# Start from ubuntu
FROM ubuntu

# Set up apt-get and install a few things we need to install other things
ENV DEBIAN_FRONTEND=noninteractive LANG=en_US.UTF-8 LC_ALL=C.UTF-8 LANGUAGE=en_US.UTF-8
RUN [ "apt-get", "-q", "update" ]
RUN [ "apt-get", "-qy", "--force-yes", "upgrade" ]
RUN [ "apt-get", "-qy", "--force-yes", "dist-upgrade" ]
RUN [ "apt-get", "install", "-qy", "--force-yes", \
      "perl", \
      "build-essential", \
      "cpanminus", \
      "wget", \
      "r-base", \
      "libxml2-dev", \
      "hmmer", \
      "python2.7", \
      "python-pip" ]
RUN [ "apt-get", "clean" ]
RUN [ "rm", "-rf", "/var/lib/apt/lists/*", "/tmp/*", "/var/tmp/*" ]

# Of course, python isn't that simple to install :/
RUN echo 'alias python="python2.7"' >> ~/.bashrc
RUN pip install biopython
#RUN pip install --upgrade pip
#RUN pip install Bio

# Install R packages needed
RUN R -e "install.packages('getopt', repos = 'http://cran.us.r-project.org')"

# Install all the perl modules needed.
RUN ["cpanm", "Capture::Tiny", "Term::ReadKey", "MLDBM", "Devel::InnerPackage", "Class::Load", "String::RewritePrefix", "XML::LibXML", "HTML::TableExtract" ]

# Install blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.2.31+-x64-linux.tar.gz
RUN rm /ncbi-blast-2.2.31+-x64-linux.tar.gz
ENV PATH=".:/ncbi-blast-2.2.31+/bin:${PATH}"

# Install the pipeline.
COPY ./pangenome /pangenome
