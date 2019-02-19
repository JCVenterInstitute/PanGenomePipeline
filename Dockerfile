# Start from centos
FROM centos

# Set up yum and install a few things we need to install other things
RUN yum update -y; 
RUN yum groupinstall -y 'Development Tools'
RUN yum install -y perl cpanminus wget libcurl-devel openssl-devel libxml2-devel epel-release centos-release-scl ruby gd gd-devel libpng12 transfig
RUN yum install -y R python27 python-pip xfig
RUN yum clean all

# Of course, python isn't that simple to install :/
RUN echo 'alias python="python2.7"' >> ~/.bashrc
RUN pip install biopython

# Install blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.2.31+-x64-linux.tar.gz; rm /ncbi-blast-2.2.31+-x64-linux.tar.gz
ENV PATH=".:/ncbi-blast-2.2.31+/bin:${PATH}"

# Install hmmer.  Wish yum install -y hmmer would work for centos, not just fedora :/
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz; tar xzf hmmer.tar.gz; cd hmmer-3.2.1; ./configure; make; make install; rm /hmmer.tar.gz

# Install R packages needed
RUN R -e "install.packages('getopt', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ape', repos = 'http://cran.us.r-project.org')"

# Install all the perl modules needed.
RUN ["cpanm", "Capture::Tiny", "Term::ReadKey", "DB_File", "MLDBM", "Devel::InnerPackage", "Class::Load", "String::RewritePrefix", "Fatal", "XML::LibXML", "HTML::TableExtract", "LWP::UserAgent", "File::Fetch", "File::Slurp@9999.19"", "Bio::SeqIO", "DBI", "GD", "Getopt::Long::Descriptive" ]

# Install the pipeline.
RUN wget https://github.com/JCVenterInstitute/PanGenomePipeline/archive/master.zip && unzip /master.zip; ln -s /PanGenomePipeline-master/pangenome /pangenome; rm /master.zip

# Retrieve the data directory
RUN wget https://sandbox.zenodo.org/record/237583/files/HMMER2GO_data.tgz?download=1 -O /HMMER2GO_data.tgz && tar -zxf /HMMER2GO_data.tgz -C /pangenome/bin/HMMER2GO/ ; rm /HMMER2GO_data.tgz

# Installing EMBOSS for a few useful binaries
RUN wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz && tar xf EMBOSS-6.6.0.tar.gz; cd /EMBOSS-6.6.0; ./configure; make; make install; rm /EMBOSS-6.6.0.tar.gz
