FROM python:3.9-alpine

RUN apk add build-base
RUN pip install biopython

# clustalw
RUN apk add g++ wget make
RUN wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz
RUN tar xvzf clustalw-2.1.tar.gz
RUN cd clustalw-2.1; ./configure; make; make install

# muscle
RUN wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz; \
    tar xvzf muscle3.8.31_i86linux64.tar.gz; \
    rm muscle3.8.31_i86linux64.tar.gz; mv muscle3.8.31_i86linux64 /bin/muscle

RUN mkdir src
WORKDIR /src