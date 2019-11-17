FROM ubuntu:16.04

RUN apt-get update \
  && apt-get install -y python-pip python-dev \
  && apt-get install -y hmmer \
  && apt-get install -y mafft \
  && apt-get autoremove \
  && apt-get clean


ADD scripts /usr/local/bin
RUN chmod +x /usr/local/bin/*
