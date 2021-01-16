FROM python:3.8-slim-buster AS base

RUN apt-get update && apt-get install -y \
  git \
  unzip \
  && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/MineralsCloud/qha.git \
  && cd qha \
  && git checkout 479bce9a16ec6cc75874b734c2f5a038acbb05d2 \
  && pip3 install --no-cache-dir . \
  && cd .. \
  && rm -r qha

RUN git clone https://github.com/MineralsCloud/cij.git \
  && cd cij \
  && git checkout a878a654158db0de5942698bd4df301b7764d812 \
  && pip3 install --no-cache-dir . \
  && cd .. \
  && rm -r cij
