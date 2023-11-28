#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }
which tar >/dev/null || { echo "tar: command not found"; exit 1; }

mkdir -p data
cd data || exit

echo "Downloading refdata-gex-GRCh38-2020-A.tar.gz"
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz  \
  && md5sum -c - <<<"dfd654de39bff23917471e7fcc7a00cd  refdata-gex-GRCh38-2020-A.tar.gz" \
  && tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
rm refdata-gex-GRCh38-2020-A.tar.gz

