#!/usr/bin/env bash

module load gcc
module load R/3.5.1

git clone --recursive https://github.com/rwdavies/STITCH.git
cd STITCH
./scripts/install-dependencies.sh
cd releases
wget https://github.com/rwdavies/stitch/releases/download/1.6.2/STITCH_1.6.2.tar.gz
R CMD INSTALL STITCH_1.6.2.tar.gz
