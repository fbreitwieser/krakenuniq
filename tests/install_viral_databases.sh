#!/bin/bash

set -xeu

mkdir -p dbs/refseq-viral/library
mkdir -p dbs/refseq-viral-plus/library

[[ -L dbs/refseq-viral/taxonomy ]] || ln -s data/taxonomy dbs/refseq-viral
[[ -L dbs/refseq-viral/library/viral ]] || ln -s data/library/viral/ dbs/refseq-viral/library
[[ -L dbs/refseq-viral-plus/library/viral ]] || ln -s data/library/viral/ dbs/refseq-viral-plus/library
[[ -L dbs/refseq-viral-plus/library/viral-neighbors ]] || ln -s data/library/viral-neighbors/ dbs/refseq-viral-plus/library

export PATH="install:$PATH"
krakenu-build --db refseq-viral --build

