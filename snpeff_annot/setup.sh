#! /usr/bin/env bash

# ------------------- Activate the miniconda hook -----------------------------------

# NOTE: SnpEff requires a newer version of java which is newer than what we have on the server.
# Thus, a conda environment is used instead of the packages in /usr/local/packages folder.
# To allow access from multiple users, a miniconda instance is installed here: 
#   /local/data/Malaria/Projects/Takala-Harrison/AFRIMS/miniconda3
#
# To activate conda's base environment in your current shell session:
eval "$(/local/data/Malaria/Projects/Takala-Harrison/AFRIMS/miniconda3/bin/conda shell.bash hook)" 

# Install mamba for better resolve dependencies if not installed
# conda install mamba -n base -c conda-forge

# ------------------ Create the environment ---------------------
mamba create -n snpeff -c bioconda snpeff
conda activate snpeff

# ------------------ Add Reference Genome to SnpEff ---------------------
# NOTE: you can the version, strain or species as needed

## curdir
CURDIR=`pwd`

## Add a genome to SnpEff config file
PREFIX=$(dirname $(realpath $(which snpEff)))
SPECIES=Pfalciparum
STRAIN=3D7
VERSION=44
echo "# Pfalciparum3D7, version PlasmoDB-${VERSION}_${SPECIES}${STRAIN}
PlasmoDB-${VERSION}_${SPECIES}${STRAIN}.genome : ${SPECIES}${STRAIN}" >> $PREFIX/snpEff.config

## Make a data folder within SnpEff folder
mkdir -p $PREFIX/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}

## Download and rename files
cd $PREFIX/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}
wget https://plasmodb.org/common/downloads/release-${VERSION}/${SPECIES}${STRAIN}/fasta/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}_Genome.fasta -O sequences.fa
wget https://plasmodb.org/common/downloads/release-${VERSION}/${SPECIES}${STRAIN}/gff/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}.gff -O genes.gff
wget https://plasmodb.org/common/downloads/release-${VERSION}/${SPECIES}${STRAIN}/fasta/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}_AnnotatedCDSs.fasta -O cds.fa
wget https://plasmodb.org/common/downloads/release-${VERSION}/${SPECIES}${STRAIN}/fasta/data/PlasmoDB-${VERSION}_${SPECIES}${STRAIN}_AnnotatedProteins.fasta -O protein.fa

## fix header in protein.fa 
#   1) change separator from | to space 
#   2) remove the '-p1' suffix in protein id and make it consistent with 
#   cds transcript ids
 sed -e 's/ | / /g;s/-p[0-9]* / /' -i protein.fa
 cd ../../

## Build the genome
snpEff build -gff3 -v PlasmoDB-${VERSION}_${SPECIES}${STRAIN} 
cd ${CURDIR}
