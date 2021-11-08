Installation
===================

Please note that users need to configure the software paths based on their own installed paths in `config.json`. All processed should be run under the root path of `LncRNA_pipeline`.

# Install or Load conda
```
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-py37_4.10.3-Linux-x86_64.sh
```
OR
```
module load miniconda/2018-02-09
```

# Create conda environment
conda create -n LncRNApipeline python=2.7.16
source activate LncRNApipeline

# CPAT installation
conda install -c bioconda cpat==1.2.4

# LGC installation
git clone https://github.com/GuangyWang/LGC

# PLEK installation
conda install -c bioconda plek=1.2

# hisat2 installation
conda install -c bioconda hisat2==2.0.5

# sambamba installation
conda install -c bioconda sambamba==0.6.6

# stringtie installation
conda install -c bioconda stringtie==1.3.4

# gffread installation
conda install -c bioconda gffread==0.9.8

# kallisto installation
conda install -c bioconda kallisto==0.44.0

# bedtools installation
conda install -c bioconda bedtools==2.28.0

# gffcompare installation
conda install -c bioconda gffcompare==0.9.8

# pysam installation
conda install -c bioconda pysam==0.15.2
conda install -c anaconda cython

# htseq python2.7 support
git clone https://github.com/simon-anders/htseq
cd htseq
python setup.py install

# Strawberry installation
wget -c https://github.com/ruolin/strawberry/releases/download/v1.1.2/strawberry
chmod 755 ./strawberry

wget -c https://ftp.gnu.org/gnu/libc/glibc-2.14.1.tar.gz

# R packages installation 
conda install -c bioconda bioconductor-genomicfeatures
conda install -c bioconda r-stringr
conda install -c r r-ggplot2


# Install / load singularity (version 3.7.1)
module load singularity/3.7.1

# Use singularity
singularity pull --arch amd64 library://zixiuli910/default/centos:0.1 (work in August/2021, don't work in Nov/2021)

for example, 
        ./centos_0.1.sif cpat.py       # 1.2.4
        ./centos_0.1.sif PLEK.py       # 1.2
        ./centos_0.1.sif hisat2        # 2.0.5
        ./centos_0.1.sif sambamba      # 0.6.6
        ./centos_0.1.sif stringtie     # 1.3.4d
        ./centos_0.1.sif gffread       # 0.9.8c
        ./centos_0.1.sif kallisto      # 0.44.0
        ./centos_0.1.sif bedtools      # v2.28.0
        ./centos_0.1.sif gffcompare    # v0.9.8
        ./centos_0.1.sif strawberry    # v1.1.2
        ./centos_0.1.sif htseq-count   # 0.11.5
```
singularity pull --arch amd64 library://zixiuli910/default/centos:0.2 (work in August/2021, don't work in Nov/2021)

for example,
        ./centos_0.2.sif R             # 3.6.0
        ./centos_0.2.sif Rscript
```

singularity pull --arch amd64 library://zixiuli910/default/centos:0.3 (work in August/2021, don't work in Nov/2021)

for example,
        ./centos_0.3.sif lgc-1.0.0.py  # 1.0
        cd DATA_LIB/CPPred_20180516/bin && /home/zl10w/Liz/LncRNA_pipeline_Python/LncRNA_pipeline/centos_0.3.sif python CPPred.py

# Usage
python src/cluster_job_submission.py --input_file INPUT_FILE --output_dir OUTPUT_DIR

for example,
python src/cluster_job_submission.py --input_file /nl/umw_chan_zhou/Liz/LncRNA_pipeline_Python/input/HSC_GSE68108_WT.txt --output_dir /home/zl10w/Liz/TEST
```        
 
Five columns are required in the INPUT_FILE
col1: 
col2:
col3:
col4:
col5:


local path/absoluate path for input_file=
absoluate path for OUTPUT_DIR
