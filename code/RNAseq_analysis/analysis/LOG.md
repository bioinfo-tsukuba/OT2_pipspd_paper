# LOG

## 2024-03-25

### Activate conda environment
```bash
conda activate labopt
```

### Mounting the external ssd drive

restart the computer

### Create metadata.tsv
```bash
ls ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv
```

### Merge fastq files
- fastq files is separated by different lanes like L001, L002, L003, L004

```bash
bash ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/00_merge_fastq.sh
```

Log file is here
"K:\20240325_pipspd2nd\LAB_465_FASTQ\log_merge_fastq_20240325-144322.txt"


### Execute 01_fastp.sh
- Activate conda environment
```bash
conda activate labopt
```

- Install fastp
```bash
conda install -c bioconda -y fastp
```

```bash
bash ..OT2_yeast/code/RNAseq_analysis/analysis_2nd/01_fastp.sh
```

Log file is here
"K:\20240325_pipspd2nd\fastp\log_fastp_20240325-151137.log"

### Execute 02_STAR.sh
- Activate conda environment
```bash
conda activate labopt
```

- Install STAR
```bash
conda install -c bioconda -y star
```

- gunzip read files
```bash
gunzip *.fastq.gz
```

- create gtf file by gffread
```bash
gffread /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff \
-g /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa \
-T -o /mnt/k/20240325_pipspd2nd/reference_genome_W303/transcripts.gtf


gffread /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff \
-T -o /mnt/k/20240325_pipspd2nd/reference_genome_W303/transcripts.gtf

rsem-gff3-to-gtf --make-genes-as-transcripts /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff /mnt/k/20240325_pipspd2nd/reference_genome_W303/rsem-gff3-to-gtf.gtf
```


- create index
- for index
```bash
	STAR --runMode genomeGenerate --genomeDir /mnt/k/20240325_pipspd2nd/STAR/index \
    --genomeFastaFiles /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa \
    --sjdbGTFfile /mnt/k/20240325_pipspd2nd/reference_genome_W303/rsem-gff3-to-gtf.gtf \
    --runThreadN 16 --sjdbOverhang 100 --sjdbGTFfeatureExon exon --genomeSAindexNbases 10
```
- for index_2
```bash
	STAR --runMode genomeGenerate --genomeDir /mnt/k/20240325_pipspd2nd/STAR/index_2 \
    --genomeFastaFiles /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa \
    --sjdbGTFfile /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff \
    --runThreadN 16 --sjdbOverhang 100 --sjdbGTFfeatureExon CDS --genomeSAindexNbases 10
```




- Execute STAR
```bash
bash ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/02_STAR.sh
```

### Execute 03_RSEM.sh
install RSEM
```bash
git clone https://github.com/deweylab/RSEM.git
cd RSEM
make
make install
```

- prepare reference
```bash
    rsem-prepare-reference --gff3 \
    /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff \
    --num-threads 16 --gff3-genes-as-transcripts \
    /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa \
    /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000
```     


        /usr/lib/rna-star/bin/STAR-avx2 --runThreadN 16 --runMode genomeGenerate --genomeDir /mnt/k/20240325_pipspd2nd/RSEM --genomeFastaFiles /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa --sjdbGTFfile /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000.gtf --sjdbOverhang 100 --outFileNamePrefix /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000
        /usr/bin/STAR  --runThreadN 16  --runMode genomeGenerate  --genomeDir /mnt/k/20240325_pipspd2nd/RSEM  --genomeFastaFiles /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa  --sjdbGTFfile /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000.gtf  --sjdbOverhang 100  --outFileNamePrefix /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000


        rsem-gff3-to-gtf --make-genes-as-transcripts /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000.gtf --extract-sequences /mnt/k/20240325_pipspd2nd/RSEM/output.fa
### install gffread
conda install -c bioconda -y gffread

gffread /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff -g /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa -E -T -o /mnt/k/20240325_pipspd2nd/reference_genome_W303/transcripts.gtf

### Create reference genome
- Download the reference genome from https://www.yeastgenome.org/strain/w303
- Create a folder named `referenge_genome_W303` in `../OT2_yeast/code/RNAseq_analysis/analysis_2nd/reference_genome_W303/`



###　違う戦略をとりはじめる。hisat2+stringtieでやってみる

#### install some tools
```bash
conda install -c bioconda -y hisat2
conda install -c bioconda -y samtools
conda install -c bioconda -y stringtie
```

#### create index
```bash
mkdir /mnt/k/20240325_pipspd2nd/HISAT2/index/ -p
hisat2-build /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa /mnt/k/20240325_pipspd2nd/HISAT2/index/ -p 16
```
#### execute hisat2
```bash
sh ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/04_hisat2.sh
```
hisat2は1サンプル当たり30分かかるのでやめた。
それよりもstarでのindex fileで5139となるような条件を探そう
https://kimbio.info/rna-seq%E3%81%AE%E3%83%9E%E3%83%83%E3%83%94%E3%83%B3%E3%82%B0%E3%83%84%E3%83%BC%E3%83%AB%E3%81%AFstar%E3%81%8Bhisat2%E3%81%8B/
あとhisat2よりstarの方が良いという記事もある


#### 違う戦略を取り始める
- RSEMのreferenceを作るコード``でSTARでアライメント（マッピング）を行うと遺伝子数が5139となるので、その条件でRSEMで定量化を行う。
```bash
rsem-gff3-to-gtf --make-genes-as-transcripts /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD.gff /mnt/k/20240325_pipspd2nd/RSEM/reference_W303_SGD_2015_JRIU00000000.gtf
```
これから得たgtfを使って、STARでアライメントを行う

```bash
STAR --runMode genomeGenerate --genomeDir /mnt/k/20240325_pipspd2nd/STAR/index \
--genomeFastaFiles /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa \
--sjdbGTFfile /mnt/k/20240325_pipspd2nd/reference_genome_W303/rsem-gff3-to-gtf.gtf \
--runThreadN 16 --sjdbOverhang 100 --sjdbGTFfeatureExon exon --genomeSAindexNbases 10
```

- このindexを用いてSTARのアライメントを行う
```bash
bash ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/02_STAR.sh
```

- その後、RSEMで定量化を行う
```bash
bash ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/03_RSEM.sh
```

- カウントデータにまとめる
```bash
cd /mnt/k/20240325_pipspd2nd/RSEM
rsem-generate-data-matrix S{01..24}.genes.results > count_numread_RSEM.csv
```


#### Install IGV to WSL2
```bash
conda install -c bioconda igv -y
```


#### ついでにsalmonもやっておく
```bash
mamba install -c bioconda -y salmon
salmon index -p 4 -t /mnt/k/20240325_pipspd2nd/reference_genome_W303/W303_JRIU00000000_SGD_orf_genomic.fsa -i /mnt/k/20240325_pipspd2nd/salmon/index
sh ../OT2_yeast/code/RNAseq_analysis/analysis_2nd/90_salmon.sh
```

## 2024-03-28
### Modify GeneExpressionMatrix generated by salmon
The Matrix is here.
"K:\20240325_pipspd2nd\salmon\count_by_NumRead_salmon_after_fastp.csv"
