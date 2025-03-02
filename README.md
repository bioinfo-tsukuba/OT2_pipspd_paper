# Investigating the Effects of Robotic OT-2 Parameter Variations on Yeast Growth and Gene Expression: A Comprehensive Study Using RNA-seq and Spot Assays 

This repository provides code and data to reproduce the experiments described in the paper "Investigating the Effects of Robotic OT-2 Parameter Variations on Yeast Growth and Gene Expression: A Comprehensive Study Using RNA-seq and Spot Assays" by Taguchi et al. [url]
DOI:

### Overview

This study comprehensively investigated the effects of pipetting speed on yeast growth and gene expression using the robotic OT-2 platform, RNA-seq, and spot assays.

### Features

* Reproduction of experimental protocols using OT-2
* Preprocessing and analysis of RNA-seq data
* Quantitative growth assay of yeast colonies using BaQFA
* Comparison of yeast growth and gene expression under different pipetting speeds

### Requirements

* Python: 3.10
* Libraries listed in `requirements.txt`:
    * pandas
    * numpy
    * STAR
    * RSEM
    * fastp
    * baQFA

### Usage for OT-2 protocol

##### Requirements
Prepare conda/mamba environment by running the following command:

```bash
mamba create --name labopt --file ./requiments.yaml
conda activate labopt
```

##### Prepare `.env` file

1. Make a file named `.env` in `experiment` directory
2. Write the following items in it
3. Make sure that .env file is placed in the same directory as `run_experiment.ipynb`.


```txt
# Path to `labopt` directory (where `chef` and `saucier` directories are located)
PATH_LABOPT=...

# Device path of Transporter GEN2 (can be checked by Arduino IDE)
GEN2_DEVICEPATH=...

# Path to SSH secret key for accessing to OT-2
OT2_SSH_KEY=...

# IP address of OT-2
OT2_IP_ADDRESS=...
```


#### Usage for RNA-seq analysis

1. Run the bash scripts step by step in `code/RNAseq_analysis/analysis/` to analyze the RNA-seq data.

2. RNA-seq analysis results are output to the `results/240328_RNAseq2nd` directory.

#### Usage for BaQFA analysis
1. Prepare the time-series image of yeast colonies obtained by the flatbed scanner.

2. Run the Jupyter Notebooks `code/BaQFA_protocols/001_rotate-crop.ipynb` and `code/BaQFA_protocols/002_run-baqfa.ipynb` to perform the quantitative growth assay of yeast colonies using BaQFA.

3. BaQFA analysis results are output to the `results/240415_baQFA` directory.


### Notes

* This repository provides code and data to reproduce the experiments in the paper and does not guarantee the results of the paper.
* Refer to the official Opentrons documentation for OT-2 operation.

### License

This repository is released under the MIT License.

### Acknowledgements

This implementation is heavily influenced by the paper "Investigating the Effects of Robotic OT-2 Parameter Variations on Yeast Growth and Gene Expression: A Comprehensive Study Using RNA-seq and Spot Assays" by Taguchi et al.

### Authors

* Shodai Taguchi
* Ryosuke Matsuzawa
* Haruka Ozaki

### Contact

Please submit any questions or comments as an Issue.
