# Investigating the Effects of a Liquid-Handling Robot - Opentrons OT-2 Parameter Variations on Yeast Growth and Gene Expression Using Growth Assays and RNA-seq

This repository provides the code, data, and analysis pipelines used in the study "Investigating the Effects of Robotic OT-2 Parameter Variations on Yeast Growth and Gene Expression: A Comprehensive Study Using RNA-seq and Spot Assays" by Taguchi et al.

## Overview

This study explores how variations in pipetting speed using the Opentrons OT-2 liquid-handling robot affect yeast growth and gene expression. The analysis integrates RNA-seq data and quantitative growth assays to provide insights into the relationship between robotic parameters and biological outcomes.

## Features

- **OT-2 Protocols**: Automated experimental protocols for yeast handling and RNA-seq preparation using the Opentrons OT-2 robot.
- **RNA-seq Analysis**: Preprocessing, alignment, and differential expression analysis pipelines for RNA-seq data.
- **BaQFA Analysis**: Quantitative growth assay of yeast colonies using BaQFA (Barcode Quantitative Fitness Analysis).
- **Integrated Analysis**: Comparison of yeast growth and gene expression under different pipetting speeds.


## Requirements

- **Python**: 3.10
- **Libraries**: Listed in `requirements.txt` and `requirements.yaml`, including:
  - `pandas`, `numpy` for data manipulation
  - `STAR`, `RSEM`, `fastp` for RNA-seq preprocessing
  - `baQFA` for growth assay analysis
- **R**: Required for RNA-seq and BaQFA data visualization and statistical analysis.

## Usage

### OT-2 Protocols

1. **Prepare Environment**:
```bash
   mamba create --name labopt --file ./requirements.yaml
   conda activate labopt
```

2. Configure .env File: Create a .env file in the experiment directory with the following: (Make sure that .env file is placed in the same directory as `run_experiment.ipynb`)

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

3. Run Protocols: Use the Python scripts in `code/OT2_protocols/` to run the OT-2 protocols. 

## RNA-seq analysis

1. Run the bash scripts step by step in `code/RNAseq_analysis/analysis/` to analyze the RNA-seq data.

2. RNA-seq analysis results are output to the `results/240328_RNAseq2nd` directory.

## BaQFA analysis
1. Prepare the time-series image of yeast colonies obtained by the flatbed scanner.

2. Run the Jupyter Notebooks `code/BaQFA_protocols/001_rotate-crop.ipynb` and `code/BaQFA_protocols/002_run-baqfa.ipynb` to perform the quantitative growth assay of yeast colonies using BaQFA.

3. BaQFA analysis results are output to the `results/240415_baQFA` directory.

## Integrated Analysis
- Combine RNA-seq and BaQFA results for integrated insights. Use R scripts in results/240415_baQFA/integrated_analysis/.

## Notes
- This repository provides the tools and data to reproduce the experiments but does not guarantee identical results.
- Refer to the official Opentrons documentation for OT-2 operation.

## License
This repository is released under the MIT License.

## Acknowledgements
This work is based on the study by Taguchi et al., "Investigating the Effects of Robotic OT-2 Parameter Variations on Yeast Growth and Gene Expression: A Comprehensive Study Using RNA-seq and Spot Assays."

## Authors
Shodai Taguchi
Ryosuke Matsuzawa
Haruka Ozaki

## Contact
For questions or comments, please submit an issue in this repository.
