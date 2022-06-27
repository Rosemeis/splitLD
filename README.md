# Optimal linkage disequilibrium splitting
Reimplementation in Python/Cython of the algorithm presented by Florian Privé in [Optimal linkage disequilibrium splitting](https://doi.org/10.1093/bioinformatics/btab519). The reimplementation is a very fast multithreadeded version and it works a standalone Python command-line program. It takes both binary PLINK files and VCF files as genotype data input and outputs the path of the optimal window splits.

### Dependencies
The program relies on the following Python packages that you can install through conda (recommended) or pip:
- numpy
- cython
- scikit-allel

You can create an environment through conda easily as follows:
```
conda env create -f environment.yml
```

## Install and build
```bash
git clone https://github.com/Rosemeis/splitLD.git
cd splitLD
python setup.py build_ext --inplace
pip3 install -e .
```

You can now run the program with the `splitLD` command.

## Usage
```bash
# VCF
splitLD --vcf data.chr1.vcf.gz --out split.chr1 --threads 32

# Binary PLINK
splitLD --plink data.chr1 --out split.chr1 --threads 32

# Outputs indices for optimal path window splits (split.chr1.win.txt) 
```

See all available options in the program as follows:
```bash
splitLD -h
```

The output file with window splits can be used directly as input to [HaploNet](https://github.com/Rosemeis/HaploNet):
```bash
haplonet --vcf data.chr1.vcf.gz --out haplonet.chr1 --windows split.chr1.win.txt
```
