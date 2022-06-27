"""
Optimal window splitting based on LD.
Fast Cython implementation of 'Optimal linkage disequilibrium splitting' by Florian Privé.
https://doi.org/10.1093/bioinformatics/btab519
"""

__author__ = "Jonas Meisner"

# Libraries
import argparse
import os
import subprocess
import sys
import numpy as np
import allel
from datetime import datetime
from math import ceil

# Import own script
from splitLD import shared_cy

# Reader help function
def extract_length(filename):
	process = subprocess.Popen(['wc', '-l', filename], stdout=subprocess.PIPE)
	result, err = process.communicate()
	return int(result.split()[0])

### Argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf",
	help="Genotype file of single chromosome in VCF format")
parser.add_argument("-p", "--plink",
	help="Prefix for binary PLINK files")
parser.add_argument("-t", "--threads", type=int, default=1,
	help="Number of threads")
parser.add_argument("-o", "--out", default="windows",
	help="Prefix for output files")
parser.add_argument("--min_length", type=int, default=500,
	help="Minimum number of SNPs in windows")
parser.add_argument("--max_length", type=int, default=5000,
	help="Maximum number of SNPs in windows")
parser.add_argument("--threshold", type=float, default=0.1,
	help="r2 threshold to be included")

##### Optimal window splitting #####
def main():
	args = parser.parse_args()
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()
	print("Optimal LD splitting")
	print("Please cite original paper: https://doi.org/10.1093/bioinformatics/btab519")
	assert args.min_length <= args.max_length, "Min length > max length!"
	assert (args.vcf is not None) or (args.plink is not None), \
		"Please provide genotype file!"
	
	### Load genotype matrix
	# Read VCF file
	if args.vcf is not None:
		print("\rLoading VCF file...", end="")
		vcf = allel.read_vcf(args.vcf)
		G = vcf['calldata/GT'].sum(axis=2, dtype=np.int8)
		m, n = G.shape
		del vcf

	# Read PLINK files
	if args.plink is not None:
		print("Reading in data matrix from PLINK files.")
		# Finding length of .fam and .bim file and read .bed file into NumPy array
		n = extract_length(args.plink + ".fam")
		m = extract_length(args.plink + ".bim")
		with open(args.plink + ".bed", "rb") as bed:
			B = np.fromfile(bed, dtype=np.uint8, offset=3)
		Bi = ceil(n/4) # Length of bytes to describe n individuals
		D = B.reshape((m, Bi))
		G = np.zeros((m, n), dtype=np.int8)
		shared_cy.convertBed(G, D, args.threads)
		del B, D
	print("\rLoaded {} samples and {} variants.".format(n, m))

	# Setup parameters
	assert args.max_length <= m, "Max length > chromosome!"
	maxW = ceil(m/args.min_length) + 1
	F = np.mean(G, axis=1, dtype=np.float32)
	V = np.std(G, axis=1, dtype=np.float32)

	# Estimating L matrix
	print("Estimating correlations and L matrix...")
	L = np.zeros((m, args.max_length), dtype=np.float32)
	shared_cy.estimateL(G, F, V, L, args.threshold, args.threads)
	del G, F, V

	# Estimating E matrix
	print("Estimating E matrix...")
	E = np.zeros((m, args.max_length), dtype=np.float32)
	shared_cy.estimateE(L, E)
	del L

	# Compute optimal path of splits
	print("Estimating cost of paths...")
	C = np.zeros((m, maxW), dtype=np.float32)
	I = np.zeros((m, maxW), dtype=np.int32)
	C.fill(np.inf)
	I.fill(-1)
	shared_cy.estimateC(E, C, I, args.min_length, args.threads)
	del E

	# Reconstruct most optimal path
	P = np.zeros(maxW, dtype=np.int32)
	minC = np.min(C[0,1:])
	for i in range(1, maxW):
		if abs(C[0,i] - minC) < 1e-6:
			optK = i
	del C
	shared_cy.reconstructPath(I, P, optK)
	P = P[:(P.shape[0]-(np.sum(P == 0)-1))]
	if P[-1] == -1:
		P = np.array([m,0], dtype=np.int32)
	print("{} optimal blocks".format(P.shape[0]-1))

	# Save matrices
	np.savetxt(args.out + ".win.txt", P[::-1], fmt="%i")
	print("Saved optimal window indices in " + args.out + ".win.txt")
