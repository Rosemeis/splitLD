# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange, parallel
from libc.math cimport sqrt

##### Cython functions for splitting chromsome in windows using DP #####
##### Based on https://doi.org/10.1093/bioinformatics/btab519 by Florian Privé #####
# Convert PLINK bed format to matrix format
cpdef convertBed(signed char[:,::1] G, unsigned char[:,::1] D, int t):
	cdef int m = G.shape[0]
	cdef int n = G.shape[1]
	cdef int Bi = D.shape[1]
	cdef unsigned char[4] recode = [0, 9, 1, 2]
	cdef unsigned char mask = 3
	cdef unsigned char byte
	cdef int i, j, b, bytepart
	with nogil:
		for j in prange(m, num_threads=t):
			i = 0
			for b in range(Bi):
				byte = D[j,b]
				for bytepart in range(4):
					G[j,i] = recode[byte & mask]
					byte = byte >> 2
					i = i + 1
					if i == n:
						break

# Estimate squared correlation between variants (r^2) and compute L matrix
cpdef estimateL(signed char[:,::1] G, float[::1] F, float[::1] V, \
				float[:,::1] L, float thr, int t):
	cdef int m = G.shape[0]
	cdef int n = G.shape[1]
	cdef int W = L.shape[1]
	cdef int i, j, k, c
	cdef float cor, r2
	with nogil:
		for i in prange(m-1, num_threads=t):
			if i > (m - W):
				c = m - i - 2
			else:
				c = W - 2
			for j in range(min(i+W, m)-1, i, -1):
				cor = 0.0
				for k in range(n):
					cor = cor + (G[i,k]-F[i])*(G[j,k]-F[j])/(V[i]*V[j])
				cor = cor/<float>n
				r2 = cor*cor
				if r2 >= thr:
					L[i,c] += r2
				L[i,c] += L[i,c+1]
				c = c - 1

# Estimate E matrix used for cost estimation
cpdef estimateE(float[:,::1] L, float[:,::1] E):
	cdef int m = E.shape[0]
	cdef int W = E.shape[1]
	cdef int i, j, k
	for i in range(m-2, -1, -1):
		for j in range(W-1, -1, -1):
			if j == 0:
				E[i,j] = L[i,j]
			else:
				E[i,j] = L[i,j] + E[i+1,j-1]

# Compute cost for different number of splits
cpdef estimateC(float[:,::1] E, float[:,::1] C, int[:,::1] I, int w0, int t):
	cdef int m = E.shape[0]
	cdef int W = E.shape[1]
	cdef int K = C.shape[1]
	cdef int c, i, j, k, w
	cdef float cost
	for c in range(w0, W+1):
		C[m-c,0] = 0
		I[m-c,0] = m
	for k in range(1, K):
		with nogil:
			for i in prange(m-(k+1)*w0, -1, -1, num_threads=t):
				cost = 0.0
				for w in range(w0-1, min(W, m-i-1)):
					cost = E[i,w] + C[i+w+1,k-1]
					if cost < C[i,k]:
						C[i,k] = cost
						I[i,k] = i+w+1

# Reconstruct path of the lowest cost
cpdef reconstructPath(int[:,::1] I, int[::1] P, int k):
	cdef int i = 0
	cdef int j = k
	while j >= 0:
		i = I[i,j]
		P[j] = i
		j -= 1
