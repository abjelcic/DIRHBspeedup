# DIRHBspeedup

## DIRHB
This code is a small modification of the original [DIRHB code](https://www.sciencedirect.com/science/article/pii/S0010465514000836).

Primary focus is on improving the routine for spectral decomposition of
HFB matrix. HFB matrix shows many good properties and they should be utilized.
There are still many more improvements that can be done, e.g. faster calculation
of densities in coordinate space from density matrix in configuration space,
faster calculation of single particle Hamiltonian matrix h and improving 
routines for calculating transformation to canonical basis and center of
mass correction.

The paper which proves correctness of the method is in progress.
Also, a complex version of the method is in progress, since quantum chemists
encounter similar matrices (Fock matrix).



## BLAS + LAPACK
BLAS and LAPACK are required.

How to install OpenBLAS on Linux machine: just follow few simple steps described in the [video](https://www.youtube.com/watch?v=85hm_kbwOJs).

[OpenBLAS](https://github.com/xianyi/OpenBLAS) is a fork of GotoBLAS, (Goto's [paper](https://dl.acm.org/doi/10.1145/1356052.1356053)) an open-source efficient implementation of BLAS fine-tuned for many
modern architectures, comparable to Intel MKL. 

Keep in mind that underlying OpenBLAS is by default automatically parallelized, and if the user wants to constrain the number of used threads to single thread, one should export the following variable in local environment: <code>OPENBLAS_NUM_THREADS=1</code>, which means, simply before running the program just type <code>export OPENBLAS_NUM_THREADS=1</code> in the command line.

Of course, you can use any other BLAS implementation (in that case a slight modification of Makefile is needed) like [ATLAS](http://math-atlas.sourceforge.net/) or [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html), but I recommend OpenBLAS since it is open source (unlike Intel MKL) and can be ready to use within
minutes (unlike ATLAS).


## HOW TO USE
Folder <code>OriginalDIRHB</code> contains origianal DIRHB code for reference.

Navigate to <code>ModifiedDIRHB/dirhbz</code> or <code>ModifiedDIRHB/dirhbt</code> directory, enter the input parameters in files <code>dirhb.dat, dirhb.par</code> and type <code>make</code> followed by <code>./run</code>. If you want to use the origianl code for comparison purpose, type <code>make original</code> followed by <code>./original</code>.


## BENCHMARK - DIRHBT
Benchmark was done on Intel® NUC Kit NUC8i7HVK machine with OpenBLAS on single-thread.

Test was performed on heavy deformed <sup>240</sup>Pu nucleus with constrained parameters <code>betac=0.600</code>, <code>gammac=10.000</code>, <code>cqad=0.010</code> and <code>n0f=18</code> shells. To reproduce the results of the
original code, it took 8 seconds in preparation phase, 16.5 minutes in iteration phase, 3.1 minutes in transformation to canonical basis and 92.4 minutes in center of mass correction calculation. For comparison, iteration phase of the original code took approximately 35 hours.

For calculating PES, one can use simpler variant of center of mass correction. Since many applications require at least 20 shells, this code makes it feasible to obtain PES of very heavy deformed nuclei in reasonable time.

In future updates, I will try to improve the transformation to canonical basis and center of mass correction calculations.
