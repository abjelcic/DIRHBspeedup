# DIRHBspeedup

## DIRHB
This code is a small modification of the original [DIRHB code](https://www.sciencedirect.com/science/article/pii/S0010465514000836).

Primary focus is on improving the routine for spectral decomposition of
HFB matrix. HFB matrix shows many good properties and they should be utilized.
There are still many more improvements that can be done, e.g. faster calculation
of densities in coordinate space from density matrix in configuration space, or
faster calculation of single particle Hamiltonian matrix h.

The paper which proves correctness of the method is in progress.
Also, a complex version of the method is in progress, since quantum chemists
encounter similar matrices (Fock matrix).



## BLAS + LAPACK
BLAS and LAPACK are required.

How to install OpenBLAS: just follow few simple steps described in the [video](https://www.youtube.com/watch?v=85hm_kbwOJs).

[OpenBLAS](https://github.com/xianyi/OpenBLAS) is a fork of GotoBLAS, (Goto's [paper](https://dl.acm.org/doi/10.1145/1356052.1356053)) an open-source efficient implementation of BLAS fine-tuned for many
modern architectures, comparable to Intel MKL. 

Keep in mind that underlying OpenBLAS is by default automatically parallelized, and if the user wants to constrain the number of used threads to single thread, one should export the following variable in local environment: <code>OPENBLAS_NUM_THREADS=1</code>, which means, simply before running the program just type <code>export OPENBLAS_NUM_THREADS=1</code> in the console.

Of course, you can use any other BLAS implementation (a slight modification of Makefile is needed) but I recommend OpenBLAS.

For testing purpose, one can easily compare outputs of <code>ModifiedDIRHB</code> and <code>OriginalDIRHB</code> by simply using the same input files <code>dirhb.dat</code> and <code>dirhb.par</code>.



dirhbt kod, make original za originalan kod, moze ATLAS MKL...

