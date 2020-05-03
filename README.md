# DIRHBspeedup
BLAS and LAPACK are required.

How to install OpenBLAS: just follow few simple steps described in video: https://www.youtube.com/watch?v=85hm_kbwOJs

OpenBLAS https://github.com/xianyi/OpenBLAS is an open-source efficient implementation of BLAS fine-tuned for many
modern achitectures, comparable to Intel MKL.

Keep in mind that underlying OpenBLAS is by default automatically parallelized, and if the user wants to constrain the number of used threads to single thread, one should export the following symbol in local environment: <code>OPENBLAS_NUM_THREADS=1</code>, which means, simply before running the program just type <code>export OPENBLAS_NUM_THREADS=1</code> in the console.

One can easily compare outputs of <code>ModifiedDIRHB</code> and <code>OriginalDIRHB</code> by simply using the same input files <code>dirhb.dat</code> and <code>dirhb.par</code>.
