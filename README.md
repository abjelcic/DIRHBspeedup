# DIRHBspeedup
BLAS and LAPACK are required.

How to install OpenBLAS: just follow few simple steps described on the video https://www.youtube.com/watch?v=85hm_kbwOJs

OpenBLAS https://github.com/xianyi/OpenBLAS is an open-source efficient implementation of BLAS fine-tuned for many
modern achitectures, comparable to Intel MKL.

Keep in mind that underlying OpenBLAS is by default automatically parallelized, and if the user wants to constrain the number of used threads to single thread, one should export the following symbol in local environment: <code>export OPENBLAS_NUM_THREADS=1</code>, which means, simply before running the program just type <code>export OPENBLAS_NUM_THREADS=1</code> in the console.
