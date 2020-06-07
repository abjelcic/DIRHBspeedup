# DIRHBspeedup
Authors: [A. Bjelcic](http://web.studenti.math.pmf.unizg.hr/~abjelcic/stranica/kontakt.html) and 
[Z. Drmac](https://web.math.pmf.unizg.hr/~drmac/).

## DIRHB
This code is a small modification of the original [DIRHB](https://www.sciencedirect.com/science/article/pii/S0010465514000836) code.

Primary focus is on improving the routine for spectral decomposition of
HFB matrix in the triaxial code. HFB matrix shows many good properties and they should be utilized.
There are still many more improvements that can be done, e.g. faster calculation
of densities in coordinate space from density matrix in configuration space or
faster calculation of single particle Hamiltonian matrix.

The paper which proves correctness of the method is in progress.
Also, a complex version of the method is in progress, since quantum chemists
encounter similar matrices (Fock matrix generated from Fock operator with Kramers partners).



## BLAS + LAPACK
BLAS and LAPACK are required.

How to install OpenBLAS on Linux machine: just follow few simple steps described in the [video](https://www.youtube.com/watch?v=85hm_kbwOJs).

How to install LAPACK on Linux machine: just follow few simple steps described in the [video](https://www.youtube.com/watch?v=5N1QYWxAAL4).


[OpenBLAS](https://github.com/xianyi/OpenBLAS) is a fork of GotoBLAS, (Goto's [paper](https://dl.acm.org/doi/10.1145/1356052.1356053)) an open-source efficient implementation of BLAS fine-tuned for many
modern architectures, comparable to Intel MKL. 

Keep in mind that underlying OpenBLAS is by default automatically parallelized, and if the user wants to constrain the number of used threads to single thread, one should export the following variable in local environment: <code>OPENBLAS_NUM_THREADS=1</code>, which means, simply before running the program just type <code>export OPENBLAS_NUM_THREADS=1</code> in the command line.

Of course, you can use any other BLAS implementation (in that case a slight modification of Makefile is needed) like [ATLAS](http://math-atlas.sourceforge.net/) or [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html), but I recommend OpenBLAS since it is open source (unlike Intel MKL) and can be ready to use within
minutes (unlike ATLAS).



## HOW TO USE
Folder <code>OriginalDIRHB</code> contains the original [DIRHB code](http://cpc.cs.qub.ac.uk/summaries/AESN_v1_0.html) for reference.

Navigate to <code>ModifiedDIRHB/dirhbz</code> or <code>ModifiedDIRHB/dirhbt</code> directory, enter the input parameters in files <code>dirhb.dat, dirhb.par</code> and type <code>make</code> followed by <code>./run</code>. If you want to use the original code for comparison purpose, type <code>make original</code> followed by <code>./original</code>.

Full paper with input parameters description can be found [here](https://github.com/abjelcic/DIRHBspeedup/blob/master/papers/DIRHB.pdf). 




## BENCHMARK - DIRHBT
Benchmark was done on Intel® NUC Kit NUC8i7HVK machine with OpenBLAS on single-thread.

Test was performed on heavy deformed <sup>240</sup>Pu nucleus with constrained parameters <code>betac=0.600</code>, <code>gammac=10.000</code>, <code>cqad=0.010</code> and <code>n0f=18</code> shells.

To reproduce the results of the original code, it took approximately 16.5 minutes.
For comparison, the original code took 35 hours.
To depict the level of agreement: the original code gave -1807.095103 MeV for total energy, while the modified code gave -1807.095038 MeV. One usually obtains the agreement within 7-8 most significant digits in total energy, which is for all practical purposes identical when calculating potential energy surfaces (PES).



## PES Tool
Requirements:
* GNU C++ Compiler (supporting at least C++17 revision of the standard)
* [Boost](https://www.boost.org/) library.

For large scale calculation of potential energy surfaces (PES) follow these steps:
* Navigate to <code>PEStool</code> folder.
* Enter input parameters into <code>pes.dat</code> file (don't disturb the original format of the file!).
* Run script: <code>bash compile.sh</code>. It will generate <code>output</code> directories system.
* Open <code>run.sh</code> script, set <code>NoSimultaneousTasks</code> and <code>OPENBLAS_NUM_THREADS</code> variables. It will allow at most <code>NoSimultaneousTasks</code> tasks to run simultaneously with allowing OpenBLAS routines to use <code>OPENBLAS_NUM_THREADS</code> threads.
* Run script: <code>nohup bash run.sh > status.out &</code>. It will run the code and dump status info into the <code>status.out</code> file.
* After the calculation is finished, run script <code>bash collect.sh</code>. It will generate (or append to) the final output file <code>pes.out</code>.
* One can use <code>bash clean.sh</code> to do the cleaning (careful, <code>output</code> folder together with the <code>pes.out</code> file will be deleted, potentially losing all calculations so far) and <code>bash killpes.sh</code> to abort all generated tasks if necessary.





