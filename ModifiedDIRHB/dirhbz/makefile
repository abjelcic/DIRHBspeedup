FC = gfortran

FLAGS = -O3 -cpp -DNDEBUG -mcmodel=large -march=native

DBGFLAGS = -g -cpp -fcheck=all -DDEBUG -mcmodel=large

ORIGFLAGS = -O3 -cpp -DORIGINAL -mcmodel=large -march=native

BLAS = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas
     # -lblas
     # /opt/OpenBLAS/lib/libopenblas.a

LINALG = $(BLAS) -llapack -lpthread

SRC =  dirhbz.f abjelcic.f

PAR = dirhb.par

run: $(SRC) $(PAR)
	$(FC) $(FLAGS) -o run $(SRC) $(LINALG)

dbg: $(SRC) $(PAR)
	$(FC) $(DBGFLAGS) -o dbg $(SRC) $(LINALG)

original: $(SRC) $(PAR)
	$(FC) $(ORIGFLAGS) -o original $(SRC) $(LINALG)

clean:
	rm -rf *.o *.wel *.del *.out *.plo run dbg original


# export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/
# export OPENBLAS_NUM_THREADS=1
