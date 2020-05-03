FC = gfortran

FLAGS = -O3 -cpp -DNDEBUG -mcmodel=large -march=native

DBGFLAGS = -Og -cpp -fcheck=all -DDEBUG -mcmodel=large

BLAS = # /opt/OpenBLAS/lib/libopenblas.a \
	     -lblas

LINALG = $(BLAS) -lblas -llapack -lpthread

SRC =  dirhbz.f abjelcic.f

PAR = dirhb.par

run: $(SRC) $(PAR)
	$(FC) $(FLAGS) -o run $(SRC) $(LINALG)

dbg: $(SRC) $(PAR)
	$(FC) $(DBGFLAGS) -o dbg $(SRC) $(LINALG)

clean:
	rm -rf *.o *.wel *.del *.out *.plo run dbg
