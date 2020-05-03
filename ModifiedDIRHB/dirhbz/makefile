FC = gfortran

FLAGS = -O3 -cpp -DNDEBUG -mcmodel=large -march=native

DBGFLAGS = -Og -cpp -fcheck=all -DDEBUG -mcmodel=large

BLAS = -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas
    #-lblas
    # /opt/OpenBLAS/lib/libopenblas.a \
	     

LINALG = $(BLAS) -llapack -lpthread

SRC =  dirhbz.f abjelcic.f

PAR = dirhb.par

run: $(SRC) $(PAR)
	$(FC) $(FLAGS) -o run $(SRC) $(LINALG)

dbg: $(SRC) $(PAR)
	$(FC) $(DBGFLAGS) -o dbg $(SRC) $(LINALG)

clean:
	rm -rf *.o *.wel *.del *.out *.plo run dbg
