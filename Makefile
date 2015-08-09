OS:=	$(shell uname)
ifeq ($(OS),Darwin)
CC		= nvcc -ccbin /usr/bin/clang --ptxas-options=-v#-save-temps
CUDA_LIB	= /usr/local/cuda/lib
CUDA_HEADER	= /usr/local/cuda/include
CFLAGS		= -g -O3 -m64 -std=c++11
GPU_ARCH	= sm_30
else
CC			= /usr/local/cuda/bin/nvcc --ptxas-options=-v --compiler-options -Wall #-save-temps
CUDA_LIB	= /usr/local/cuda/lib64
CUDA_HEADER	= /usr/local/cuda/include
CHOSTFLAGS	= #-fopenmp
CFLAGS		= -g -std=c++11#-malign-double
GPU_ARCH	= sm_30
endif

CLINKER		= $(CC)
RM		= /bin/rm
INCFLAGS	= -I$(CUDA_HEADER) -I/opt/builds/include 
LDFLAGS		= -L$(CUDA_LIB) -L/opt/builds/lib 
EXECS		= gpue # BINARY NAME HERE

gpue: fileIO.o kernels.o split_op.o tracker.o minions.o ds.o 
#node.o edge.o lattice.o
	g++ *.o $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS) -lm -lcufft -lcudart -o gpue 
	#rm -rf ./*.o

split_op.o: ./src/split_op.cu ./include/split_op.h ./include/kernels.h ./include/constants.h ./include/fileIO.h ./include/minions.h Makefile
	$(CC) -c  ./src/split_op.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH) 
	
kernels.o: ./include/split_op.h Makefile ./include/constants.h ./include/kernels.h ./src/kernels.cu
	$(CC) -c  ./src/kernels.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -arch=$(GPU_ARCH) 
	
fileIO.o: ./include/fileIO.h ./src/fileIO.cc Makefile
	g++ -c ./src/fileIO.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) 

tracker.o: ./src/tracker.cc ./include/tracker.h ./include/fileIO.h
	g++ -c ./src/tracker.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

minions.o: ./src/minions.cc ./include/minions.h
	g++ -c ./src/minions.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS) 

ds.o: ./src/ds.cc ./include/ds.h
	g++ -c ./src/ds.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

node.o: ./src/node.cc ./include/node.h
	$(CC) -c ./src/node.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

edge.o: ./src/edge.cc ./include/edge.h
	$(CC) -c ./src/edge.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

lattice.o: ./src/lattice.cc ./include/lattice.h
	$(CC) -c ./src/lattice.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

graphtest.o: ./src/graphtest.cc ./include/latticegraph.hh
	$(CC) -c ./src/graphtest.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

gtest:  edge.o node.o lattice.o graphtest.o
	$(CC) $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS) edge.o node.o lattice.o graphtest.o -o gtest

minions: ./src/minions.cc ./include/minions.h minions.o
	$(CC) minions.o -o mintest $(INCFLAGS) $(CFLAGS) $(LDFLAGS) 

tracker_test: tracker.o fileIO.o ./src/tracker.cc ./include/fileIO.h ./src/fileIO.cc ./include/tracker.h
	gcc ./tracker.o ./fileIO.o -o tracker_test $(INCFLAGS) $(CFLAGS) $(LDFLAGS)
	
default:	gpue
all:		gpue test

.c.o:
	$(CC) $(INCFLAGS) $(CFLAGS) -c $<

clean:
	@-$(RM) -f r_0 Phi_0 E* px_* py_0* xPy* xpy* ypx* x_* y_* yPx* p0* p1* p2* EKp* EVr* gpot wfc* Tpot 0* V_* K_* Vi_* Ki_* 0i* k s_* si_* *.o *~ PI* $(EXECS) $(OTHER_EXECS) *.dat *.png *.eps *.ii *.i *cudafe* *fatbin* *hash* *module* *ptx test* vort* v_opt*; 
