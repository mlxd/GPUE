CUDA_HOME = /opt/cuda/
GPU_ARCH	= sm_30
OS:=	$(shell uname)
ifeq ($(OS),Darwin)
CUDA_LIB	= $(CUDA_HOME)/lib
CUDA_HEADER	= $(CUDA_HOME)/include
CC		= $(CUDA_HOME)/bin/nvcc -ccbin /usr/bin/clang --ptxas-options=-v#-save-temps
CFLAGS		= -g -std=c++11
else
CUDA_LIB	= $(CUDA_HOME)/lib64
CUDA_HEADER	= $(CUDA_HOME)/include
CC		= $(CUDA_HOME)/bin/nvcc --ptxas-options=-v --compiler-options -Wall #-save-temps
CHOSTFLAGS	= #-fopenmp
CFLAGS		= -g -std=c++11 -Xcompiler '-std=c++11' -Xcompiler '-fopenmp' #-malign-double -lboost_math
endif

CLINKER		= $(CC) 
RM		= /bin/rm
INCFLAGS	= -I$(CUDA_HEADER) 
LDFLAGS		= -L$(CUDA_LIB) 
EXECS		= gpue # BINARY NAME HERE

gpue: fileIO.o kernels.o split_op.o tracker.o minions.o ds.o edge.o node.o lattice.o manip.o vort.o parser.o evolution.o init.o unit_test.o operators.o
	$(CC) *.o $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS) -lm -lcufft -lcudart -o gpue
	#rm -rf ./*.o

init.o: ./src/init.cu ./include/split_op.h ./include/kernels.h ./include/constants.h ./include/fileIO.h ./include/minions.h ./include/parser.h ./include/evolution.h Makefile
	$(CC) -c  ./src/init.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH)

split_op.o: ./src/split_op.cu ./include/split_op.h ./include/kernels.h ./include/constants.h ./include/fileIO.h ./include/minions.h
	$(CC) -c  ./src/split_op.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH)

kernels.o: ./include/split_op.h Makefile ./include/constants.h ./include/kernels.h ./src/kernels.cu
	$(CC) -c  ./src/kernels.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -arch=$(GPU_ARCH)

fileIO.o: ./include/fileIO.h ./src/fileIO.cc Makefile
	$(CC) -c ./src/fileIO.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS)

tracker.o: ./src/tracker.cc ./include/tracker.h ./include/fileIO.h
	$(CC) -c ./src/tracker.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

operators.o: ./src/operators.cc ./include/operators.h
	$(CC) -c ./src/operators.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

minions.o: ./src/minions.cc ./include/minions.h
	$(CC) -c ./src/minions.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

parser.o: ./src/parser.cc ./include/parser.h
	 $(CC) -c ./src/parser.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

unit_test.o: ./src/unit_test.cu ./include/unit_test.h
	$(CC) -c ./src/unit_test.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

evolution.o: ./src/evolution.cu ./include/evolution.h ./include/split_op.h ./include/constants.h ./include/kernels.h ./include/fileIO.h 
	$(CC) -c ./src/evolution.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH)

ds.o: ./src/ds.cc ./include/ds.h ./include/operators.h
	$(CC) -c ./src/ds.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

node.o: ./src/node.cc ./include/node.h
	$(CC) -c ./src/node.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

edge.o: ./src/edge.cc ./include/edge.h
	$(CC) -c ./src/edge.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

lattice.o: ./src/lattice.cc ./include/lattice.h
	$(CC) -c ./src/lattice.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

manip.o: ./src/manip.cu ./include/manip.h
	$(CC) -c ./src/manip.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

vort.o: ./src/vort.cc ./include/vort.h
	$(CC) -c ./src/vort.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

graphtest.o: ./src/graphtest.cc
	$(CC) -c ./src/graphtest.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS)

gtest:  edge.o node.o lattice.o graphtest.o
	$(CC) $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CHOSTFLAGS) edge.o node.o lattice.o graphtest.o -o gtest

minions: ./src/minions.cc ./include/minions.h minions.o
	$(CC) minions.o -o mintest $(INCFLAGS) $(CFLAGS) $(LDFLAGS)

tracker_test: tracker.o fileIO.o ./src/tracker.cc ./include/fileIO.h ./src/fileIO.cc ./include/tracker.h
	$(CC) ./tracker.o ./fileIO.o -o tracker_test $(INCFLAGS) $(CFLAGS) $(LDFLAGS)

default:	gpue
all:		gpue test

.c.o:
	$(CC) $(INCFLAGS) $(CFLAGS) -c $<

clean:
	@-$(RM) -f r_0 Phi_0 E* px_* py_0* xPy* xpy* ypx* x_* y_* yPx* p0* p1* p2* EKp* EVr* gpot wfc* Tpot 0* V_* K_* Vi_* Ki_* 0i* k s_* si_* *.o *~ PI* $(EXECS) $(OTHER_EXECS) *.dat *.png *.eps *.ii *.i *cudafe* *fatbin* *hash* *module* *ptx test* vort* v_opt*;
