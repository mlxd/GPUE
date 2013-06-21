OS:=	$(shell uname)
ifeq ($(OS),Darwin)
CUDA_LIB	= /usr/local/cuda/lib
CFLAGS		= -g -G -m64
else
CUDA_LIB	= /usr/local/cuda/lib64
CFLAGS		= -g -G
endif

CC		= nvcc#-save-temps
CLINKER		= $(CC)
RM		= /bin/rm
CUDA_HEADER	= /usr/local/cuda/include
INCFLAGS	= -I$(CUDA_HEADER)
CFLAGS		= -g -G -m64
LDFLAGS		= -L$(CUDA_LIB) 
EXECS		= gpue # BINARY NAME HERE

gpue: fileIO.o kernels.o split_op.o
	gcc fileIO.o split_op.o kernels.o  $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -lm -lcufft -lcudart -o gpue
	rm -rf ./*.o

split_op.o: ./src/split_op.cu ./include/split_op.h ./include/kernels.h ./include/constants.h ./include/fileIO.h Makefile
	$(CC) -c ./src/split_op.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -arch=sm_30
	
kernels.o: ./include/split_op.h Makefile ./include/constants.h ./include/kernels.h ./src/kernels.cu
	$(CC) -c ./src/kernels.cu -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -arch=sm_30
	
fileIO.o: ./include/fileIO.h ./src/fileIO.cc Makefile
	gcc -c ./src/fileIO.cc -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS)
	
default:	gpue
all:		gpue

.c.o:
	$(CC) $(INCFLAGS) $(CFLAGS) -c $<

clean:
	@-$(RM) -f r_0 Phi_0 E* px_* py_0* xPy* xpy* ypx* x_* y_* yPx* p0* p1* p2* EKp* EVr* gpot wfc* Tpot 0* V_* K_* Vi_* Ki_* 0i* k s_* si_* *.o *~ PI* $(EXECS) $(OTHER_EXECS) *.dat *.png *.eps *.ii *.i *cudafe* *fatbin* *hash* *module* *ptx; 
