F90 = gfortran
###CFLAGS = -Wall -g -O3 -I../include
CFLAGS = -Wall -g -O3 -I$(SUPERLU_PATH)/SRC
FFLAGS = -Wall -g -O3
###LDFLAGS = -L../lib -lsuperlu_5.0 -lblas
LDFLAGS = -L$(SUPERLU_PATH)/lib -lsuperlu_5.0 -lblas


OBJS = slu_f90.o slu.o
MODS = slu.mod

matrixtest: matrixtest.f90 $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^ $(LDFLAGS) 

slu.o: slu.f90
	$(F90) $(FFLAGS) -c -o $@ $^

clean:
	rm -f matrixtest $(OBJS) $(MODS)
