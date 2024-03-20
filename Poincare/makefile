# This is a makefile for FORTRAN program Poincare
LIB = -L /usr/lib/
LIBB = -llapack -lblas

CC = gfortran
CFLAGS = -g -Wall 
OBJS = systemparameter.o integral.o  poincare.o
EXECNAME = poincare
  
$(EXECNAME) : $(OBJS)
		$(CC) -o $(EXECNAME) $(OBJS) $(LIB) $(LIBB)
			
systemparameter.o : systemparameter.f90
		$(CC) $(CFLAGS)-c systemparameter.f90
integral.o : integral.f90
		$(CC) $(CFLAGS) -I. -c integral.f90
poincare.o : poincare.f90
		$(CC) $(CFLAGS) -I. -c poincare.f90
			 
.PHONY : clean exe
exe :
		./$(EXECNAME)
clean :
		rm -f $(OBJS);
		rm -f $(EXECNAME)
		rm -f *.mod
		rm -f *~

