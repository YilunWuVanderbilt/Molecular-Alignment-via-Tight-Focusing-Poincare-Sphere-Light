# This is a makefile for FORTRAN program chaos
LIB = -L /usr/lib/
LIBB = -llapack -lblas

 CC = gfortran
 CFLAGS = -g -Wall 
  
  OBJS = systemparameter.o initialize.o wigner.o interaction.o diagonalize.o motion.o angle.o partition.o main.o
  EXECNAME = main
   
$(EXECNAME) : $(OBJS)
	    $(CC) -o $(EXECNAME) $(OBJS) $(LIB) $(LIBB)
		 
systemparameter.o : systemparameter.f90
	    $(CC) $(CFLAGS)-c systemparameter.f90
initialize.o : initialize.f90
	    $(CC) $(CFLAGS) -I. -c initialize.f90
wigner.o : wigner.f90
		$(CC) $(CFLAGS) -I. -c wigner.f90
interaction.o : interaction.f90
		$(CC) $(CFLAGS) -I. -c interaction.f90
diagonalize.o : diagonalize.f90
		$(CC) $(CFLAGS) -I. -c diagonalize.f90
motion.o : motion.f90
		$(CC) $(CFLAGS) -I. -c motion.f90
angle.o : angle.f90
		$(CC) $(CFLAGS) -I. -c angle.f90
partition.o : partition.f90
		$(CC) $(CFLAGS) -I. -c partition.f90
main.o : main.f90
	    $(CC) $(CFLAGS) -I. -c main.f90
		 
.PHONY : clean exe
exe :
	./$(EXECNAME)
		 
clean : 
	rm -f $(OBJS);
	rm -f $(EXECNAME) 
	rm -f *.mod
	rm -f *~
