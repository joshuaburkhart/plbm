#to make jlapack, uncomment below lines and uncomment main function in jlapack.c

CXX = gcc

#jlapack: jlapack.o
#	$(LINK.cc) -o a.jlapack.out jlapack.o

ig: jlapack.o ig.o
	$(LINK.cc) -o a.ig.out jlapack.o ig.o -lm -llapack -lblas

clean:
	rm -rf *.out *.o
