#to make jlapack, uncomment below lines and uncomment main function in jlapack.c

CXX = gcc

#jlapack: jlapack.o
#	$(LINK.cc) -o a.jlapack.out jlapack.o

ig: asa047.o jlapack.o ig.o
	$(LINK.cc) -o a.ig.out asa047.o jlapack.o ig.o -lm -llapack -lblas

clean:
	rm -rf *.out *.o
