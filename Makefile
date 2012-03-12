CXX = gcc

jlapack: jlapack.o
	$(LINK.cc) -o a.jlapack.out jlapack.o

ig: ig.o jlapack.o
	$(LINK.cc) -o a.ig.out ig.o jlapack.o 

clean:
	rm -rf *.out *.o
