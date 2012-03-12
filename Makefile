CXX = gcc

jlapack: jlapack.o
	$(LINK.cc) -o a.jlapack.out jlapack.o

ig: jlapack.o ig.o
	$(LINK.cc) -o a.ig.out jlapack.o ig.o

clean:
	rm -rf *.out *.o
