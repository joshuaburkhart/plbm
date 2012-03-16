CXX = gcc

minimize: arraylib.o minimize.o
	$(LINK.cc) -o ./out/min arraylib.o minimize.o -lm -llapack -lblas

clean:
	rm -rf *.out *.o *.orig
