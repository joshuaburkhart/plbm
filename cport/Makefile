CC = mex

minimize: arraylib.c minimize.c
	$(CC)  minimize.c arraylib.c -f ./mexopts.sh

clean:
	rm -rf *.out *.o *.orig
