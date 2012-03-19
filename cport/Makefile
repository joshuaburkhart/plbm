DEFAULT_CXX = g++
MATLAB_CXX = mex

default: glob_init.cpp arraylib.cpp minimize.cpp main.cpp
	$(DEFAULT_CXX) main.cpp minimize.cpp arraylib.cpp glob_init.cpp -lm -llapack -lblas

matlab: mex_gateway.cpp arraylib.cpp minimize.cpp
	bash -c "module load matlab && \
	$(MATLAB_CXX) minimize.cpp mex_gateway.cpp arraylib.cpp -f ./mexopts.sh"

clean:
	rm -rf *.out *.o *.orig *.mexa64
