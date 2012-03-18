DEFAULT_CXX = g++
MATLAB_CXX = mex

default: arraylib.cpp minimize.cpp
	$(DEFAULT_CXX) minimize.cpp arraylib.cpp -lm -llapack -lblas

matlab: arraylib.cpp minimize.cpp
	bash -c "module load matlab && \
	$(MATLAB_CXX)  minimize.cpp arraylib.cpp -f ./mexopts.sh"

clean:
	rm -rf *.out *.o *.orig
