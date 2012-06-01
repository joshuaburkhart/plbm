# distributed under the terms of the GNU General Public License
# Copyright 2012 Joshua Burkhart

DEFAULT_CXX = g++
CUDA_CXX = nvcc
MATLAB_CXX = mex
MAC_CXX = /opt/local/bin/g++-mp-4.7

default: glob_init.cpp arraylib.cpp minimize.cpp main.cpp
	bash -c "\
	module load gcc & \
	$(DEFAULT_CXX) -O3 main.cpp minimize.cpp arraylib.cpp glob_init.cpp -lm -llapack -lblas -fopenmp"

mac: glob_init.cpp arraylib.cpp minimize.cpp main.cpp
	bash -c "\
	module load gcc & \
	$(MAC_CXX) -O3 main.cpp minimize.cpp arraylib.cpp glob_init.cpp -lm -llapack -lblas -fopenmp"

cuda: glob_init.cpp arraylib.cu minimize.cpp main.cpp
	bash -c "\
	module load gcc/4.3.6 && \
	module load cuda && \
	$(CUDA_CXX) -O3 main.cpp minimize.cpp arraylib.cu glob_init.cpp -lm -llapack -lblas"

matlab: mex_gateway.cpp arraylib.cpp minimize.cpp
	bash -c "\
	module load matlab && \
	$(MATLAB_CXX) mex_gateway.cpp minimize.cpp arraylib.cpp -f ./mexopts.sh"

test: glob_init.cpp arraylib.cpp minimize.cpp main.cpp
	bash -c "\
	module load gcc & \
	$(DEFAULT_CXX) test_main.cpp minimize.cpp arraylib.cpp glob_init.cpp -lm -llapack -lblas -fopenmp"

test_mac: glob_init.cpp arraylib.cpp minimize.cpp main.cpp
	bash -c "\
	module load gcc & \
	$(MAC_CXX) test_main.cpp minimize.cpp arraylib.cpp glob_init.cpp -lm -llapack -lblas -fopenmp"

test_cuda: glob_init.cpp arraylib.cu minimize.cpp main.cpp
	bash -c "\
	module load gcc/4.3.6 && \
	module load cuda && \
	$(CUDA_CXX) test_main.cpp minimize.cpp arraylib.cu glob_init.cpp -lm -llapack -lblas"

clean:
	rm -rf *.out *.o *.orig *.mexa64
