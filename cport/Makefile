CXX = mex

minimize: arraylib.cpp minimize.cpp
	bash -c "module load matlab && \
	$(CXX)  minimize.cpp arraylib.cpp -f ./mexopts.sh"

clean:
	rm -rf *.out *.o *.orig
