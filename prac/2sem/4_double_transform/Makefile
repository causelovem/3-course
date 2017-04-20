TARGET = transformation
GEN = gen
CHECK = check
CXX = mpixlcxx_r
CXXFLAGS = -g -qsmp=omp

all: $(TARGET)

run:
	rm -rf *.out
	rm -rf *.err

	for i in 1 2 4 ; do \
		mpisubmit.bg -np $$i -e "OMP_NUM_THREADS=1" -w 00:45:00 -m SMP ./transformation -- vectors/vecA_23 vectors/vecB_23 1 2; \
		mpisubmit.bg -np $$i -e "OMP_NUM_THREADS=2" -w 00:45:00 -m SMP ./transformation -- vectors/vecA_23 vectors/vecB_23 1 2; \
		mpisubmit.bg -np $$i -e "OMP_NUM_THREADS=4" -w 00:45:00 -m SMP ./transformation -- vectors/vecA_23 vectors/vecB_23 1 2; \
	done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

$(GEN): $(GEN).cpp
	$(CXX) $(CXXFLAGS) $(GEN).cpp -o $(GEN)

	mpisubmit.bg -np 8 -w 00:30:00 ./$(GEN) -- vectors/vecA_23 23

$(CHECK): $(CHECK).cpp
	$(CXX) $(CXXFLAGS) $(CHECK).cpp -o $(CHECK)

clean:
	rm -rf $(TARGET)
	rm -rf $(GEN)
	rm -rf $(CHECK)
	rm -rf *.o