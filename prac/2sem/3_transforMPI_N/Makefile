TARGET = transformation
GEN = gen
CHECK = check
CXX = mpixlcxx_r
CXXFLAGS = -g

all: $(TARGET)
	 
gen: $(GEN)

check: $(CHECK)

run:
	rm -rf *.out
	rm -rf *.err

	#for i in 1 2 4 ; do \

	for (( count = 0; count < 20; count++ )) do \
		for j in 1 2 4 ; do \
			mpisubmit.bg -np 64 -e OMP_NUM_THREADS=$$j -w 00:15:00 -m SMP ./transformation -- vectors/vecA_23 vectors/vecB_23 0.001 $$j ; \
		done \
	done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

$(GEN): $(GEN).cpp
	$(CXX) $(CXXFLAGS) $(GEN).cpp -o $(GEN)

	mpisubmit.bg -np 16 -w 00:30:00 ./$(GEN) -- vectors/vecA_24 24
	mpisubmit.bg -np 16 -w 00:30:00 ./$(GEN) -- vectors/vecA_25 25
	mpisubmit.bg -np 16 -w 00:30:00 ./$(GEN) -- vectors/vecA_26 26
	mpisubmit.bg -np 16 -w 00:30:00 ./$(GEN) -- vectors/vecA_27 27
	mpisubmit.bg -np 16 -w 00:30:00 ./$(GEN) -- vectors/vecA_28 28

$(CHECK): $(CHECK).cpp
	$(CXX) $(CXXFLAGS) $(CHECK).cpp -o $(CHECK)

clean:
	rm -rf $(TARGET)
	rm -rf $(GEN)
	rm -rf $(CHECK)
	rm -rf *.o