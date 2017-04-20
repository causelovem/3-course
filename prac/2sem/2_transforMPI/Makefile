TARGET = transformation
GEN = gen
CHECK = check
CXX = mpicxx
CXXFLAGS = -g -O3

all: $(TARGET)
	 
gen: $(GEN)

check: $(CHECK)

run:
	rm -rf *.out
	rm -rf *.err

	for i in 1 2 4 8 ; do \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_25 vector_res/vecB_25 1; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_26 vector_res/vecB_26 1; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_27 vector_res/vecB_27 1; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_25 vector_res/vecB_25 5; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_26 vector_res/vecB_26 5; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_27 vector_res/vecB_27 5; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_25 vector_res/vecB_25 25; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_26 vector_res/vecB_26 26; \
		mpisubmit -np $$i -w 00:30:00 ./$(TARGET) -- vectors1/vecA_27 vector_res/vecB_27 27; \
	done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

$(GEN): $(GEN).cpp
	$(CXX) $(CXXFLAGS) $(GEN).cpp -o $(GEN)

	mpisubmit -np 8 -w 00:30:00 ./$(GEN) -- vectors/vecA_25 25
	mpisubmit -np 8 -w 00:30:00 ./$(GEN) -- vectors/vecA_26 26
	mpisubmit -np 8 -w 00:30:00 ./$(GEN) -- vectors/vecA_27 27

$(CHECK): $(CHECK).cpp
	$(CXX) $(CXXFLAGS) $(CHECK).cpp -o $(CHECK)

clean:
	rm -rf $(TARGET)
	rm -rf $(GEN)
	rm -rf $(CHECK)
	rm -rf *.o