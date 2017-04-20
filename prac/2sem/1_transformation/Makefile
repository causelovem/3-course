TARGET = transformation
GEN = gen
CHECK = check
CXX = g++
CXXFLAGS = -g -O3 -fopenmp

all: $(TARGET)
	 
gen: $(GEN)

check: $(CHECK)

run:
	rm -rf *.out
	rm -rf *.err

	for i in 1 2 4 8 ; do \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_20 vecB_20 1 $$i; \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_24 vecB_24 1 $$i; \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_20 vecB_20 5 $$i; \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_24 vecB_24 5 $$i; \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_20 vecB_20 20 $$i; \
		ompsubmit -n $$i -w 00:30:00 ./$(TARGET) -- vectors/vecA_24 vecB_24 24 $$i; \
	done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

$(GEN): $(GEN).cpp
	$(CXX) $(CXXFLAGS) $(GEN).cpp -o $(GEN)

	./$(GEN) vectors/vecA_20 20 10
	./$(GEN) vectors/vecA_24 24 10
	./$(GEN) vectors/vecA_28 28 10

$(CHECK): $(CHECK).cpp
	$(CXX) $(CXXFLAGS) $(CHECK).cpp -o $(CHECK)

clean:
	rm -rf $(TARGET)
	rm -rf $(GEN)
	rm -rf $(CHECK)
	rm -rf *.o