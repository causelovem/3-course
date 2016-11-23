TARGET = primes
CXX = mpicxx
CXXFLAGS = -g -O0

all: $(TARGET)

test:
	mpirun -np 4 ./$(TARGET) out 1 100000000

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

clean:
	rm -rf $(TARGET)
	rm -rf *.o