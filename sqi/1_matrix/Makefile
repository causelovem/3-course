TARGET = matrix
CXX = g++
CXXFLAGS = -Wall -g -O0

all: $(TARGET)

report:
	cp /dev/null plot.dat
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 1
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 2
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 3
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 4
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 5
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 6
	gnuplot plot.gnu

test:
	mkoctfile --link-stand-alone check.cpp -o check
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 1
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 2
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 3
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 4
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 5
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 6


$(TARGET): matrix.cpp
	$(CXX) $(CXXFLAGS) matrix.cpp -o matrix

clean:
	rm -rf $(TARGET)
	cp /dev/null plot.dat
