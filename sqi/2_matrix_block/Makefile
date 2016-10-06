TARGET = matrix_block
CXX = g++
CXXFLAGS = -Wall -g -O0

all: $(TARGET)

report:
	#cp /dev/null plot.dat
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 0
	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 1
	#gnuplot plot.gnu

test:
	mkoctfile --link-stand-alone check.cpp -o check
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 30 0
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 30 1
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300


$(TARGET): matrix_block.cpp
	$(CXX) $(CXXFLAGS) matrix_block.cpp -o matrix_block -lpapi

clean:
	rm -rf $(TARGET)
	rm -rf check
	cp /dev/null plot.dat
