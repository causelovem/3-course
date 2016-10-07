TARGET = matrix_block
CXX = g++
CXXFLAGS = -Wall -g -O0
FILES = plot_data/l1_data.dat plot_data/l2_data.dat plot_data/cycles_data.dat plot_data/tlb_data.dat plot_data/time_data.dat

all: $(TARGET)

report:
	cp /dev/null plot_data/l1_data.dat
	cp /dev/null plot_data/l2_data.dat
	cp /dev/null plot_data/cycles_data.dat
	cp /dev/null plot_data/tlb_data.dat
	cp /dev/null plot_data/time_data.dat
	
	for fil in $(FILES) ; do \
		echo -n 1 >> $$fil; \
	done \

	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 32 0

	for fil in $(FILES) ; do \
		echo -n 2 >> $$fil; \
	done \

	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 32 1

	for fil in $(FILES) ; do \
		echo -n 3 >> $$fil; \
	done \

	./$(TARGET) test_matrix/matrixA_500 test_matrix/matrixB_500 test_matrix/matrixC_500 -1 1

	gnuplot plot_data/l1_plot.gnu
	gnuplot plot_data/l2_plot.gnu
	gnuplot plot_data/cycles_plot.gnu
	gnuplot plot_data/tlb_plot.gnu
	gnuplot plot_data/time_plot.gnu

test:
	mkoctfile --link-stand-alone check.cpp -o check
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 -1 0
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300
	./$(TARGET) test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300 -1 1
	./check test_matrix/matrixA_300 test_matrix/matrixB_300 test_matrix/matrixC_300


$(TARGET): matrix_block.cpp
	$(CXX) $(CXXFLAGS) matrix_block.cpp -o matrix_block -lpapi

clean:
	rm -rf $(TARGET)
	rm -rf check
	rm -rf gen
	rm -rf *.o
	cp /dev/null plot_data/l1_data.dat
	cp /dev/null plot_data/l2_data.dat
	cp /dev/null plot_data/cycles_data.dat
	cp /dev/null plot_data/tlb_data.dat
	cp /dev/null plot_data/time_data.dat