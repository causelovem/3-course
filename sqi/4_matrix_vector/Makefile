TARGET = matrix_vector
CXX = mpicxx
CXXFLAGS = -g -O0

all: $(TARGET)

run:
	rm -rf *.out
	rm -rf *.err
	
	for i in 1 32 64 128 256 512 ; do \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_512x512 data_B/matrixB_512 data_C/matrixC_512; \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_1024x1024 data_B/matrixB_1024 data_C/matrixC_1024; \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_2048x2048 data_B/matrixB_2048 data_C/matrixC_2048; \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_4096x4096 data_B/matrixB_4096 data_C/matrixC_4096; \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_4096x1024 data_B/matrixB_1024 data_C/matrixC_1024x1; \
		mpisubmit.bg -n $$i -w 00:02:00 ./$(TARGET) -- data_A/matrixA_1024x4096 data_B/matrixB_4096 data_C/matrixC_4096x1; \
	done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

clean:
	rm -rf $(TARGET)
	rm -rf *.o