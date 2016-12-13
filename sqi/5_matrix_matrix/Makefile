TARGET = matrix_matrix
CXX = mpicxx
CXXFLAGS = -g -O5

all: $(TARGET)

run:
	rm -rf *.out
	rm -rf *.err

	mpisubmit.bg -n 1 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_1024 data_B/matrixB_1024 data_C/matrixC_1024
	mpisubmit.bg -n 1 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_2048 data_B/matrixB_2048 data_C/matrixC_2048
	mpisubmit.bg -n 1 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_4096 data_B/matrixB_4096 data_C/matrixC_4096
	mpisubmit.bg -n 1 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_8192 data_B/matrixB_8192 data_C/matrixC_8192

	mpisubmit.bg -n 8 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_1024 data_B/matrixB_1024 data_C/matrixC_1024
	mpisubmit.bg -n 8 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_2048 data_B/matrixB_2048 data_C/matrixC_2048
	mpisubmit.bg -n 8 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_4096 data_B/matrixB_4096 data_C/matrixC_4096
	mpisubmit.bg -n 8 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_8192 data_B/matrixB_8192 data_C/matrixC_8192

	mpisubmit.bg -n 64 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_1024 data_B/matrixB_1024 data_C/matrixC_1024
	mpisubmit.bg -n 64 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_2048 data_B/matrixB_2048 data_C/matrixC_2048
	mpisubmit.bg -n 64 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_4096 data_B/matrixB_4096 data_C/matrixC_4096
	mpisubmit.bg -n 64 -w 02:00:00 ./$(TARGET) -- data_A/matrixA_8192 data_B/matrixB_8192 data_C/matrixC_8192

	mpisubmit.bg -n 125 -w 00:15:00 ./$(TARGET) -- data_A/matrixA_1024 data_B/matrixB_1024 data_C/matrixC_1024
	mpisubmit.bg -n 125 -w 00:15:00 ./$(TARGET) -- data_A/matrixA_2048 data_B/matrixB_2048 data_C/matrixC_2048
	mpisubmit.bg -n 125 -w 00:15:00 ./$(TARGET) -- data_A/matrixA_4096 data_B/matrixB_4096 data_C/matrixC_4096
	mpisubmit.bg -n 125 -w 00:15:00 ./$(TARGET) -- data_A/matrixA_8192 data_B/matrixB_8192 data_C/matrixC_8192

	llsubmit prog1024.jcf
	llsubmit prog2048.jcf
	llsubmit prog4096.jcf
	llsubmit prog8192.jcf

	llq -W | grep $USER
	
	# for i in 1 8 64 125 ; do \
		# mpisubmit.bg -n $$i -w 00:15:00 ./$(TARGET) -- data_A/matrixA_1024 data_B/matrixB_1024 data_C/matrixC_1024; \
		# mpisubmit.bg -n $$i -w 00:15:00 ./$(TARGET) -- data_A/matrixA_2048 data_B/matrixB_2048 data_C/matrixC_2048; \
		# mpisubmit.bg -n $$i -w 00:15:00 ./$(TARGET) -- data_A/matrixA_4096 data_B/matrixB_4096 data_C/matrixC_4096; \
		# mpisubmit.bg -n $$i -w 00:15:00 ./$(TARGET) -- data_A/matrixA_8192 data_B/matrixB_8192 data_C/matrixC_8192; \
	# done \

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)

clean:
	rm -rf $(TARGET)
	rm -rf *.o
	rm -rf src/*.o