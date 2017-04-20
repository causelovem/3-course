import sys


file_out = open("/home/causelove/mishenka/prog/sqi-solutions/3_transforMPI_N/23_res/23_resfile", 'a')

mse = 0

for j in range(1, 60):
    file_in = open(sys.argv[j], 'r')

    for i in range(0, 6):
        tmp = file_in.readline().split()

    res = tmp[4].split(",")
    lol = res[0].split("(")

    file_out.write(str(j) + ' ' + lol[1] + '\n')

    mse += float(lol[1])

    file_in.close()

mse = mse / 59

print mse

file_out.close()
