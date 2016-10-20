import sys


def explicit_metod(data, num, dist, a, step, fin):
    for i in xrange(0, fin):
        new_data = []
        new_data.extend(data)
        print(new_data[0])
        for j in xrange(1, num - 2):
            tmp = (data[j + 1] - 2 * data[j] + data[j - 1])
            new_data[j] = step * a * tmp / (dist ** 2) + data[j]
            print(new_data[j])
        print (new_data[num - 1])
        print ("\n")
        data = new_data
    return


if len(sys.argv) != 3:
    print ">Unexpected quantity of arguments, check your comand string.\n"

in_file = open(sys.argv[1])

file = in_file.readlines()

tmp = file[0].split()

num = int(tmp[0])
dist = float(tmp[1])
a = float(tmp[2])
step = float(tmp[3])
fin = int(tmp[4])

data = file[1].split()

for i in range(len(data)):
    data[i] = float(data[i])

if tmp:
    explicit_metod(data, num, dist, a, step, fin)

in_file.close()
