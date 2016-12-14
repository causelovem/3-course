from sys import stdout


def solve(st, num, matrix):
    pos = [-1] * num
    weight = [float("inf")] * num
    checked = [0] * num
    weight[st] = 0
    pos[st] = st
    while True:
        min = float("inf")
        min_pos = -1
        for i in range(num):
            if (min > weight[i]) and (not checked[i]):
                min_pos = i
                min = weight[i]

        if min_pos == -1:
            break

        for i in range(num):
            if (not checked[i]) and (matrix[min_pos][i] != 0) and (matrix[min_pos][i] + weight[min_pos] < weight[i]):
                weight[i] = matrix[min_pos][i] + weight[min_pos]
                pos[i] = min_pos
        checked[min_pos] = True

    return (weight, pos)


file = open("map", 'r')
tmp = file.readline().split()

# number
num = int(tmp[0])
# start
st = int(tmp[1])
# finish
fin = int(tmp[2])

matrix = [[0] * num for i in range(num)]
for i in range(num):
    tmp = file.readline().split()
    for j in range(num):
        matrix[i][j] = int(tmp[j])
    matrix[i][i] = 0

(weight, pos) = solve(st, num, matrix)

for i in range(num):
    stdout.write(str(i) + " ")
stdout.write("P")
stdout.write("\n")

for i in range(num):
    stdout.write(str(weight[i]) + " ")
stdout.write("W")
stdout.write("\n")

for i in range(num):
    stdout.write(str(pos[i]) + " ")
stdout.write("R")
stdout.write("\n")

if(pos[fin] == -1):
    stdout.write(">NO WAY!\n")
    exit()

stdout.write(">ROOT CALCULATED!\n")

while(fin != st):
    stdout.write(str(fin) + " <- ")
    fin = pos[fin]
stdout.write(str(st) + "\n")
