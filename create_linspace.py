from sys import argv
from numpy import linspace,savetxt
bi=int(argv[2])
bs=int(argv[3])
num=int(argv[4])
arr=linspace(bi,bs,num)
file = open(argv[1], "w+")
savetxt(file,arr)


