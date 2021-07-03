import os,sys

f = open("typing_tree.log")
lines = f.readlines()
f.close()

for line in lines:
  terms = line.split()
  num_1 = int(terms[0])
  num_2 = int(terms[1])
  num_3 = int(terms[2])
  num_4 = int(terms[3])
  if(num_2 >= 14):
    num_2 += 1
  if(num_3 >= 14):
    num_3 += 1
  if(num_4 >= 14):
    num_4 += 1
  print("%5d %5d %5d %5d" %(num_1, num_2, num_3, num_4))
