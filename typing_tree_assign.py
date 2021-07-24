import os,sys
from typing_tree import typing_tree
import numpy as np

def typing_tree_assign(typing_tree, term, comb, classParameterDict):

  # Bond stretching terms
  if(term == 'b'):
    atoms = comb.split('_')
    lists_1 = typing_tree.search_in_tree(max(int(atoms[0]),int(atoms[1])))
    lists_2 = typing_tree.search_in_tree(min(int(atoms[0]),int(atoms[1])))
    
    # check if there is H (HCH, HCH2, HCH3 and so on)
    if(int(atoms[0]) == 1 or int(atoms[0]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[1]
        comb_2 = atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        if(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])
    if(int(atoms[1]) == 1 or int(atoms[1]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])

    # check if there is any canonical sp3 C
    collect = []
    collect_1 = []
    if(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j)
          comb_2 = str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1]
        comb_2 = atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[1]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:,0].mean()
      ave1 = np.array(collect_1)[:,1].mean()
      average = (ave0, ave1)
      return (comb, average)
    else:
      pass

    # canonical searching
    collect = []
    collect_1 = []
    for i_1 in lists_1:
      for i_2 in lists_2:
        comb_1 = str(i_1) + '_' + str(i_2)
        comb_2 = str(i_2) + '_' + str(i_1)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        if(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:5,0].mean()
      ave1 = np.array(collect_1)[:5,1].mean()
      average = (ave0, ave1)
      return (comb, average)
    else:
      pass

  # Angle bending terms
  if(term == 'a'):
    atoms = comb.split('_')
    lists_1 = typing_tree.search_in_tree(int(atoms[0]))
    lists_2 = typing_tree.search_in_tree(int(atoms[1]))
    lists_3 = typing_tree.search_in_tree(int(atoms[2]))

    # check if there is H (HCH, HCH2, HCH3 and so on)
    if(int(atoms[0]) == 1 or int(atoms[0]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[2]
        comb_2 = atoms[2] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        if(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])
    if(int(atoms[2]) == 1 or int(atoms[2]) in range(16,29)):
      for i in lists_3:
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])

    # check if there is any canonical sp3 C
    collect = []
    collect_1 = []
    if(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          for k in range(32,36):
            comb_1 = str(i) + '_' + str(j) + '_' + str(k)
            comb_2 = str(k) + '_' + str(j) + '_' + str(i)
            if(comb_1 in classParameterDict):
              collect.append((comb_1,classParameterDict[comb_1]))
              collect_1.append(classParameterDict[comb_1])
            elif(comb_2 in classParameterDict):
              collect.append((comb_2,classParameterDict[comb_2]))
              collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j) + '_' + atoms[2]
          comb_2 = atoms[2] + '_' + str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j) + '_' + atoms[0]
          comb_2 = atoms[2] + '_' + str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + atoms[1] + '_' + str(j)
          comb_2 = str(j) + '_' + atoms[1] + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[2]
        comb_2 = atoms[2] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        comb_1 = atoms[2] + '_' + str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i) + '_' + atoms[2]
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:,0].mean()
      ave1 = np.array(collect_1)[:,1].mean()
      average = (ave0, ave1)
      return (comb, average)
    else:
      pass

    # canonical searching
    collect = []
    collect_1 = []
    for i_2 in lists_2:
      for i_1 in lists_1:
        for i_3 in lists_3:
          comb_1 = str(i_1) + '_' + str(i_2) + '_' + str(i_3)
          comb_2 = str(i_3) + '_' + str(i_2) + '_' + str(i_1)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          if(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:5, 0].mean()
      ave1 = np.array(collect_1)[:5, 1].mean()
      average = (ave0, ave1)
      return (comb, average)
    else:
      pass

  # strbnd terms
  if(term == 'ba'):
    atoms = comb.split('_')
    lists_1 = typing_tree.search_in_tree(int(atoms[0]))
    lists_2 = typing_tree.search_in_tree(int(atoms[1]))
    lists_3 = typing_tree.search_in_tree(int(atoms[2]))

    # check if there is H (HCH, HCH2, HCH3 and so on)
    if(int(atoms[0]) == 1 or int(atoms[0]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[2]
        comb_2 = atoms[2] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        if(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])
    if(int(atoms[2]) == 1 or int(atoms[2]) in range(16,29)):
      for i in lists_3:
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, classParameterDict[comb_1])
        if(comb_2 in classParameterDict):
          return (comb, classParameterDict[comb_2])

    # check if there is any canonical sp3 C
    collect = []
    collect_1 = []
    if(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          for k in range(32,36):
            comb_1 = str(i) + '_' + str(j) + '_' + str(k)
            comb_2 = str(k) + '_' + str(j) + '_' + str(i)
            if(comb_1 in classParameterDict):
              collect.append((comb_1,classParameterDict[comb_1]))
              collect_1.append(classParameterDict[comb_1])
            if(comb_2 in classParameterDict):
              collect.append((comb_2,classParameterDict[comb_2]))
              collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j) + '_' + atoms[2]
          comb_2 = atoms[2] + '_' + str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j) + '_' + atoms[0]
          comb_2 = atoms[2] + '_' + str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + atoms[1] + '_' + str(j)
          comb_2 = str(j) + '_' + atoms[1] + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[2]
        comb_2 = atoms[2] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) in range(32,36) and int(atoms[1]) not in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1] + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    elif(int(atoms[0]) not in range(32,36) and int(atoms[1]) in range(32,36) and int(atoms[2]) not in range(32,36)):
      for i in range(32,36):
        comb_1 = atoms[2] + '_' + str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i) + '_' + atoms[2]
        if(comb_1 in classParameterDict):
          collect.append((comb_1,classParameterDict[comb_1]))
          collect_1.append(classParameterDict[comb_1])
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,classParameterDict[comb_2]))
          collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:,0].mean()
      ave1 = np.array(collect_1)[:,1].mean()
      average = (ave0, ave1) 
      return (comb, average)
    else:
      pass

    # canonical searching
    collect = []
    collect_1 = []
    for i_2 in lists_2:
      for i_1 in lists_1:
        for i_3 in lists_3:
          comb_1 = str(i_1) + '_' + str(i_2) + '_' + str(i_3)
          comb_2 = str(i_3) + '_' + str(i_2) + '_' + str(i_1)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,classParameterDict[comb_1]))
            collect_1.append(classParameterDict[comb_1])
          if(comb_2 in classParameterDict):
            collect.append((comb_2,classParameterDict[comb_2]))
            collect_1.append(classParameterDict[comb_2])
    if(len(collect) != 0):
      ave0 = np.array(collect_1)[:5,0].mean()
      ave1 = np.array(collect_1)[:5,1].mean()
      average = (ave0, ave1)
      return (comb, average)
    else:
      pass

  # Opbend terms
  if(term == 'o'):
    atoms = comb.split('_')
    lists_1 = typing_tree.search_in_tree(max(int(atoms[0]),int(atoms[1])))
    lists_2 = typing_tree.search_in_tree(min(int(atoms[0]),int(atoms[1])))

    # check if there is H (HCH, HCH2, HCH3 and so on)
    if(int(atoms[0]) == 1 or int(atoms[0]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[1]
        comb_2 = atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, float(classParameterDict[comb_1]))
        elif(comb_2 in classParameterDict):
          return (comb, float(classParameterDict[comb_2]))
    if(int(atoms[1]) == 1 or int(atoms[1]) in range(16,29)):
      for i in lists_1:
        comb_1 = str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i)
        if(comb_1 in classParameterDict):
          return (comb, float(classParameterDict[comb_1]))
        elif(comb_2 in classParameterDict):
          return (comb, float(classParameterDict[comb_2]))

    # check if there is any canonical sp3 C
    collect = []
    collect_1 = []
    if(int(atoms[0]) in range(32,36) and int(atoms[1]) in range(32,36)):
      for i in range(32,36):
        for j in range(32,36):
          comb_1 = str(i) + '_' + str(j)
          comb_2 = str(j) + '_' + str(i)
          if(comb_1 in classParameterDict):
            collect.append((comb_1,float(classParameterDict[comb_1])))
            collect_1.append(float(classParameterDict[comb_1]))
          elif(comb_2 in classParameterDict):
            collect.append((comb_2,float(classParameterDict[comb_2])))
            collect_1.append(float(classParameterDict[comb_2]))
    elif(int(atoms[0]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[1]
        comb_2 = atoms[1] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,float(classParameterDict[comb_1])))
          collect_1.append(float(classParameterDict[comb_1]))
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,float(classParameterDict[comb_2])))
          collect_1.append(float(classParameterDict[comb_2]))
    elif(int(atoms[1]) in range(32,36)):
      for i in range(32,36):
        comb_1 = str(i) + '_' + atoms[0]
        comb_2 = atoms[0] + '_' + str(i)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,float(classParameterDict[comb_1])))
          collect_1.append(float(classParameterDict[comb_1]))
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,float(classParameterDict[comb_2])))
          collect_1.append(float(classParameterDict[comb_2]))
    if(len(collect) != 0):
      average = np.array(collect_1).mean()
      return (comb, average)
    else:
      pass

    # canonical searching
    collect = []
    collect_1 = []
    for i_1 in lists_1:
      for i_2 in lists_2:
        comb_1 = str(i_1) + '_' + str(i_2)
        comb_2 = str(i_2) + '_' + str(i_1)
        if(comb_1 in classParameterDict):
          collect.append((comb_1,float(classParameterDict[comb_1])))
          collect_1.append(float(classParameterDict[comb_1]))
        elif(comb_2 in classParameterDict):
          collect.append((comb_2,float(classParameterDict[comb_2])))
          collect_1.append(float(classParameterDict[comb_2]))
    if(len(collect) != 0):
      average = np.array(collect_1[:5]).mean()
      return (comb, average)
    else:
      pass


#a = typing_tree()
#a.read_ranking_file('typing_tree.log','amoebaplusBondedType.dat')
#a.sorting_tree()
#class1, class2, kbonds = np.loadtxt("opbend.prm", usecols=(1,2,5), unpack=True, dtype="str",skiprows=1)
#classes = []
#for c1, c2 in zip(class1, class2):
#  classes.append(c1 + "_" + c2)
#classParameterDict = dict(zip(classes, kbonds))
#
#comb = []
#f = open("undefined_2")
#lines = f.readlines()
#f.close()
#for line in lines:
#  terms = line.split()
#  comb.append(terms[0] + '_' + terms[1])
#for k in comb:
#  i = typing_tree_assign(a, 'o', k, classParameterDict)
#  atoms = i[0].split('_')
#  print("opbend %5d %6d %6d %6d %15.4f" %(int(atoms[0]), int(atoms[1]), 0, 0, float(i[1])))
