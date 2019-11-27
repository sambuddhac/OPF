import sys, os
import numpy as np
import matplotlib.pyplot as plt


presFn= r'primresult.txt'
dresFn= r'dualresult.txt'
oresFn= r'objective.txt'

def readRes(fn):
	with open(fn) as f:
		r= [l.strip() for l in f]
	
	return [float(r[i]) for i in range(2, len(r))]
	
pe= readRes(presFn)
de= readRes(dresFn)
ob= readRes(oresFn)
fig1 = plt.figure()

pep, = plt.plot(pe, label= "Primal Error")
dep, =plt.plot(de, label= "Dual Error")
#oep, =plt.plot(ob, label= "Objective")
plt.legend(handles=[pep, dep])#, oep])
#plt.xlim(0, 1)
#plt.ylim(0, 1)
#plt.xlabel('x')
#plt.title('test')
plt.show()
