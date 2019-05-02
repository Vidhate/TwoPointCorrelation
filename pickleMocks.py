import pickle
import numpy as np

#path="./../../DLA Mock Catalogue/"
path="./Data/"
names=["b1.txt","b1alpha.txt","b1T10.txt"]

for name in names:
	pos=[]
	fname=path+name
	with open(fname) as f:
		for line in f:
			if line[0]!='#':
				col=line.split(' ')
				cood=np.array([float(col[1]),float(col[2]),float(col[3])])
				pos.append(cood)
	pos=np.array(pos)
	pickle.dump(pos,open(fname+"_pickled","wb"))
