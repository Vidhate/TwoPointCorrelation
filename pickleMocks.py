import pickle
import numpy as np

#path="./../../DLA Mock Catalogue/"
path="./Data/"
names=["DLAhost_snap49_r1_b1 (copy).txt","DLAhost_snap49_r1_b1alpha (copy).txt","DLAhost_snap49_r1_b1T10 (copy).txt"]

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
