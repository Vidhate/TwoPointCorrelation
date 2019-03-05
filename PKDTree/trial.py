import scipy.spatial as sp
import numpy as np

# Boundaries (0 or negative means open boundaries in that dimension)
bounds = np.array([30.0, 30.0, 30.0])   # xy periodic, open along z

# Points
n = 1000
data = 30.0 * np.random.randn(n, 3)

# Build kd-tree
#T = PeriodicCKDTree(bounds, data)
T_self=sp.cKDTree(data)
T_other=sp.cKDTree(data)
neighbors=T_self.count_neighbors(T_other,10)
'''
#build other (periodic imaged) dataset
T_other=#build tree

corr=[]
for i in range(0,len(bins)):
	corr.append(T_self.count_neighbors(T_other,bins[i]))

corr=np.array(corr[1:])-np.array(corr[:-1])
neighbors=
'''
print(neighbors)
