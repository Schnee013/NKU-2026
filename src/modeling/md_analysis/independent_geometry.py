"""Independent NumPy SVD reference fit."""
import numpy as np

def independent_fit(mobile, reference, xyz):
    mc=mobile.mean(0);rc=reference.mean(0)
    U,_,V=np.linalg.svd((mobile-mc).T @ (reference-rc))
    D=np.diag([1,1,np.linalg.det(U@V)])
    return (xyz-mc) @ U @ D @ V + rc

