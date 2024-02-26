

def mesh_nodes(NV,nelr,nelt,axisymmetric,R1,R2):


    xV=np.zeros(NV,dtype=np.float64) 
    yV=np.zeros(NV,dtype=np.float64) 
    rad=np.zeros(NV,dtype=np.float64)  
    theta=np.zeros(NV,dtype=np.float64) 

    if axisymmetric:



    else:



    return xV,yV,rad,theta
