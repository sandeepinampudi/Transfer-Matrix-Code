# Transfer matrix code for light propagation through planar uniaxial layers
# Current example computes the transmission coeffcicnt of light through
# distributed bragg mirrors with a defect layer
 
from TransferMatrix import  TMM
import cmath as math 
import numpy as np

n1=1;
n2=6;

lam00=1.55;

d1=lam00/(4*n1)
d2=lam00/(4*n2)


eps1=n1**2;
eps2=n2**2;

# provide a list of permittivity of each layer here
epsxx=[1,eps2,eps1,eps2,eps1,eps2,eps2,eps1,eps2,eps1,eps2,1] 
epsyz=[1,eps2,eps1,eps2,eps1,eps2,eps2,eps1,eps2,eps1,eps2,1]
# provide a list of thickness of each layer here
xd=[0,d2,d1,d2,d1,d2,d2,d1,d2,d1,d2,0]
pol='TE' # Polarization of light
kz=0 # Incident wavevector =k0*sin(theta) for incident angle theta  
lamLst=np.linspace(0.5,4.5,16001) # List of wavelenghts to consider
R=[]
T=[]
for lam0 in lamLst:
    k0=2*math.pi/lam0
    myf = TMM(epsxx,epsyz,xd,kz,k0,pol)
    R.append(myf.Ref())
    T.append(myf.Tran())

import matplotlib.pyplot as plt
plt.plot(lamLst,T)
plt.xlabel('Wavelength (microns)')
plt.ylabel('Transmission')


