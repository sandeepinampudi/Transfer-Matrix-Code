class TMM:
  
    def __init__(self,epsxx,epsyz,xd,kz,k0,pol):
        import cmath as math 
        self.epsxx=epsxx
        self.epsyz=epsyz
        self.xd=xd
        self.kz=kz
        self.k0=k0
        self.pol=pol
        self.kx=[]
        self.R=[]
        self.T=[]        
        
        for epsx,epsy in zip(self.epsxx,self.epsyz):                     
            if self.pol=='TE':
               self.kx.append(math.sqrt(k0**2*epsy-kz**2))
            else:
               self.kx.append(math.sqrt(k0**2*epsy-kz**2*epsy/epsx))
               
        n=len(self.epsxx)
        for epsx,epsy in zip(reversed(self.epsxx),reversed(self.epsyz)):                     
            if n==len(self.epsxx):
                R=0
                T=1                
            else:            
                kx1=self.kx[n-1]
                kx2=self.kx[n]
                phi1=math.exp((0+1j)*kx1*self.xd[n-1])
                phi2=math.exp((0+1j)*kx2*self.xd[n])
        
                if self.pol=='TE':
                    D=(kx1+kx2)
                    Rpl=(kx1-kx2)/D
                    Tpl=(2*kx1)/D
                    Tmi=(2*kx2)/D
                    Rmi=(kx2-kx1)/D                    
                else:                                      
                    k1e2=kx1*self.epsyz[n]
                    k2e1=kx2*self.epsyz[n-1]
                    
                    D=(k2e1+k1e2)                    
                    Rpl=-(k2e1-k1e2)/D
                    Tpl=(2*k1e2)/D
                    Tmi=(2*k2e1)/D
                    Rmi=(k2e1-k1e2)/D
                
                T=(Tpl*phi1)/(1-Rmi*phi2*R)
                R=Rpl*phi1+Tmi*phi2*R*T 
            self.R.append(R)
            self.T.append(T)            
            n=n-1
    def Ref(self):
        return abs(self.R[len(self.R)-1])**2 
    def Tran(self):
        Tend=1
        for T in self.T:
            Tend=Tend*T 
        if self.pol=='TE':
           return abs(Tend)**2*self.kx[-1].real/self.kx[0].real
        else:
           a1=self.kx[0]/self.epsyz[0]
           an=self.kx[-1]/self.epsyz[-1]
           return abs(Tend)**2*an.real/a1.real
    
