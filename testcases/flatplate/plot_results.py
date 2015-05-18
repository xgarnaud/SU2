import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.close('all')

# reference data
tmp = np.loadtxt('cfl3d/sa-upyp_cfl3d.txt',skiprows=3)
yPlus_cfl3d = 10**tmp[:,0]#np.exp(tmp[:,0])
uPlus_cfl3d = tmp[:,1]
plt.figure('U - wall units')

plt.semilogx(yPlus_cfl3d,uPlus_cfl3d,'k--',label = 'cfl3d')

# Spalding formula
kappa  = .41
B      = 5.5
E      = np.exp(kappa*B)
plt.semilogx(uPlus_cfl3d + 1./E*(np.exp(kappa*uPlus_cfl3d)-1.-(kappa*uPlus_cfl3d)-(kappa*uPlus_cfl3d)**2/2.-(kappa*uPlus_cfl3d)**3/6.),uPlus_cfl3d,'r--',label = 'Spalding')

h5f = h5py.File('results.h5','r')
names = ["h_1e-5","h_1e-4","h_1e-3"]

for name in names:
    grp  = h5f[name]
    x1   = grp['x=1']
    y    = x1['y']
    U    = x1['U']
    nu   = x1['nu']
    nut  = x1['nut']
    uTau = x1.attrs['uTau']
    
    uPlus = U / uTau
    yPlus = y*uTau / nu
    nuPlus = nut / nu.value


    plt.figure('U')
    plt.plot(y,U,label = name)
    plt.xlabel('y')
    plt.ylabel('U')

    plt.figure('U - wall units')
    plt.semilogx(yPlus,uPlus,label = name)
    plt.xlabel('y+')
    plt.ylabel('U+')

    plt.figure('nut')
    plt.plot(y,nut)#,label = name)
    plt.xlabel('y')
    plt.ylabel('nut')

    plt.figure('nut - wall units')
    plt.semilogx(yPlus,nuPlus)#,label = name)
    plt.xlabel('y+')
    plt.ylabel('nut/nu')    

    wall = grp['wall']
    x = wall['x']
    Cf = wall['Cf']
    plt.figure('Cf')
    plt.plot(x,Cf)

    plt.ylim([0,.005])
    
plt.figure('U')
plt.legend(loc = 'best')

plt.figure('U - wall units')
plt.legend(loc = 'best')

for name in names:
    tmp = np.loadtxt('%s/surface_flow.csv'%name,delimiter = ',',skiprows=1)
    x = tmp[:,1]
    Cf = tmp[:,5]

    idx = x.argsort()
    x = x[idx]
    Cf = Cf[idx]

        
    
plt.show()
