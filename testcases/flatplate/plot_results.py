import numpy as np
import matplotlib.pyplot as plt


plt.close('all')

# reference data
tmp = np.loadtxt('cfl3d/sa-upyp_cfl3d.dat',skiprows=3)
yPlus_cfl3d = 10**tmp[:,0]#np.exp(tmp[:,0])
uPlus_cfl3d = tmp[:,1]
plt.figure('U - wall units')

plt.semilogx(yPlus_cfl3d,uPlus_cfl3d,'k--',label = 'cfl3d')

# Spalding formula
kappa  = .41
B      = 5.5
E      = np.exp(kappa*B)
plt.semilogx(uPlus_cfl3d + 1./E*(np.exp(kappa*uPlus_cfl3d)-1.-(kappa*uPlus_cfl3d)-(kappa*uPlus_cfl3d)**2/2.-(kappa*uPlus_cfl3d)**3/6.),uPlus_cfl3d,'r--',label = 'Spalding')

names = ["h_1e-6","h_1e-5","h_1e-4","h_1e-3"]

for name in names:
    tmp = np.loadtxt('%s/slice_x=1.dat'%name)
    plt.figure('U')
    plt.plot(tmp[:,0],tmp[:,1],label = name)
    plt.xlabel('y')
    plt.ylabel('U')

    plt.figure('U - wall units')
    plt.semilogx(tmp[:,3],tmp[:,4],label = name)
    plt.xlabel('y+')
    plt.ylabel('U+')

    plt.figure('nut')
    plt.plot(tmp[:,0],tmp[:,2])#,label = name)
    plt.xlabel('y')
    plt.ylabel('nut')

    plt.figure('nut - wall units')
    plt.semilogx(tmp[:,3],tmp[:,5])#,label = name)
    plt.xlabel('y+')
    plt.ylabel('nut/nu')

    # plt.figure('sa_prod')
    # plt.plot(tmp[:,0],tmp[:,6])

    # plt.figure('sa_dest')
    # plt.plot(tmp[:,0],tmp[:,7])

    # plt.figure('sa_crossprod')
    # plt.plot(tmp[:,0],tmp[:,8])
    
    
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

    plt.figure('Cf')
    plt.plot(x,Cf)

    plt.ylim([0,.005])
        
    
plt.show()
