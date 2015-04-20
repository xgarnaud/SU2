import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from vtk import vtkUnstructuredGridReader,vtkProbeFilter,vtkLineSource
from vtk.util import numpy_support as VN

def readVTK(filename):
    #read the vtk file with an unstructured grid
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    return reader

def createLine(p1,p2,numPoints):
    # Create the line along which you want to sample
    line = vtkLineSource()
    line.SetResolution(numPoints)
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    line.Update()
    return line

def probeOverLine(line,reader,varname):
    #Interpolate the data from the VTK-file on the created line.
    data = reader.GetOutput()
    # vtkProbeFilter, the probe line is the input, and the underlying dataset is the source.
    probe = vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(data)
    probe.Update()
    #get the data from the VTK-object (probe) to an numpy array
    try:
        q=VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray(varname))
    except AttributeError:
        q=VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray('"'+varname+'"'))
    numPoints = probe.GetOutput().GetNumberOfPoints() # get the number of points on the line
    #intialise the points on the line    
    x = np.zeros(numPoints)
    y = np.zeros(numPoints)
    z = np.zeros(numPoints)
    points = np.zeros((numPoints , 3))
    #get the coordinates of the points on the line
    for i in range(numPoints):
        x[i],y[i],z[i] = probe.GetOutput().GetPoint(i)
        points[i,0]=x[i]
        points[i,1]=y[i]
        points[i,2]=z[i]
    return points,q

def readvars(reader,variables):

    data = reader.GetOutput()

    npts   = data.GetNumberOfPoints()
    ncells = data.GetNumberOfCells()

    x = np.zeros(npts)
    y = np.zeros(npts)
    z = np.zeros(npts)

    for i in range(npts):
        x[i],y[i],z[i] = data.GetPoint(i)

    idx = x.argsort()
    idx = idx[:-1]
    x = x[idx]

    pdata = data.GetPointData()

    na = pdata.GetNumberOfArrays()

    res = []
    for vname in variables:
        arr = pdata.GetArray(vname)
        if arr == None:
            arr = pdata.GetArray('"'+vname+'"')
        nparr = VN.vtk_to_numpy(arr)
        res.append(nparr[idx])

    return x,res


extra_out = False

d    = 'h_1e-3_tri'
h0   = 2e-3

# d    = 'h_1e-6'
# h0   = 1e-6
# name = 'h1e-6_nowm'

p1=[1.,0.0 ,0.0]
p2=[1.,0.05,0.0]

reader2d = readVTK('%s/flow.vtk'%d)
reader1d = readVTK('%s/surface_flow.vtk'%d)

x,[yPlus,Cf] = readvars(reader1d,['Y_Plus','Skin_Friction_Coefficient'])
Cf_i    = interp1d(x,Cf)
yPlus_i = interp1d(x,yPlus)

Cf = Cf_i(1.0)

# print x,yPlus,Cf
# line=createLine(p1,p2,numPoints = 2)
# points,Cf  =  probeOverLine(line,reader1d,'Skin_Friction_Coefficient')
# points,yPlus  =  probeOverLine(line,reader1d,'Y_Plus')
# points,mu  =  probeOverLine(line,reader1d,'Laminar_Viscosity')

# print points,Cf
# Cf     = Cf[0]
# mu     = mu[0]
# print Cf,mu,yPlus

p1=[1.,0.  ,0.0]
p2=[1.,0.05,0.0]

line=createLine(p1,p2,numPoints = int(0.05/h0))
points,rho  =  probeOverLine(line,reader2d,'Conservative_1')
points,rhoU =  probeOverLine(line,reader2d,'Conservative_2')
points,mut  =  probeOverLine(line,reader2d,'Eddy_Viscosity')
points,mu   =  probeOverLine(line,reader2d,'Laminar_Viscosity')

if extra_out:
    points,prod       =  probeOverLine(line,reader2d,'Production')
    points,dest       =  probeOverLine(line,reader2d,'Destruction')
    points,crossprod  =  probeOverLine(line,reader2d,'CrossProduction')

y1    = points[1,1]
U1    = rhoU[1] / rho[1]
rho0  = rho[0]

mu0   = mu[0]
mut0  = mut[0]
nu    = mu0 /rho0

wallShearStress = (mu0 + mut0) * U1 / y1

rhoRef   = 1.32905
Uref     = 69.4448
tauw     = .5*rhoRef*Uref**2 * Cf
uTau     = (tauw/rho0)**.5

wallShearStress = rho0 * uTau**2

uPlus = rhoU/rho / uTau
yPlus = points[:,1]*uTau / nu

npts = len(points[:,1])
f = open('%s/slice_x=1.dat'%d,'w')
for i in range(npts):
    if extra_out:
        f.write('%g %g %g %g %g %g %g %g %g\n'%(points[i,1],rhoU[i]/rho[i],mut[i],yPlus[i],uPlus[i],mut[i]/mu[i],prod[i],dest[i],crossprod[i]))
    else:
        f.write('%g %g %g %g %g %g\n'%(points[i,1],rhoU[i]/rho[i],mut[i],yPlus[i],uPlus[i],mut[i]/mu[i]))
f.close()
