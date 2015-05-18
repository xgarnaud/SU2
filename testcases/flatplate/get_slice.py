import sys,os
sys.path.append('./xdmf')

import numpy as np
from scipy.interpolate import interp1d
import create_xdmf as cx
from vtk import vtkXdmfReader,vtkUnstructuredGridReader,vtkLineSource,vtkPoints,vtkProbeFilter
from vtk.util import numpy_support as VN

def readXdmf(filename):
    #read the xdmf file
    print 'Reading ',filename,
    reader = vtkXdmfReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    print 'done'

    print 'Point arrays:'
    narrays = reader.GetNumberOfPointArrays()
    for i in range(narrays):
        print '  ',reader.GetPointArrayName(i).strip()

    iout = -1
    print 'Grids:'
    ngrids = reader.GetNumberOfGrids()
    for i in range(ngrids):
        gname = reader.GetGridName(i)
        print '  ',gname,', status = ',reader.GetGridStatus(gname)

    return reader

def getXdmfGrid(reader,blockname):
    
    iout = -1
    ngrids = reader.GetNumberOfGrids()
    for i in range(ngrids):
        gname = reader.GetGridName(i)
        if gname == blockname:
           iout = i
           break

    if iout < 0:
        raise ValueError('Block %s not found'%blockname)

    return reader.GetOutputDataObject(0).GetBlock(iout)

def readVTK(filename):
    #read the vtk file with an unstructured grid
    print 'Reading ',filename,
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    print 'done'

    return reader

def readVTKGrid(reader):

    return reader.GetOutput()

def createLine(p1,p2,h1,stretch):
    
    print 'Create a line:'
    print '  p1    = ',p1
    print '  p2    = ',p2
    print '  h(p1) = ',h1

    L = np.linalg.norm(p2 - p1)
    l = [0]; dl = h1 / L
    while l[-1] < 1.:
        l.append(l[-1] + dl)
        dl *= stretch
    l[-1] = 1.
    n = len(l)
    
    print '  hmax = ',dl
    print '  n    = ',n
    
    sourcePoints = vtkPoints()
    for i in range(n):
        p = p1 + l[i]*(p2 - p1)
        sourcePoints.InsertNextPoint(p[0],p[1],p[2])
        
    # Create the line along which you want to sample
    line = vtkLineSource()
    line.SetPoints(sourcePoints)
    line.Update()
    return line

def probeOverLine(line,data,varnames):

    probe = vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(data)
    probe.Update()

    return getArrays(probe.GetOutput(),varnames)

def getArrays(data,varnames,sortx = False,extract = False):
    
    npts   = data.GetNumberOfPoints()
    ncells = data.GetNumberOfCells()
    
    x = np.zeros(npts)
    y = np.zeros(npts)
    z = np.zeros(npts)

    for i in range(npts):
        x[i],y[i],z[i] = data.GetPoint(i)

    pdata = data.GetPointData()

    res = [VN.vtk_to_numpy(pdata.GetArray(varname)) for varname in varnames]

    if extract:
        idx = []
        for icell in range(ncells):
            cell = data.GetCell(icell)
            ids  = cell.GetPointIds()
            npcell = ids.GetNumberOfIds()
            for ip in range(npcell):
                ii = ids.GetId(ip)
                if not ii in idx:
                    idx.append(ii)       
        x   = x[idx]
        y   = y[idx]
        z   = z[idx]
        res = [tmp[idx] for tmp in res]

    if sortx:
        idx = x.argsort()
        idx = idx[:-1]
        x   = x[idx]
        y   = y[idx]
        z   = z[idx]
        res = [tmp[idx] for tmp in res]
        
        
    return x,y,z,res


if __name__=="__main__":

    import matplotlib.pyplot as plt
    import h5py as h5

    h5f = h5.File('results.h5','a')

    case      = "h_1e-3"
    mesh      = "flatplate_1e-3.cgns"
    h         = 1e-3

    dir_name  = os.path.join(os.getcwd(),case)
    mesh_name = os.path.join(dir_name,mesh)
    res_name  = os.path.join(dir_name,"restart_flow.dat")
    xdmf_name = os.path.join(dir_name,"restart.xdmf")

    elzones = cx.get_cgmesh_connectivity(mesh_name)
    cx.add_zones_h5(elzones,res_name)
    cx.create_xdmf(res_name, xdmf_name )

    reader  = readXdmf(xdmf_name)
    fluid   = getXdmfGrid(reader,'FLUID')
    wall    = getXdmfGrid(reader,'WALL')

    p1      = np.array([1,0 ,0])
    p2      = np.array([1,.5,0])
    line    = createLine(p1,p2,h,1.025)

    # wall data
    x,y,z,Q = getArrays(wall,['Y_Plus','Skin_Friction_Coefficient'],sortx = True,extract = True)
    yPlus = Q[0]
    Cf    = Q[1]

    if case in h5f:
        del h5f[case]

    grp      = h5f.create_group(case)
    h5f_wall = grp.create_group('wall')
    dset = h5f_wall.create_dataset('x'     ,data = x     )
    dset = h5f_wall.create_dataset('Cf'    ,data = Cf    )
    dset = h5f_wall.create_dataset('yPlus' ,data = yPlus )
        
    plt.close('all')

    plt.subplot(221)
    plt.plot(x,Cf)

    plt.subplot(222)
    plt.plot(x,yPlus)

    Cf_i    = interp1d(x,Cf)
    yPlus_i = interp1d(x,yPlus)    
    Cf      = Cf_i(1.0)

    # Fluid data
    x,y,z,Q = probeOverLine(line,fluid,['Conservative_1','Conservative_2','Laminar_Viscosity','Eddy_Viscosity'])
    rho  = Q[0]
    rhoU = Q[1]; U = rhoU / rho
    mu   = Q[2]; nu = mu / rho
    nut  = Q[3]
    
    rhoRef   = 1.32905
    Uref     = 69.4448
    tauw     = .5*rhoRef*Uref**2 * Cf
    uTau     = (tauw/rhoRef)**.5

    wallShearStress = rhoRef * uTau**2

    uPlus = U / uTau
    yPlus = y*uTau / nu

    h5f_x1 = grp.create_group('x=1')
    h5f_x1.create_dataset('y'  ,data = y)
    h5f_x1.create_dataset('rho',data = rho)
    h5f_x1.create_dataset('U'  ,data = U)
    h5f_x1.create_dataset('nu' ,data = nu)
    h5f_x1.create_dataset('nut',data = nut)
    h5f_x1.attrs['uTau'] = uTau

    h5f.close()

    plt.subplot(223)
    plt.semilogx(yPlus,uPlus)

    plt.subplot(224)
    plt.semilogx(yPlus,nut)
    
    plt.show()
