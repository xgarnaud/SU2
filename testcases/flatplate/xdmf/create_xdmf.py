#!/bin/python
import h5py
import CGNS
import CGNS.PAT.cgnsutils as CGU
import CHLone
from argparse import ArgumentParser
import numpy as np
import os


#              name       npoints ndims  cgnsid xdmf_type         xfmf_id
ELTYPES = {3 :['line'    ,2       ,1     ,3     ,'Polyline'      ,2],
           5 :['triangle',3       ,2     ,5     ,'Triangle'      ,4],
           9 :['quad'    ,4       ,2     ,7     ,'Quadrilateral' ,5],
           10:['tetra'   ,4       ,3     ,-1    ,'Tetrahedron'   ,6],
           12:['hexa'    ,8       ,3     ,-1    ,'Hexahedron'    ,9],
           13:['wedge'   ,6       ,3     ,-1    ,'Wedge'         ,8],
           14:['pyramid' ,5       ,3     ,-1    ,'Pyramid'       ,7],
           -1:['mixed'   ,-1      ,3    ,20     ,'Mixed'         ,-1],
       }
CGNSID_MIXED = 20

CG2SU = {}
for i in ELTYPES:
    CG2SU[ELTYPES[i][3]] = i

def get_cgmesh_connectivity(fname):

    (t,l,p)=CHLone.load(fname,lksearch=['.'])

    bases = CGU.hasChildType(t,'CGNSBase_t')
    if not len(bases) == 1:
        raise ValueError('there should be 1 base')
    base = bases[0]

    zones = CGU.hasChildType(base,'Zone_t')
    if not len(zones) == 1:
        raise ValueError('there should be 1 zone')
    zone = zones[0]

    # get the coordinates
    coords = CGU.hasChildType(zone,'GridCoordinates_t')
    if not len(coords) == 1:
        raise ValueError("invalid coordinates")
    x = CGU.getValueByPath(coords[0],'CoordinateX')
    y = CGU.getValueByPath(coords[0],'CoordinateY')
    z = CGU.getValueByPath(coords[0],'CoordinateZ')


    if abs(z).max() > 1e-10:
        ndims = 3
    else:
        ndims = 2

    # Get the connectivity
    elzones = {}
    conns  = CGU.hasChildType(zone,'Elements_t')
    for conn in conns:
        eltype_cg =  CGU.getValueByPath(conn,'.')[0]
        eltype_su = CG2SU[eltype_cg]

        elrange   = CGU.getValueByPath(conn,'ElementRange')
        elconn    = CGU.getValueByPath(conn,'ElementConnectivity')

        if not eltype_cg == CGNSID_MIXED:
            npts_per_el = ELTYPES[eltype_su][1]
            nel = len(elconn) / npts_per_el
            connectivity = elconn.reshape((nel,npts_per_el)) - 1
        else:
            nread = 0; ntot = len(elconn)
            while nread < ntot:
                tmp_elid_cg = elconn[nread];
                tmp_elid_su = CG2SU[tmp_elid_cg]
                tmp_elid_xf = ELTYPES[tmp_elid_su][5]
                elconn[nread] = tmp_elid_xf
                nread +=1
                tmp_npts    = ELTYPES[tmp_elid_su][1]
                elconn[nread:nread+tmp_npts] -= 1
                nread += tmp_npts
            connectivity = elconn

        elzones[conn[0]] = [eltype_su,elrange,connectivity]

    return elzones

def add_zones_h5(elzones,fname):

    h5f = h5py.File(fname,'a')

    for zone in elzones:
        if not zone in h5f:
            eltype,elrange,elconn = elzones[zone]
            dset = h5f.create_dataset(zone,data = elconn)
            dset.attrs['dim']  = ELTYPES[eltype][2]
            dset.attrs['type'] = ELTYPES[eltype][4]
            dset.attrs['nel'] = elrange[1] - elrange[0] +1

    h5f.close()

def create_xdmf(ofname,xfname):

    import xdmf
    h5f = h5py.File(ofname,'r')

    data   = h5f['PointData']

    if 'z' in data:
        ndims = 3
        coords_paths = ['PointData/x','PointData/y','PointData/z']
    else:
        ndims = 2
        coords_paths = ['PointData/x','PointData/y']

    fields = [x for x in data if not x in ['x','y','z']]
                        
    zones  = [x for x in h5f if not x in ['PointData']]

    nv = h5f['PointData/x'].shape[0]

    xf = xdmf.XDMFFile(xfname)
    for zone in zones:
        zdim  = h5f[zone].attrs['dim']
        zType = h5f[zone].attrs['type']
        s = h5f[zone].shape
        if len(s) == 2:
            nc,np_p_c = s
            mixed = False
        else:
            ntot = s
            mixed = True

        if h5f[zone].attrs['dim'] > 0:
            print zone
            xzone = xdmf.XDMFZone(zone)

            xzone.set_coords(ofname,coords_paths,nv)
            if mixed:
                nc = h5f[zone].attrs['nel']
                xzone.set_connectivity_mixed(zType,ofname,zone,nc,ntot)
            else:
                xzone.set_connectivity(zType,ofname,zone,nc,np_p_c)
            for field in fields:
                xzone.add_vfield(ofname,field,'PointData/%s'%field)
            xf.add_zone(xzone)

    xf.write()
    # print zones

    # print fields

    h5f.close()
    

if __name__=="__main__":
    
    parser=ArgumentParser()

    parser.add_argument("-m","--mesh", dest='mfile',
                        help='cgns mesh file',metavar='FILE')

    parser.add_argument("-o","--output", dest='ofile',
                        help='hdf5 restart file',metavar='FILE')

    parser.add_argument("-x","--xdmf", dest='xfile',
                        help='xdmf file',metavar='FILE')

    cwd = os.getcwd()

    options = parser.parse_args()

    elzones = get_cgmesh_connectivity(os.path.join(cwd,options.mfile))
    add_zones_h5(elzones,os.path.join(cwd,options.ofile))
    create_xdmf(os.path.join(cwd,options.ofile),os.path.join(cwd,options.xfile))
