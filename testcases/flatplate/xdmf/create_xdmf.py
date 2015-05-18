#!/bin/python
import h5py
import CGNS
import CGNS.PAT.cgnsutils as CGU
import CHLone
from argparse import ArgumentParser
import numpy as np


#              name       npoints ndims  cgnsid xdmf_type
ELTYPES = {3 :['line'    ,2       ,1     ,3     ,'Polyline'],
           5 :['triangle',3       ,2     ,5     ,'Triangle'],
           9 :['quad'    ,4       ,2     ,7     ,'Quadrilateral'],
           10:['tetra'   ,4       ,3     ,-1    ,'Tetrahedron'],
           12:['hexa'    ,8       ,3     ,-1    ,'Hexahedron'],
           13:['wedge'   ,6       ,3     ,-1    ,'Wedge'],
           14:['pyramid' ,5       ,3     ,-1    ,'Pyramid'],
       }

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
        eltype = -1
        for el in ELTYPES:
            if ELTYPES[el][3] == eltype_cg: eltype = el
        elrange   = CGU.getValueByPath(conn,'ElementRange')
        elconn    = CGU.getValueByPath(conn,'ElementConnectivity') - 1
        nel = len(elconn) / ELTYPES[eltype][1]
        elconn = elconn.reshape((nel,ELTYPES[eltype][1]))
        elzones[conn[0]] = [eltype,elrange,elconn]

    return elzones

def add_zones_h5(elzones,fname):

    h5f = h5py.File(fname,'a')

    for zone in elzones:
        if not zone in h5f:
            eltype,elrange,elconn = elzones[zone]
            dset = h5f.create_dataset(zone,data = elconn)
            dset.attrs['dim']  = ELTYPES[eltype][2]
            dset.attrs['type'] = ELTYPES[eltype][4]
        
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
        nc,np_p_c =  h5f[zone].shape
        if h5f[zone].attrs['dim'] > 0:
            xzone = xdmf.XDMFZone(zone)

            xzone.set_coords(ofname,coords_paths,nv)
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



    options = parser.parse_args()

    elzones = get_cgmesh_connectivity(options.mfile)
    add_zones_h5(elzones,options.ofile)
    create_xdmf(options.ofile,options.xfile)
