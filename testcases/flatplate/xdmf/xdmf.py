import h5py
import numpy as np
import os

class XDMFZone(object):

    def __init__(self,name):

        self.name      = name
        self.nc      = 0
        self.nv      = 0
        self.vfields = []
        self.cfields = []

    def set_coords(self,coords_file,coords_path,nv):

        self.coords_file = coords_file
        self.coords_path = coords_path
        self.ndims = len(coords_path)
        self.nv = nv

    def set_connectivity(self,conn_type,conn_file,conn_path,nc,np_p_c):
        
        self.conn_type = conn_type
        self.conn_file = conn_file
        self.conn_path = conn_path
        self.nc = nc
        self.np_p_c = np_p_c
        self.mixed = False

    def set_connectivity_mixed(self,conn_type,conn_file,conn_path,nc,ntot):
        
        self.conn_type = conn_type
        self.conn_file = conn_file
        self.conn_path = conn_path
        self.nc = nc
        self.ntot = ntot
        self.mixed = True

    def add_vfield(self,field_file,field_name,field_path):

        self.vfields.append([field_file,field_name,field_path])

    def add_cfield(self,field_file,field_name,field_path):

        self.cfields.append([field_file,field_name,field_path])

    def write_simple(self,f):

        ie = 0
        f.write('         <Grid Name="%s">\n'%self.name)
        f.write('           <Topology TopologyType="%s" NumberOfElements="%d" >\n'%(self.conn_type,self.nc))
        if self.mixed:
            f.write('             <DataItem Format="HDF" DataType="Int" Dimensions="%d">\n'%(self.ntot))
        else:
            f.write('             <DataItem Format="HDF" DataType="Int" Dimensions="%d %d">\n'%(self.nc,self.np_p_c))
        f.write('               %s:%s\n'%(self.conn_file,self.conn_path))
        f.write('             </DataItem>\n')
        f.write('           </Topology>\n')
        f.write('           \n')
        if self.ndims == 2:
            f.write('           <Geometry GeometryType="X_Y">\n')
        else:
            f.write('           <Geometry GeometryType="X_Y_Z">\n')
        f.write('             <DataItem Format="HDF" Dimensions="%d">\n'%self.nv)
        f.write('               %s:%s\n'%(self.coords_file,self.coords_path[0]))
        f.write('             </DataItem>\n')
        f.write('             <DataItem Format="HDF" Dimensions="%d">\n'%self.nv)
        f.write('               %s:%s\n'%(self.coords_file,self.coords_path[1]))
        f.write('             </DataItem>\n')
        if self.ndims == 3:
            f.write('             <DataItem Format="HDF" Dimensions="%d">\n'%self.nv)
            f.write('               %s:%s\n'%(self.coords_file,self.coords_path[2]))
            f.write('             </DataItem>\n')
        f.write('           </Geometry>\n')
        f.write('           \n')
        for field in self.cfields:
            field_file,field_name,field_path = field
            f.write('           <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n'%field_name)
            f.write('             <DataItem Dimensions="%d 1" Format="HDF">\n'%self.nc)
            f.write('               %s:%s\n'%(field_file,field_path))
            f.write('             </DataItem>\n')
            f.write('           </Attribute>\n')
        for field in self.vfields:
            field_file,field_name,field_path = field
            f.write('           <Attribute Name="%s" AttributeType="Scalar" Center="Node">\n'%field_name)
            f.write('             <DataItem Dimensions="%d 1" Format="HDF">\n'%self.nv)
            f.write('               %s:%s\n'%(field_file,field_path))
            f.write('             </DataItem>\n')
            f.write('           </Attribute>\n')

        f.write('         </Grid>\n')


class XDMFFile(object):

    def __init__(self,fname):

        self.fname = fname
        self.zones = []

    def add_zone(self,zone):

        self.zones.append(zone)

    def write(self):

        f = open(self.fname,'w')
        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf Version="2.0">\n')
        f.write('  <Domain>\n')
        f.write('    <Grid CollectionType="Spatial" GridType="Collection">\n')

        for z in self.zones:
            z.write_simple(f)

        f.write('    </Grid>\n')
        f.write('  </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()
