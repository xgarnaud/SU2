import sys
import CGNS
import CGNS.PAT.cgnsutils as CGU
import CHLone
from argparse import ArgumentParser
import numpy as np

#              name       npoints ndims  cgnsid
ELTYPES = {3 :['line'    ,2       ,1     ,3 ],
           5 :['triangle',3       ,2     ,5 ],
           9 :['quad'    ,4       ,2     ,7 ],
           10:['tetra'   ,4       ,3     ,-1],
           12:['hexa'    ,8       ,3     ,-1],
           13:['wedge'   ,6       ,3     ,-1],
           14:['pyramid' ,5       ,3     ,-1],
       }

def write_su2(ndims,x,y,z,elzones,ofname):

    npts = len(x)

    f = open(ofname,'w')
    f.write("NDIME= %d\n"%ndims)

    # write fluid elements
    nel_fluid = 0
    for zonename in elzones:
        eltype,elrange,elconn = elzones[zonename]
        eldim = ELTYPES[eltype][2]
        if eldim == ndims:
            nel = len(elconn) / ELTYPES[eltype][1]
            nel_fluid += nel
            print '%d elements of type %s'%(nel,ELTYPES[eltype][0])
            
    f.write("NELEM= %d\n"%nel_fluid)
    iel = 0
    for zonename in elzones:
        eltype,elrange,elconn = elzones[zonename]
        eldim = ELTYPES[eltype][2]
        if eldim == ndims:
            nel = len(elconn) / ELTYPES[eltype][1]
            elconn = elconn.reshape((nel,ELTYPES[eltype][1]))
            for i in range(nel):
                f.write('%d \t '%eltype)
                for j in range(ELTYPES[eltype][1]):
                    f.write('%d \t '%elconn[i,j])
                f.write('%d\n'%iel)
                iel +=1
        
    # write the mesh points
    f.write("NPOIN= %d\n"%npts)
    if ndims == 2:
        for i in range(npts):
            f.write('%g \t %g \t %d\n'%(x[i],y[i],i))
    else:
        for i in range(npts):
            f.write('%g \t %g \t %g \t %d\n'%(x[i],y[i],z[i],i))
    
    # write boundaries
    nmark = 0
    for zonename in elzones:
        eltype,elrange,elconn = elzones[zonename]
        eldim = ELTYPES[eltype][2]
        if eldim == ndims-1:
            nmark +=1

    f.write('NMARK= %d\n'%nmark)
    for zonename in elzones:
        eltype,elrange,elconn = elzones[zonename]
        eldim = ELTYPES[eltype][2]            
        if eldim == ndims-1:
            f.write('MARKER_TAG= %s\n'%zonename)
            nel = len(elconn) / ELTYPES[eltype][1]
            f.write('MARKER_ELEMS= %d\n'%nel)
            elconn = elconn.reshape((nel,ELTYPES[eltype][1]))
            for i in range(nel):
                f.write('%d \t '%eltype)
                for j in range(ELTYPES[eltype][1]):
                    f.write('%d \t '%elconn[i,j])
                f.write('%d\n'%iel)
                iel +=1
        
    f.close()

def convert_uns(ifname,ofname):

    (t,l,p)=CHLone.load(ifname,lksearch=['.'])

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
        elzones[conn[0]] = [eltype,elrange,elconn]

    # write the su2 file
    
    write_su2(ndims,x,y,z,elzones,ofname)

def slice_array(arr,pr):

    i0 = min(pr[0,0],pr[0,1])
    i1 = max(pr[0,0],pr[0,1])
    j0 = min(pr[1,0],pr[1,1])
    j1 = max(pr[1,0],pr[1,1])
    k0 = min(pr[2,0],pr[2,1])
    k1 = max(pr[2,0],pr[2,1])

    return arr[i0:i1+1,j0:j1+1,k0:k1+1]

def convert_struct(ifname,ofname):

    (t,l,p)=CHLone.load(ifname,lksearch=['.'])

    bases = CGU.hasChildType(t,'CGNSBase_t')
    if not len(bases) == 1:
        raise ValueError('there should be 1 base')
    base = bases[0]

    zones = CGU.hasChildType(base,'Zone_t')
    
    X = np.array([])
    Y = np.array([])
    Z = np.array([])

    added_zones = {}
    to_remove = []
    for zone in zones:
        npts = len(X)

        # get the coordinates
        coords = CGU.hasChildType(zone,'GridCoordinates_t')
        if not len(coords) == 1:
            raise ValueError("invalid coordinates")
        x = CGU.getValueByPath(coords[0],'CoordinateX')
        y = CGU.getValueByPath(coords[0],'CoordinateY')
        z = CGU.getValueByPath(coords[0],'CoordinateZ')

        n = 1
        for l in x.shape: n *= l
    
        idx = npts + np.arange(n)
        idx = idx.reshape(x.shape)

        els_0 = idx[:-1,:-1,:-1]
        els_1 = idx[:-1,:-1,:-1]
        els_2 = idx[:-1,:-1,:-1]
        els_3 = idx[:-1,:-1,:-1]
        

        added_zones[zone[0]] = [idx,x.reshape(n),y.reshape(n),z.reshape(n)]

        # get the block interfaces
        gridconnectivities = CGU.hasChildType(zone,'ZoneGridConnectivity_t')
        for gridconnectivity in gridconnectivities:

            # matching interfaces
            interfaces = CGU.hasChildType(gridconnectivity,'GridConnectivity1to1_t')
            for interface in interfaces:
                opp    = CGU.getValueByPath(interface,'.')
                opp    = ''.join(opp)
                pr     = CGU.getValueByPath(interface,'PointRange') - 1
                pr_opp = CGU.getValueByPath(interface,'PointRangeDonor') - 1
                
                periodic = CGU.hasChildType(interface,'GridConnectivityProperty_t')

                
                if opp in added_zones and periodic == None:
                    idx_opp = added_zones[opp][0]
                    i1 = slice_array(idx,pr)
                    i2 = slice_array(idx_opp,pr_opp)
                    i1 = np.squeeze(i1)
                    i2 = np.squeeze(i2)
                    
                    n = 1
                    for l in i1.shape: n *= l

                    ok = False
                    for T in [0,1]:
                        for i in [-1,1]:
                            for j in [-1,1]:
                                idx1 = i1.reshape(n)
                                if T:
                                    idx2 = i2[::i,::j].reshape(n)
                                else:
                                    idx2 = i2.T[::i,::j].reshape(n)

                                dst = (added_zones[zone[0]][1][idx1] - added_zones[opp][1][idx2])**2 + \
                                      (added_zones[zone[0]][2][idx1] - added_zones[opp][2][idx2])**2 + \
                                      (added_zones[zone[0]][3][idx1] - added_zones[opp][3][idx2])**2
                                if (dst.max() < 1e-20):
                                    to_remove.append([idx1,idx2])
                                    ok = True
                                    break                                
                    assert(ok)
            
        
    # if abs(z).max() > 1e-10:
    #     ndims = 3
    # else:
    #     ndims = 2

    # # Get the connectivity
    # elzones = {}
    # conns  = CGU.hasChildType(zone,'Elements_t')
    # for conn in conns:
    #     eltype_cg =  CGU.getValueByPath(conn,'.')[0]
    #     eltype = -1
    #     for el in ELTYPES:
    #         if ELTYPES[el][3] == eltype_cg: eltype = el
    #     elrange   = CGU.getValueByPath(conn,'ElementRange')
    #     elconn    = CGU.getValueByPath(conn,'ElementConnectivity') - 1
    #     elzones[conn[0]] = [eltype,elrange,elconn]

    # # write the su2 file

if __name__=='__main__':

    # Parse the script arguments

    parser=ArgumentParser()

    parser.add_argument("-i","--inputfile", dest='inputfile',
                        default='mesh.cgns',
                        help='cgns input file',metavar='FILE')

    parser.add_argument("-o","--outputfile", dest='outputfile',
                        default='mesh.su2',
                        help='cgns output file',metavar='FILE')

    parser.add_argument("-s","--structured", dest='structured',
                        action='store_true',
                        help='block structured mesh?')

    options = parser.parse_args()

    print 'Input file           :',options.inputfile
    print 'Output file          :',options.outputfile

    if options.structured:
        convert_struct(options.inputfile,options.outputfile)
    else:
        convert_uns(options.inputfile,options.outputfile)
