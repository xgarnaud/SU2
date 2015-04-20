import sys
import CGNS
import CGNS.PAT.cgnsutils as CGU
import CHLone
from argparse import ArgumentParser

#              name       npoints ndims  cgnsid
ELTYPES = {3 :['line'    ,2       ,1     ,3 ],
           5 :['triangle',3       ,2     ,5 ],
           9 :['quad'    ,4       ,2     ,7 ],
           10:['tetra'   ,4       ,3     ,-1],
           12:['hexa'    ,8       ,3     ,-1],
           13:['wedge'   ,6       ,3     ,-1],
           14:['pyramid' ,5       ,3     ,-1],
       }

def convert(ifname,ofname):

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

    npts = len(x)
    print '%d points'%npts

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
    

if __name__=='__main__':

    # Parse the script arguments

    parser=ArgumentParser()

    parser.add_argument("-i","--inputfile", dest='inputfile',
                        default='mesh.cgns',
                        help='cgns input file',metavar='FILE')

    parser.add_argument("-o","--outputfile", dest='outputfile',
                        default='mesh.su2',
                        help='cgns output file',metavar='FILE')
    
    options = parser.parse_args()

    print 'Input file           :',options.inputfile
    print 'Output file          :',options.outputfile

    convert(options.inputfile,options.outputfile)
