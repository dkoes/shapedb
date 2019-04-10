#!/usr/bin/python
import sys, os, numpy, math
from optparse import OptionParser
from pybel import *
import openbabel
from numpy.linalg import eig,det

    # Create the inertia tensor
    # m_i = mass of atom i
    # (x_i, y_i, z_i) = pos of atom i
    # Ixx = sum(m_i*(y_i^2+z_i^2)); Iyy = sum(m_i*(x_i^2+z_i^2)); Izz = sum(m_i*(x_i^2+y_i^2))
    # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
    # Ixz = Izx = -1*sum(m_i*x_i*z_i)
    # Iyz = Izy = -1*sum(m_i*y_i*z_i)
    #this is shamelessly stolen from the mdanalysis package

#calculate the moment of inertia for axis C on both sides
def calcSplitMoments(coordinates,C):
    above = 0
    below = 0
    A = (C+1)%3
    B = (C+2)%3
    for pt in coordinates:
        val = pt[A]*pt[A] + pt[B]*pt[B]
        if pt[C] < 0:
            below += val
        else:
            above += val
    return (below,above)

def calcI(masses, coordinates):
    values = zip(masses, coordinates)
    # Create the inertia tensor
    # m_i = mass of atom i
    # (x_i, y_i, z_i) = pos of atom i
    # Ixx = sum(m_i*(y_i^2+z_i^2)); Iyy = sum(m_i*(x_i^2+z_i^2)); Izz = sum(m_i*(x_i^2+y_i^2))
    # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
    # Ixz = Izx = -1*sum(m_i*x_i*z_i)
    # Iyz = Izy = -1*sum(m_i*y_i*z_i)
    Ixx = reduce(lambda t, a: t + a[0] * (a[1][1] * a[1][1] + a[1][2] * a[1][2]), values, 0.)
    Iyy = reduce(lambda t, a: t + a[0] * (a[1][0] * a[1][0] + a[1][2] * a[1][2]), values, 0.)
    Izz = reduce(lambda t, a: t + a[0] * (a[1][0] * a[1][0] + a[1][1] * a[1][1]), values, 0.)
    Ixy = Iyx = -1 * reduce(lambda t, a: t + a[0] * a[1][0] * a[1][1], values, 0.)
    Ixz = Izx = -1 * reduce(lambda t, a: t + a[0] * a[1][0] * a[1][2], values, 0.)
    Iyz = Izy = -1 * reduce(lambda t, a: t + a[0] * a[1][1] * a[1][2], values, 0.)

    return numpy.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])
#read in molecules and align them to their moments of inertia



def alignmol(mol):
    masses = []    
    coordinates = [ a.coords for a in mol.atoms]
    totalMass = 0
    for atom in mol:
        m = 1
        if options.useweights:
            m = atom.atomicmass
        masses.append(m)
        totalMass += m
    
    coordinates = numpy.array(coordinates)
    masses = numpy.array(masses)
    trans = -numpy.sum(coordinates * masses[:, numpy.newaxis], axis=0) / totalMass
    mol.OBMol.Translate(openbabel.vector3(*trans))
    coordinates += trans
    #this is shamelessly stolen from the mdanalysis package
    inertia = calcI(masses, coordinates)

    eigenval, eigenvec = eig(inertia)
        # Sort
    indices = numpy.argsort(eigenval)
        # Return transposed in more logical form. See Issue 33.
    principalAxes = eigenvec[:,indices].T
    if det(principalAxes) < 0:
        principalAxes = -principalAxes #remove reflection !!Is this correct?? Kabasch wiki says negate just one row, but this is Kabasch...

    xrot = openbabel.double_array([1,0,0,0,-1,0,0,0,-1])
    yrot = openbabel.double_array([-1,0,0,0,1,0,0,0,-1])
    zrot = openbabel.double_array([-1,0,0,0,-1,0,0,0,1])
    ident = openbabel.double_array([1,0,0,0,1,0,0,0,1])
    
    rotmat = openbabel.double_array(principalAxes.flatten())
    mol.OBMol.Rotate(rotmat)

    coordinates = [ a.coords for a in mol.atoms] #should be modified
    (a,b) = calcSplitMoments(coordinates, 1)
    if(a < b):
        mol.OBMol.Rotate(xrot)
#    print a,b
    (a,b) = calcSplitMoments(coordinates, 0)
    if(a < b):
        mol.OBMol.Rotate(yrot)  
#    print a,b
    return



parser = OptionParser()
parser.add_option("-w", "--weights", help="use weights for inertial mass", action="store_true", dest="useweights")

(options, args) = parser.parse_args()

if len(args) != 2:
    print "Need input and output files"
    sys.exit(-1)
    
ofs = Outputfile("sdf", args[1])
for mol in readfile("sdf",args[0]):
    assert mol.OBMol.NumConformers() == 1
    alignmol(mol)
    ofs.write(mol)
