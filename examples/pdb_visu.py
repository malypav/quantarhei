# -*- coding: utf-8 -*-
import numpy

from quantarhei import *

from quantarhei.builders.pdb import PDBFile
from quantarhei.builders import pdb
from quantarhei.models.chlorophylls import ChlorophyllA, ChlorophyllB

from quantarhei import Aggregate
from quantarhei import energy_units
import matplotlib.pyplot as plt

Pi=3.14159
Temperature=15.0

writefiles = True
center_only = False
disorder_averaged = True
weakTEJ = False
#
# Read a PDB file
#
file = PDBFile("1rwt-cut.pdb")
#file = PDBFile("3eni.pdb")
print("Loaded", file.linecount, "lines")

#
# Bacteriochlorophyll model will be used to extract molecules from PDB
#
cl_model = ChlorophyllA(model_type="PDB")
molecules_A = file.get_Molecules(model=cl_model)
cl_model = ChlorophyllB(model_type="PDB")
molecules_B = file.get_Molecules(model=cl_model)
names = []
molecules=[]
for m in molecules_A:
    molecules.append(m)
for m in molecules_B:
    molecules.append(m)


#Optional renaming
#naming_map = {"A371":"BChl1", "A372":"BChl2",
#              "A373":"BChl3", "A374":"BChl4", "A375":"BChl5", 
#              "A376":"BChl6", "A377":"BChl7", "A378":"BChl8"}
#
#for name in names:
#    for m in molecules:
#        if m.name == name:
#            m.set_name(naming_map[name])

#sorting by name    
molecules.sort(key=lambda molecule: molecule.name)

for m in molecules:
    names.append(m.name)

#names.sort()
print(names)

#add coordinates of Nitrogen atoms
for m in molecules:
    for line in m.data:
                if pdb.line_matches(line, by_atmName="NA"):
                    xyz = pdb.line_xyz(line)
                    m.NAxyz=xyz
                if pdb.line_matches(line, by_atmName="NB"):
                    xyz = pdb.line_xyz(line)
                    m.NBxyz=xyz
                if pdb.line_matches(line, by_atmName="NC"):
                    xyz = pdb.line_xyz(line)
                    m.NCxyz=xyz
                if pdb.line_matches(line, by_atmName="ND"):
                    xyz = pdb.line_xyz(line)
                    m.NDxyz=xyz

#
#
    
#
# Create an new aggregate of the Bacteriochlorophylls without BChl8
#
agg = Aggregate(name="LHCII", molecules=molecules)


# Setting site energies according to literature
energies_lit={
        "A602": 15160.0,
        "A603": 15283.0,
        "A610": 15073.0,
        "A611": 15115.0,
        "A612": 15097.0,
        "A613": 15763.0,
        "A614": 15721.0,
        "A604": 15890.0,
        "A601": 15175.0,
        "A608": 15264.0,
        "A609": 15460.0,
        "A605": 15679.0,
        "A606": 15851.0,
        "A607": 15712.0
        }

with energy_units("1/cm"):
    for na in names:
        m=agg.get_Molecule_by_name(na)
        m.set_energy(1,energies_lit[na])

#with energy_units("1/cm"):
#    m = agg.get_Molecule_by_name("A602")
#    m.set_energy(1, 15160.0)
#    m = agg.get_Molecule_by_name("A603")
#    m.set_energy(1, 15283.0)
#    m = agg.get_Molecule_by_name("A610")
#    m.set_energy(1, 15073.0)    
#    m = agg.get_Molecule_by_name("A611")
#    m.set_energy(1, 15115.0)
#    m = agg.get_Molecule_by_name("A612")
#    m.set_energy(1, 15097.0)
#    m = agg.get_Molecule_by_name("A613")
#    m.set_energy(1, 15763.0)
#    m = agg.get_Molecule_by_name("A614")
#    m.set_energy(1, 15721.0)
#    m = agg.get_Molecule_by_name("A604")
#    m.set_energy(1, 15890.0)
#    m = agg.get_Molecule_by_name("A601")
#    m.set_energy(1, 15175.0)
#    m = agg.get_Molecule_by_name("A608")
#    m.set_energy(1, 15264.0)
#    m = agg.get_Molecule_by_name("A609")
#    m.set_energy(1, 15460.0)
#    m = agg.get_Molecule_by_name("A605")
#    m.set_energy(1, 15679.0)
#    m = agg.get_Molecule_by_name("A606")
#    m.set_energy(1, 15851.0)
#    m = agg.get_Molecule_by_name("A607")
#    m.set_energy(1, 15712.0)


#delete dipole moment of 611,612 for no coupling (but still considered in equilibrium)
if weakTEJ:
    m=agg.get_Molecule_by_name("A610")
    m.set_dipole(0,1,m.get_dipole(0,1)/10.0)
    m=agg.get_Molecule_by_name("A611")
    m.set_dipole(0,1,m.get_dipole(0,1)/10.0)
    m=agg.get_Molecule_by_name("A612")
    m.set_dipole(0,1,m.get_dipole(0,1)/10.0)
    
# Set resonance coupling by dipole-dipole method
# !!Screening Value to agree with Novoderezhkin!!! - could also be included in the molecule dipole in Debye!!
agg.set_coupling_by_dipole_dipole(epsr=1.21/1.25)


#
# Build the aggregate
#
agg.build()

numpy.set_printoptions(precision=2, linewidth=100,
                           formatter={'all':lambda x: "%8.1f" % x})
#
# Now we can start simulations
#
H = agg.get_Hamiltonian()

sb_reference = BasisReferenceOperator(H.dim,
                                      name="site basis reference")

manag = Manager()
manag.warn_about_basis_change = False 

with energy_units("1/cm"):
    print("\nExcited state Hamiltonian (energies in 1/cm):\n")
#    for i in range(1, H.dim):
#        for j in range(i,H.dim):
#            print(i,j,":",H.data[i,j])
    with eigenbasis_of(sb_reference):
        print(H.data[1:,1:])
        f=open(agg.name+'_H.dat','w')
        f.write(str(H.data[1:,1:]))
        f.close()

with eigenbasis_of(H):
    with energy_units("1/cm"):
        print(H.data.diagonal())
#H.undiagonalize()
#with energy_units("1/cm"):
#    print(H.data.diagonal())

#print(H.SS)
        
with eigenbasis_of(H):
            #rho=agg.get_DensityMatrix(condition_type="thermal_excited_state", temperature=300.0) 
            rho=agg.get_DensityMatrix(condition_type="impulsive_excitation",temperature=300.0)
        #    print(rho.data.diagonal())    
        #with eigenbasis_of(sb_reference):
        #    print(rho.data.diagonal())
        
with eigenbasis_of(sb_reference):
        pops=rho.data.diagonal().real
        with energy_units("1/cm"):
            Hrho=0.5*(numpy.dot(H.data,rho.data.real)+numpy.dot(rho.data.real,H.data)).diagonal()
            Hrho=(Hrho-min(Hrho[1:]))/(max(Hrho)-min(Hrho[1:]))

#Gaussian disorder averaging
if disorder_averaged:
    sigma=110.0
    ensemble=1000
else:
    sigma=0.0
    ensemble=1
Hrho_av=numpy.zeros(H.dim)
pops_av=numpy.zeros(H.dim)
TEpops=numpy.zeros(ensemble)

for din in range(0,ensemble):
    with energy_units("1/cm"):
        for na in names:
            m=agg.get_Molecule_by_name(na)
            m.set_energy(1,energies_lit[na]+numpy.random.normal(0,sigma))
    agg.build()
    H = agg.get_Hamiltonian()
    #is T in K?
    with eigenbasis_of(H):
        #rho=agg.get_DensityMatrix(condition_type="thermal_excited_state", temperature=300.0) 
        rho=agg.get_DensityMatrix(condition_type="thermal_excited_state",temperature=Temperature)
#        with energy_units("1/cm"):
#            print('XX:' + str(H.data.diagonal()))
    
    with eigenbasis_of(sb_reference):
            pops=rho.data.diagonal().real
            with energy_units("1/cm"):
                Hrho=0.5*(numpy.dot(H.data,rho.data.real)+numpy.dot(rho.data.real,H.data)).diagonal()
                Hrho=(Hrho-min(Hrho[1:]))/(max(Hrho)-min(Hrho[1:]))
    Hrho_av+=Hrho
    pops_av+=pops
    TEpops[din]=pops[10]+pops[11]+pops[12]
Hrho_av=Hrho_av/ensemble
pops_av=pops_av/ensemble    



#plotting the histogram of Terminal Emitter population
weights=numpy.ones_like(TEpops)/float(len(TEpops))*100
plt.hist(TEpops, bins=15,weights=weights)
plt.show()
#writing the density files for VMD

def gauss(x,x0=numpy.array([0,0,0]),alpha=0.025):
    #return pow(alpha/Pi,1.5)*numpy.exp(-alpha*numpy.dot(x-x0,x-x0))
    return numpy.exp(-alpha*numpy.dot(x-x0,x-x0))
#print(gauss(x=[1,1,1]))


def dens(x,n,alpha=0.025):
    mymol=agg.get_Molecule_by_name(names[n])
    if center_only:
        return(gauss(x=x,x0=mymol.position,alpha=alpha))
    return (gauss(x=x,x0=mymol.NAxyz,alpha=alpha)+gauss(x=x,x0=mymol.NBxyz,alpha=alpha)+
           gauss(x=x,x0=mymol.NCxyz,alpha=alpha)+gauss(x=x,x0=mymol.NDxyz,alpha=alpha))





if writefiles:
    xl=60
    yl=40
    zl=40
    x0=-30
    y0=0
    z0=85
    
    fdens=open(agg.name+'_dens_eq4K-dis110.dx', 'w')
    fedens=open(agg.name+'_edens_eq4K-dis110.dx', 'w')
    
    fdens.write('# Data calculated by QuantaRhei\n') 
    fdens.write("object 1 class gridpositions counts %s %s %s \n" %(xl,yl,zl))
    fdens.write("origin %s %s %s \n" %(x0,y0,z0) )
    fdens.write('delta 1 0 0\n')
    fdens.write('delta 0 1 0\n')
    fdens.write('delta 0 0 1\n')
    fdens.write("object 2 class gridconnections counts %s %s %s\n"%(xl,yl,zl))
    fdens.write("object 3 class array type double rank 0 items %s data follows\n\n"%(xl*yl*zl))
    
    fedens.write('# Data calculated by QuantaRhei\n') 
    fedens.write("object 1 class gridpositions counts %s %s %s \n" %(xl,yl,zl))
    fedens.write("origin %s %s %s \n" %(x0,y0,z0) )
    fedens.write('delta 1 0 0\n')
    fedens.write('delta 0 1 0\n')
    fedens.write('delta 0 0 1\n')
    fedens.write("object 2 class gridconnections counts %s %s %s\n"%(xl,yl,zl))
    fedens.write("object 3 class array type double rank 0 items %s data follows\n\n"%(xl*yl*zl))
    
    
    #        print(H.data.diagonal())
    #        print(Hrho)
    #1/2 above for fun:)
        
    #for n in range(0,H.dim-1):
    #    print(pops[n+1])
        
    for x in range(x0,x0+xl):
        for y in range(y0,y0+yl):
            for z in range(z0,z0+zl):
                totdens=0.0
                totedens=0.0
                for n in range(0,H.dim-1):
                    locdens=dens(x=[x,y,z],n=n,alpha=0.025)
                    locedens=dens(x=[x,y,z],n=n,alpha=0.01)
                    totdens+=locdens*pops_av[n+1]
                    totedens+=locedens*Hrho_av[n+1]
                fdens.write(str(totdens)+' ')
                fedens.write(str(totedens)+' ')
                if (z-z0) % 3 == 0:
                    fdens.write('\n')
                    fedens.write('\n')
    
    fdens.write('\nobject "density (element N) [A^-3]" class field')           
    fdens.close()
    fedens.write('\nobject "density (element N) [A^-3]" class field')           
    fedens.close()
#now print the density


    
    