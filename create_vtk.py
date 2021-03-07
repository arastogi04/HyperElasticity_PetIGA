from igakit.io import PetIGA,VTK
from numpy import linspace
import glob
#uniform = lambda U: linspace(U[0], U[-1], 100)
j = 0
nrb = PetIGA().read("geometry"+str(j)+".dat")

for infile in glob.glob("disp*.dat"):
    sol = PetIGA().read_vec(infile,nrb)
    outfile = infile.split(".")[0] + ".vtk"
    VTK().write(outfile,
                nrb,
                fields=sol,
                #sampler=uniform,
                vectors={'solution':[0,1]})
