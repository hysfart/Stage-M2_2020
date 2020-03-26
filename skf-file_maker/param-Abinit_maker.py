import os
import sys
import numpy
import json

from ase import __version__
import ase.io.jsonio as jsonio



# INITIALISATION DES PARAMETRES ABINIT
energy    = -1.26252524021298E+02
energies  = -1.2887971288E+00

forces    = [( 0.22920362060196, -0.00000000000000, -0.00000000000000),
             (-0.22920362060196, -0.00000000000000, -0.00000000000000)]
cell      = (2.0000000000E+01, 2.0000000000E+01, 2.0000000000E+01)
numbers   = 2
pbc       = [False,False,False]
positions = [(-2.32837971779600, 0.00000000000000, 0.00000000000000),
             (2.32837971779600, 0.00000000000000, 0.00000000000000)]

abinit_parameters = []


def write_file():
    f_ = open("tmp.traj","a")
    f_.write("{}".format(spositions))

    return 0











if __name__ == '__main__':
    
    #write_file()    
    jsonio.write_json("tmp.traj",positions)

    
    
    sys.exit()




