#import numpy as np

#from ase import Atoms
#from ase.io.trajectory import Trajectory
#from ase.calculators.emt import EMT


"""
a = 4.0  # approximate lattice constant
b = a / 2
ag = Atoms('C2',positions=[(0,0,0),(0,0.5,0)],
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=EMT())  # use EMT potential
cell = ag.get_cell()
traj = Trajectory('Ag.traj', 'w')
for x in np.linspace(0.95, 2.05, 10):
    ag.set_cell(cell * x, scale_atoms=True)
    ag.get_potential_energy()
    traj.write(ag)
"""
import numpy
import ase.io as io
import sys
import os

def main(argv):
    tr = io.Trajectory(argv)
    return tr

if __name__=="__main__" :

    tr = main(sys.argv[1])
    print("Longueur de la chaine :: {}".format(len(tr)))
    print("Périodic Boundary Condition (pbe) :: {}".format(tr.pbc))
    print("Backend   :: \n{}".format(tr.backend))
    print("Masses    :: \n{}".format(tr.masses))
    #print("Energies  :: \n{}".format(tr.backend.energies))
    print("Backend   :: \n{}".format(tr.backend))
    print("Calculator.forces :: \n{}".format(tr.backend.calculator.forces))
    print("Backend   :: \n{}".format(tr.backend))
    print("Positions :: \n{}".format(tr.backend.positions))
    print("Forces    :: \n{}".format(tr.backend.calculator.forces))
    





"""
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
configs = read('Ag.traj@0:10')  # read 5 configurations
# Extract volumes and energies:
volumes = [ag.get_volume() for ag in configs]
energies = [ag.get_potential_energy() for ag in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Ag-eos.png')
"""

