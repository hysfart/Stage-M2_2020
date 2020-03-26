import os
import sys
from numpy import *
from hotbit import * 
from ase import *
from box.data import data
from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.calculators.abinit import Abinit
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from hotbit.parametrization.atom import KSAllElectron
from hotbit.parametrization.util import IP_EA
from time import asctime

class skf_file:
        def __init__(self, atom1_='Al', atom2_='Al', rcut_ = 5.0, R_ = 2.0):
                """
                On definit ici les différentes variables qui seront utile pour creer les 
                potentiels de la DFTB. On pose que :

                - atom1 / atom2 : 
                les atomes qui vont être nécessaires pour le potntiel, 
                ils sont des objets de classe KSAllElectron(**kargs), 
                ils ont besoins d'un rayon de confinement, et d'être déjà "run()"

                - rcut : 
                c'est le rayon de coupure 

                - name_1 / name_2 : 
                c'est la référence aux noms des elements chimique que l'on prend
                
                """
                
                self.atom1   = atom1_
                self.atom2   = atom2_
                self.rcut    = rcut_
                self.Req     = R_

                Bohr = 0.52917725750691647

                elem1 = atom1_
                r0_cov_1 = data[elem1]['R_cov']/Bohr
                r0_factor_1 = 1.85
                #r0_cov_1 = 75/Bohr
                r0_1 = r0_cov_1 * r0_factor_1
                
                elem2 = atom2_
                r0_cov_2 = data[elem2]['R_cov']/Bohr
                r0_factor_2 = 1.85
                r0_2 = r0_cov_2 * r0_factor_2
                
                # Calculate wave functions for the confined isolated atoms
                self.atom1 = KSAllElectron(elem1, confinement={'mode':'quadratic', 'r0':r0_1})
                self.atom1.run()
                
                self.atom2 = KSAllElectron(elem2, confinement={'mode':'quadratic', 'r0':r0_2})
                self.atom2.run()

                el = self.atom1
                fractional_limit = 0.000001
                r=max( [el.get_wf_range(nl,fractional_limit) for nl in el.get_valence_orbitals()] )
                print("--debug :: La range MAX == {}".format(r))
                
                # Nom des elements
                self.name_1   = elem1
                self.name_2   = elem2
                # Nom de la molécule en étude selon le bon format
                if elem1 == elem2:
                        self.molecule_name = self.name_1+"2"
                else :
                        self.molecule_name = self.name_1+self.name_2
                print("--debug :: La molécule traité ici est : {}".format(self.molecule_name))
                

                # On vérifie ici si les fichiers '.elm' sont présents dans le répertoire
                try:
                        var = 0
                        with open("%s.elm"%self.name_1) : var=1;pass
                        with open("%s.elm"%self.name_2) : var=2;pass
                except IOError:
                        if var ==1:
                                print("Erreur! Le fichier \"{}\" n'a pas pu être ouvert"
                                      .format("%s.elm"%self.name_1))
                        else :
                                print("Erreur! Le fichier \"{}\" n'a pas pu être ouvert"
                                      .format("%s.elm"%self.name_2))

                # Initialisation des titres
                self.title_no_rep = '%s_%s_no_repulsion.par'%(self.name_1,self.name_2)
                self.title_slako  = '%s_%s_slako.pdf'%(self.name_1,self.name_2)
                

        def creatTable(self, rmin_=2, rmax_=14):
                """
                On créé la table Slater-Koster, sur un certain nombres de points, défini
                par N ou plutot par rmin et rmax.
                Les fichiers sont updates.
                """
                table = SlaterKosterTable(self.atom1,self.atom2)
                rmin, rmax, dr = rmin_, rmax_, 0.1
                N = int( (rmax-rmin)/dr )
                
                print('--debug :: There are {} points in the table'.format(N))

                # Remplissage de la table
                table.run(rmin,rmax,N)
                table.write(self.title_no_rep)
                
                
                table.plot(self.title_slako)
        
        def creatRepulsivePotential(self):
                """
                Potentiel repulsif :
                On tente de tracer le potentiel repulsif, à partir du calculateur adéquate.
                   - Le mixer sert à definir les propriétés des calculateurs
                   - Les calculateurs sont décrit en fin de programme
                On fait appel aux fonctions qui servent à trouver les fichiers '.traj'.
                On fit le potentiel selon ces valeurs V(r)

                """
                tab = {"%s%s"%(self.name_1,self.name_2):self.title_no_rep}
                elm = {self.name_1:'%s.elm'%self.name_1,self.name_2:'%s.elm'%self.name_2}
                mixer_ = {'name':'Anderson','mixing_constant':0.2, 'memory':5}
                mixer = mixer_
                calc0 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True)
                calc1 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=False,kpts=(6,6,6))
                calc2 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True,charge=-1.0)
                
                
                # Appel de la classe
                rep = RepulsiveFitting(self.name_2,self.name_2,r_cut=self.rcut)

                """
                Ajout de données - Tracés des points
                -> on decide de tracer les points du potentiel selon : 
                   - un calcul d'énergie de DFTB par le calculateur Hobbit ;
                   - penser à changer le calculateur pour de la DFT (ce serait mieux) quite 
                à inclure MOLCAS directement au projet ou alors récupérer les résultats 
                'manuellement" ;
                   - Chaque calculateur prend des caractéristiques bien spécifique décrit 
                en fin de programme (copie de la doc) ;
                """

                # Ces points sont ABSOLUMENT n'importe quoi !!!
                rep.append_point(weight=1, R=2.9, dvrep=-0.2)
                rep.append_energy_curve(weight=1.0,calc=calc0,traj='Al2.traj',
                                        comment='dimer curve')
                rep.append_scalable_system(weight=1.0,calc=calc1,atoms='Al.traj',
                                           comment='bulk Al')
                
                # Réalisation du fit et de l'output
                rep.fit()
                rep.write_par(self.title_no_rep,
                              filename = '%s_%s_repulsion.par'%(self.name_1,self.name_2))
                rep.plot('%s_%s_repulsion.pdf'%(self.name_1,self.name_2))
                
                print("--debug :: Le fichier {} a bien été créé."
                      .format('%s_%s_repulsion.pdf'%(self.name_1,self.name_2)))
        





        """
        FONCTION DE TEST DE L'ECRITURE
        
        La plupart des fonction ne fonctionne pas ici, mais on servit pour tester les morceaux 
        de code. Sinon, il faut aller voir les fichier "test_0n.py" dans le même répertoire.
        
        """
                
        def trajectoryFile(self, cell_ = [(10, 0, 0), (0, 10, 0), (0, 0, 10)],
                           calculateur_ = EMT(), borne=[1,4.0], N=10):

                """
                atom = Atoms(self.molecule_name, positions=[(0, 0, 0), (0, 0, 3.0)],
                             cell=cell_,
                             pbc=1,
                             calculator=calculateur_)  # use EMT potential
                a1 = Atoms('C', [(0, 0, 0)])
                a2 = Atoms('C', [(0, 0, 3)])
                """
                
                #atoms1=Atoms(symbols='C',positions=[(1,1,1)],cell=(47,48,49),pbc=True,calculator=calculateur_)
                #atoms2=Atoms(symbols='C',positions=[(1,1,1.1)],cell=(17,18,19),pbc=True,calculator=calculateur_)
                
                #bothatoms = atoms1 + atoms2
                """
                cell = atoms1.get_cell()
                traj = Trajectory('%s.traj'%(self.molecule_name),'w', bothatoms)
                for x in linspace(borne[0], borne[1], N):
                        bothatoms.set_positions(bothatoms.positions + 0.25*x)
                        bothatoms.get_potential_energy()
                        print(atoms1.get_potential_energy())
                        traj.write()
                """

                """
                h = Atoms('H', calculator=Abinit(ecut=200, toldfe=0.001))
                h.center(vacuum=3.0)
                e = h.get_potential_energy()
                print(e)
                """
                """
                a0     = 5.43
                dimerH = Atoms([Atom('H', ( 0.7,   0,   0)),
                                Atom('H', (-0.7,   0,   0))],
                                pbc=True)
                b = a0 / 2
                dimerH.set_cell([(10.7,  0,  0),
                                 ( 0, 10.7,  0),
                                 ( 0,  0, 10.7)], scale_atoms=True)
                dimerH.set_celldisp([(2.,2.,2.)])
                
                calc = Abinit(label  ='H2' ,
                              kptopt = 0   ,
                              nkpt   = 1   ,
                              ecut   = 100.0,
                              nstep  = 10  ,
                              toldfe = 1e-6,
                              prtvol = 2    )
                calc.set(diemac = 2.0)
                """
                # - Appel du calculateur par défault
                #calc = Abinit(label  ='H2')
                #calc.set(diemac = 2.0)

                # - Pour renomer les fichier - ne fonctionne pas
                #os.rename("H2.in", "H2_prev.in")
                #os.rename("tbase2_4.in", "H2.in")
                """
                dimerH.set_calculator(calc)
                
                e = dimerH.get_potential_energy()
                print("Energie potentiel : {}".format(e))

                traj = Trajectory('dimerH2.traj', 'w')
                traj.write(dimerH)
                """
                """
                import numpy as np
                a = 4.0  # approximate lattice constant
                b = a / 2
                calc=Abinit(ecut=500,
                                  toldfe=0.001)
                ag = Atoms('C',
                           calculator=calc)  # use pseudo-potential from Abinit
                #cell = ag.get_cell()
                ag.center(vacuum=3.0)
                calc.set(ecut=50)
                """
                """
                import numpy as np
                
                #traj = Trajectory('Ag.traj', 'w')

                for x in np.linspace(0.5, 1.05, 5):
                        dimerH.set_cell(dimerH.cell * x, scale_atoms=True)
                        #print(dimerH.cell * x)
                        print(dimerH.positions)
                        print(dimerH.get_potential_energy())
                        traj.write(dimerH)

                 """


        """
        Fonction d'appel de la class
        """
        def __call__(self):
                self.creatTable()
                self.creatRepulsivePotential()


atom1 = Atoms

g = skf_file()
g()
#g.trajectoryFile()



"""
DOCUMENTATION HOTBIT CALCULATEUR

Hotbit -- density-functional tight-binding calculator
          for atomic simulation environment (ASE).
        Parameters:
        -----------
        parameters:       The directory for parametrization files.
                          * If parameters==None, use HOTBIT_PARAMETERS environment variable.
                          * Parametrizations given by 'elements' and 'tables' keywords
                            override parametrizations in this directory.
        elements:         Files for element data (*.elm).
                          example: {'H':'H_custom.elm','C':'/../C.elm'}
                          * If extension '.elm' is omitted, it is assumed.
                          * Items can also be elements directly: {'H':H} (H is type Element)
                          * If elements==None, use element info from default directory.
                          * If elements['rest']=='default', use default parameters for all other
                            elements than the ones specified. E.g. {'H':'H.elm','rest':'default'}
                            (otherwise all elements present have to be specified explicitly).
        tables:           Files for Slater-Koster tables.
                          example: {'CH':'C_H.par','CC':'C_C.par'}
                          * If extension '.par' is omitted, it is assumed.
                          * If tables==None, use default interactions.
                          * If tables['rest']='default', use default parameters for all other
                            interactions, e.g. {'CH':'C_H.par','rest':'default'}
                          * If tables['AB']==None, ignore interactions for A and B
                            (both chemical and repulsive)
        mixer:            Density mixer.
                          example: {'name':'Anderson','mixing_constant':0.2, 'memory':5}.
        charge:           Total charge for system (-1 means an additional electron)
        width:            Width of Fermi occupation (eV)
        SCC:              Self-Consistent Charge calculation
                          * True for SCC-DFTB, False for DFTB
        kpts:             Number of k-points.
                          * For translational symmetry points are along the directions
                            given by the cell vectors.
                          * For general symmetries, you need to look at the info
                            from the container used
        rs:               * 'kappa': use kappa-points
                          * 'k': use normal k-points. Only for Bravais lattices.
        physical_k        Use physical (realistic) k-points for generally periodic systems.
                          * Ignored with normal translational symmetry
                          * True for physically allowed k-points in periodic symmetries.
        maxiter:          Maximum number of self-consistent iterations
                          * only for SCC-DFTB
        coulomb_solver:   The Coulomb solver object. If None, a DirectCoulomb
                          object will the automatically instantiated.
                          * only for SCC-DFTB
        charge_density:   Shape of the excess charge on each atom. Possibilities
                          are:
                          * 'Gaussian': Use atom centered Gaussians. This is the
                            default.
                          * 'Slater': Slater-type exponentials as used in the
                            original SCC-DFTB scheme.
                          * only for SCC-DFTB
        gamma_cut:        Range for Coulomb interaction if direct summation is
                          selected (coulomb_solver = None).
                          * only for SCC-DFTB
        vdw:              Include van der Waals interactions
        vdw_parameters:   Dictionary containing the parameters for the van-der-Waals
                          interaction for each element.
                          i.e. { el: ( p, R0 ), ... }
                          where *el* is the element name, *p* the polarizability and
                          *R0* the radius where the van-der-Waals interaction starts.
                          Will override whatever read from .elm files.
        txt:              Filename for log-file.
                          * None: standard output
                          * '-': throw output to trash (/null)
        verbose_SCC:      Increase verbosity in SCC iterations.
        internal:         Dictionary for internal variables, some of which are set for
                          stability purposes, some for quick and dirty bug fixes.
                          Use these with caution! (For this reason, for the description
                          of these variables you are forced to look at the source code.)
"""
