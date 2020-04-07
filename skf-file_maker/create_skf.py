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
from fitting_2 import *

from param_Abinit_maker import *


class skf_file:
        def __init__(self, atom1_='Al', atom2_='Al', rcut_ = 5.4, R_ = 2.0):
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
                self.atom1.write_unl('Al.elm')
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
                calc1 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True,gamma_cut=2.0)
                calc2 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True,charge=-1)
                calc3 = Hotbit(txt='-',elements=elm,mixer=mixer,tables=tab,SCC=True,charge=+1) 
               

                # save : [(-0.7296118704,0.8012224477,-1.1920228795),(0.5749930354,-1.4822803034,-0.2596875043),(0.1546188350,0.6810578558,1.4517103838)]
                pbc = False
                pos    = [(-0.7296118704,0.8012224477,-1.1920228795),(0.5749930354,-1.4822803034,-0.2596875043),(0.1546188350,0.6810578558,1.4517103838)]
                acell  = 50
                cell   = [( acell, 0, 0),
                        ( 0, acell, 0), 
                        ( 0, 0, acell)]

                c3_LDA = abinitio_potential()

                atoms_ = ['Al','Al','Al', pos, cell, pbc]
                c3_LDA.initial_parameters(atoms_)
                param = c3_LDA.atoms
              
                # Appel de la classe
                rep = RepulsiveFitting(self.name_2,self.name_2,r_cut=self.rcut)

                from ase import Atoms

                t = Atoms(symbols='Al10',
                        positions=[(-1.9186181624, -1.2461903716,  2.0949120187),
                        (1.7189795173,  0.3209741318, -2.7757087630),
                        (3.4905846244, -0.3068711577, -0.6988022955),
                        (1.8411920351,  1.9299250914, -0.5489679385),
                        (2.2864531542,  0.4070430211,  1.7067325251),
                        (1.4718330447, -2.2539373690,  1.0382629850),
                        (-2.2450606096,  1.3551110704,  1.1606376673),
                        (-3.8084873252, -0.6509295625,  0.1149179121),
                        (-1.1906453574, -0.7975179805, -0.5007932389),
                        (-2.6462309211,  1.2423931268, -1.5911908723)],
                        cell=cell)

                # Points pour le cluster Al10
                rep.append_homogeneous_cluster(weight=0.01,calc=calc1,atoms=t,comment='cluster Al_10 curve',label="cluster Al10")
                #rep.append_homogeneous_cluster(weight=1,calc=calc2,atoms=t,comment='cluster Al_10 curve',label="cluster Al10 -")
                #rep.append_homogeneous_cluster(weight=1,calc=calc3,atoms=t,comment='cluster Al_10 curve',label="cluster Al10 +")


                # Points pour le trimer de Al
                #rep.append_homogeneous_cluster(weight=1,calc=calc1,atoms=param,comment='cluster Al_3 curve',label="cluster Al3")
                #rep.append_homogeneous_cluster(1.0,calc1,'Al-3_LDA.traj',comment='cluster Al_3 curve',label="cluster Al3")


                e = 1.0
                pos    = [(-0.7296118704*e,0.8012224477*e,-1.1920228795*e),(0.5749930354*e,-1.4822803034*e,-0.2596875043*e),(0.1546188350*e,0.6810578558*e,1.4517103838*e)]
                acell  = 50
                cell   = [( acell, 0, 0),
                        ( 0, acell, 0),
                        ( 0, 0, acell)]
                q = Atoms('Al3', positions=pos, cell=cell)
                rep.append_homogeneous_cluster(weight=0.05,calc=calc1,atoms=q,comment='cluster Al_3 curve')

                e = 0.85
                pos    = [(-0.7296118704*e,0.8012224477*e,-1.1920228795*e),(0.5749930354*e,-1.4822803034*e,-0.2596875043*e),(0.1546188350*e,0.6810578558*e,1.4517103838*e)]
                q = Atoms('Al3', positions=pos, cell=cell)
                #rep.append_homogeneous_cluster(weight=0.5,calc=calc1,atoms=q,comment='cluster Al_3 curve')

                e = 0.95
                pos    = [(-0.7296118704*e,0.8012224477*e,-1.1920228795*e),(0.5749930354*e,-1.4822803034*e,-0.2596875043*e),(0.1546188350*e,0.6810578558*e,1.4517103838*e)]
                q = Atoms('Al3', positions=pos, cell=cell)
                #rep.append_homogeneous_cluster(weight=0.5,calc=calc1,atoms=q,comment='cluster Al_3 curve')




                t = Atoms(symbols='Al3',
                        positions=[(-1.7296118704,1.8012224477,-2.1920228795),(1.5749930354,-2.4822803034,-1.2596875043),(1.1546188350,1.6810578558,2.4517103838)],cell=cell)
                #rep.append_homogeneous_cluster(weight=0.3,calc=calc1,atoms=t,comment='cluster Al_3 curve',label="cluster Al3")
                
                
                
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
                o = ['FCC', -52.21372200113491, -52.464288211025725, -52.695972314978036, -52.130185729801426, -51.99540489548358,
                        'CC', -52.262294090116185, -52.13562106779265, -51.921780044220846, -52.40378588439806, -52.52957128932968, -51.574129460197526]
                l = ['FCC', 3.985, 3.7, 3.5, 4.1, 4.3, 'CC', 3.13, 3.25, 3.5, 3.0, 2.9]
                
                # Les points de FCC
                rep.append_point(weight=0.7, R=2.863, dvrep=(o[1]+16.5017),label="FCC")
                rep.append_point(weight=0.7, R=3.985, dvrep=(o[1]+16.5017)/12,label="FCC")
                #rep.append_point(weight=0.99, R=3.985*np.sqrt(2), dvrep=(o[1]+16.5017)/18,label="FCC")
                
                #rep.append_point(weight=1, R=l[2]/2, dvrep=(o[2]+16.6177)/12)
                #rep.append_point(weight=1, R=l[3]/2, dvrep=(o[3]+16.7734)/12)
                #rep.append_point(weight=1, R=l[4]/2, dvrep=(o[4]+16.4747)/12)
                #rep.append_point(weight=1, R=l[5]/2, dvrep=(o[5]+16.4474)/12)
                

                # Points suplémentaire AU FCC bulk
                """
                d = ['FCC', -52.21372200113491, -52.695972314978036, -53.58063130613739, -52.02871140432845, -54.692280873780156, -55.31389106044733, -54.139118607676934, 'CC']
                g = ['FCC', 3.985, 3.5, 3.0, 4.25, 2.5, 2.25, 2.75, 'CC']
                rep.append_point(weight=1, R=g[1]/1, dvrep=(d[1]+16.5017)/12)
                rep.append_point(weight=1, R=g[2]/1, dvrep=(d[2]+16.7734)/12)
                rep.append_point(weight=1, R=g[3]/1, dvrep=(d[3]+17.1653)/12)
                rep.append_point(weight=1, R=g[5]/1, dvrep=(d[5]+19.3642)/12)
                rep.append_point(weight=1, R=g[6]/1, dvrep=(d[6]+20.6727)/12)
                rep.append_point(weight=1, R=g[7]/1, dvrep=(d[7]+18.4393)/12)
                """

                # Les points CC
                rep.append_point(weight=0.7, R=2.802, dvrep=(o[7]+16.5195),label="CC")
                rep.append_point(weight=0.7, R=3.185, dvrep=(o[7]+16.5195)/8,label="CC")
                rep.append_point(weight=0.7, R=3.185*np.sqrt(2), dvrep=(o[7]+16.5195)/14,label="CC")


                # Courbe pour le dimer Al2
                rep.append_energy_curve(weight=0.1,calc=calc0,traj='Al2.traj',comment='dimer curve',label="Dimere Al2")


                rep.append_dimer(weight=0.1,calc=calc2,R=2.0,comment='dimer curve',label="Dimere Al2-")                
                
                #rep.append_dimer(weight=1.0,calc=calc3,R=2.0,comment='dimer curve',label="Dimere Al2+")
                
                print(rep.deriv) 
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
                #self.creatTable()
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
