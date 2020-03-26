import os
import sys
import numpy as np
from numpy import *
from ase import *
from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.calculators.abinit import Abinit
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from time import asctime


class potential_initialisation :
    def __init__(self, calculateur_ = None, name_ = "Création de fichier *.traj", repertory_ = None):
        
        #print("-- debug :: enter into the class initialisation method")
        self.calculateur = calculateur_
        self.name        = name_
        self.repertory   = repertory_
        
        
    def initiate_parameter(self, atoms_properties_, traj, **kwargs):
        """
        Ici, on place les différentes propriétés des arguments de la fonction d'écriture :
                
        - calculateur_ : objet de calculteur (de type EMT() ou Abinit() par exemple)
        - atoms_properties : argument qui prend un tableau bien précis d'objets comme suit:
            -> un string "Atoms" de comprennant deux atomes au format 'C2' ou CH'
            -> les positions de ces atomes
            -> une cellule des atomes en question
        - calcultateur_properties : argument qui prend un tableau bien précis de 
        valeur comme suit :
            -> un certain nombre d'élements, qui seront ajouté au calculateur, comme 
            propriétés de calcul 
        - name : nom du fichier qui sera écrit
                
        """
        #print("-- ALERT :: CLASS potential_initialisation : write_trajectory_bulk_files")
        atoms_properties = atoms_properties_
        
        """
        Paramétrisation à 2 atomes (exemples) :
        atoms_properties = ['CH',
                            [(pos1)(pos2)], 
                            [(cell1), (cell2), (cell3)], 
                            pbc=False
                            scale_atoms=True]
        """
        nb_atoms = len(atoms_properties[0])
        
        i = 0.0
        traj_atoms = Atoms()
        symbols = ['H','He','Li','Be','B','C','N','O','F',
                   'Ne','Na','Mg','Al','Si','P','S','Cl','Ar']

        print("\n******************************* ")
        print(" Initialisation de(s) atome(s) ")
        print("******************************* ")
        self.print(sys.stdout, title_ = "Paramètres ABINIT : ",**kwargs)

        i  = len(atoms_properties)
        ic = 0
        print(" Les élements présents :")
        while i > 2:
            if atoms_properties[ic] in symbols:
                traj_atoms.append(Atom(atoms_properties[ic],
                                       atoms_properties[len(atoms_properties)-2][ic]))
                print("\t{}, {}".format(atoms_properties[ic],
                                    atoms_properties[len(atoms_properties)-2][ic]))
            else :
                print("La chaine {} ne correspond à aucun atomes connus".format(ic))
            i  -= 1
            ic += 1

        # AFFICHAGE DES DETAILS DE CALCULS
        print("--debug :: Le système se compose avec {}".format(traj_atoms.symbols))

        traj_atoms.set_cell(atoms_properties[len(atoms_properties)-1])
        
        calc = self.calculateur(label = "%s"%traj_atoms.symbols)

        calc.set(**kwargs)
        traj_atoms.set_calculator(calc)

        # Affichage de l'énergie potentiel qui sera écrit dans le fichier ".traj"
        energy = traj_atoms.get_potential_energy()
        print("--debug :: Energie potentiel : {}".format(energy))
        
        if self.name == None or self.name != "%s"%traj_atoms.symbols:
            self.name = "%s"%traj_atoms.symbols
        name_file = self.name+".traj"
            
        # On écrit dans le fichier choisit
        traj.write(traj_atoms)
        print("--debug :: Le fichier {} a bien été créé".format(name_file))


    """
    Fonction d'affichages des variables d'environnement de calcul ABINIT
    """
    def print(self, out_, title_ = "Paramètres ABINIT ", indent_ = " ", **kwargs) :
        if len(title_) :
            out_.write("{}{} {}\n".format(indent_, title_, self.name))
        out_.write("".format())
        for key, value in kwargs.items(): 
            out_.write("\t{} {} = {}\n".format(indent_,key, value)) 

    """
    Fonction d'affichages des variables résultant des calculs ABINITS
        -> print du fichier *.traj
    """
    def print_trajectory_file(self, file_ = None):
        
        if file_ == None:
            try :
                file_ = self.name+".traj"
                tr = Trajectory(file_)
                print(tr.backend)
            except IOError:
                print("--WARNING :: Le fichier n'existe pas ou n'est pas dans le bon répertoire")
        else :
            try :
                tr = Trajectory(file_)
                print(tr.backend)
            except IOError:
                print("--WARNING :: Le fichier n'existe pas ou n'est pas dans le bon répertoire")

        return 0




if __name__ == "__main__" :

    # #################################
    # Trajectory File pour DIMERE LDA
    # ################################# 
    pos    = [(2.2,0.0,0.0),(-2.2,0.0,0.0)]
    cell   = [(20,  0,  0),
              ( 0, 20,  0),
              ( 0,  0, 20)]
    atoms_ = ['Al','Al',
              pos, 
              cell]
    traj = Trajectory("Al2.traj", 'w')
    dimer = potential_initialisation(Abinit)
    
    dimer.initiate_parameter(atoms_, traj,
                            ixc = 1,
                            kptopt = 1   ,
                            kptrlen = 8.0e1 ,
                            kptrlatt = "8 0 0  0 1 0  0 0 1" ,
                            nshiftk = 1,
                            shiftk = "0.5 0 0",
                            occopt = 5,
                            ecut   = 40.0,
                            toldfe = 1e-6,
                            tsmear = 0.05)
    
    pos    = [(2.0,0.0,0.0),(-2.0,0.0,0.0)]
    atoms_ = ['Al','Al',
              pos, 
              cell]
    dimer.initiate_parameter(atoms_, traj,
                            ixc = 1,
                            kptopt = 1   ,
                            kptrlen = 8.0e1 ,
                            kptrlatt = "8 0 0  0 1 0  0 0 1" ,
                            nshiftk = 1,
                            shiftk = "0.5 0 0",
                            occopt = 5,
                            ecut   = 40.0,
                            toldfe = 1e-6,
                            tsmear = 0.05)


    # #################################
    # Trajectory File pour FCC LDA
    # #################################    
    pos    = [(0.0,0.0,0.0)]
    acell  = 7.25
    cell   = [( 0,  acell,  acell),
              ( acell, 0,  acell),
              ( acell,  acell, 0)]
    
    traj = Trajectory("Al-fcc_LDA.traj", 'w') 
    fcc_LDA = potential_initialisation(Abinit)
    
    atoms_ = ['Al',
              pos, 
              cell]
    fcc_LDA.initiate_parameter(atoms_, traj,
                                nband = 6,
                                kptopt = 1,
                                kptrlen = 6.4488e1,
                                kptrlatt = "-12 12 0  0 12 0  0 12 -12",
                                nshiftk = 1,
                                occopt = 5,
                                shiftk ="0 0 0",
                                ecut   = 40.0,
                                toldfe = 1e-6,
                                tsmear = 0.05)

    acell  = 7.0
    cell   = [( 0,  acell,  acell),
              ( acell, 0,  acell),
              ( acell,  acell, 0)]
    atoms_ = ['Al',
              pos, 
              cell]
    fcc_LDA.initiate_parameter(atoms_, traj,
                                nband = 6,
                                kptopt = 1,
                                kptrlen = 6.4488e1,
                                kptrlatt = "-12 12 0  0 12 0  0 12 -12",
                                nshiftk = 1,
                                occopt = 5,
                                shiftk ="0 0 0",
                                ecut   = 40.0,
                                toldfe = 1e-6,
                                tsmear = 0.05)

    acell  = 7.5
    cell   = [( 0,  acell,  acell),
              ( acell, 0,  acell),
              ( acell,  acell, 0)]
    atoms_ = ['Al',
              pos, 
              cell]
    fcc_LDA.initiate_parameter(atoms_, traj,
                                nband = 6,
                                kptopt = 1,
                                kptrlen = 6.4488e1,
                                kptrlatt = "-12 12 0  0 12 0  0 12 -12",
                                nshiftk = 1,
                                occopt = 5,
                                shiftk ="0 0 0",
                                ecut   = 40.0,
                                toldfe = 1e-6,
                                tsmear = 0.05)


    fcc_LDA.print_trajectory_file('Al-fcc_LDA.traj')
    
    sys.exit()

