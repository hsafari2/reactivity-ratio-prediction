from rdkit import *
import pandas as pd
import numpy as np
from abc import ABC, abstractmethod
from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdPartialCharges
from rdkit.Chem import AllChem, rdMolDescriptors





class MinimizedEnergyGenerator:
    @staticmethod
    def generate(mol: Chem.Mol):
        # Make sure there is no added conformers before
        mol.RemoveAllConformers()
        # Add hydrogen
        hydrogenatedMolecule = Chem.AddHs(mol)
        # Generate conformers
        AllChem.EmbedMolecule(hydrogenatedMolecule)


        
        if hydrogenatedMolecule.GetNumConformers() == 0:  # Check if embedding was successful
            # Raise an exception if embedding failed
            
            print(f"There is no Confomer in this molecule!")
            print(str(Chem.MolToSmiles(mol)))

        # MMFF Optimization
        AllChem.MMFFOptimizeMolecule(hydrogenatedMolecule)
         
        return hydrogenatedMolecule






#parent
class getGroup(ABC):
    @staticmethod
    @abstractmethod
    def getSpecificGroup(self, mol: Chem.Mol):
        pass


#child 
class getVinyls(getGroup):
    
    @staticmethod
    def getSpecificGroup (mol: Chem.Mol):
        bonds = mol.GetBonds()
        vinyls = []
        for i in bonds:
            if  (i.GetBeginAtom().GetSymbol() == 'C' and i.GetEndAtom().GetSymbol() == 'C' and (i.GetBondTypeAsDouble() == 2 or i.GetBondTypeAsDouble() == 3 )):
                vinyls.append(i)

        if vinyls == []:
            print('vinyl is empty for this molecule')
            print(str(Chem.MolToSmiles(mol)))
            
        return vinyls


#parent
class Descriptor(ABC):
    @staticmethod
    @abstractmethod
    def calculate(self, mol: Chem.Mol):
        pass


#child
class MolecularWeightCalculator(Descriptor):
    @staticmethod
    def calculate(mol: Chem.Mol):
        
        AllChem.AddHs(mol)

        # Calculate the molecular weight
        MW = Descriptors.ExactMolWt(mol)

        
        return MW
        

#child
class VinylLoopIndicator(Descriptor):
    @staticmethod
    def calculate(mol: Chem.Mol):
        vinylGroup = getVinyls.getSpecificGroup(mol)
        for vinyl  in vinylGroup:
            if vinyl.IsInRing():
                return True
        return False
        


#child        
class PolarChargeCarbonVinyl(Descriptor):
    
    @staticmethod
    def calculate(mol: Chem.Mol):

        # Minimize energy
        minimizedEnergyMolecule = MinimizedEnergyGenerator.generate(mol)

        # Compute Gasteiger charges
        AllChem.ComputeGasteigerCharges(minimizedEnergyMolecule)

        #get vinyl Group
         
        atomCarbonBond = getVinyls.getSpecificGroup(minimizedEnergyMolecule)[0]
        



        #carbon 1
        carbon1 = atomCarbonBond.GetBeginAtom()
        #carbon2
        carbon2 = atomCarbonBond.GetEndAtom()

        carbonChargeList = []

        charge1 = float(carbon1.GetProp('_GasteigerCharge'))
        carbonChargeList.append(charge1)

        charge2 = float(carbon2.GetProp('_GasteigerCharge'))
        carbonChargeList.append( charge2)

        return carbonChargeList


        
        
#child        
class DipoleMomentCalculator(Descriptor):
    
    @staticmethod
    def calculate(mol: Chem.Mol):
        # Minimize energy for the molecule
        minimizedEnergyMolecule = MinimizedEnergyGenerator.generate(mol)

        #compute the charge over the molecule 
        rdPartialCharges.ComputeGasteigerCharges(minimizedEnergyMolecule)

        #initialize the dipole for the molecule
        dipole_moment = [0.0, 0.0, 0.0]

        #go thorgh the all atoms in the molecule
        for atom in minimizedEnergyMolecule.GetAtoms():
            #compute charge in the atom
            charge = float(atom.GetProp('_GasteigerCharge'))
            
            #get vector of the position of the atoms
            pos = minimizedEnergyMolecule.GetConformer().GetAtomPosition(atom.GetIdx())
            
            for i in range(3):
                dipole_moment[i] += charge * pos[i]
        
        # Calculate magnitude of the dipole moment
        magnitude = sum(d**2 for d in dipole_moment)**0.5
        return magnitude


    

#child        
class MolecularVolumeCalculator(Descriptor):
    
    @staticmethod
    def calculate(mol: Chem.Mol) -> float:
        
        # Minimize energy for the molecule
        minimizedEnergyMolecule = MinimizedEnergyGenerator.generate(mol)

        MV = AllChem.ComputeMolVolume(minimizedEnergyMolecule)
        
        return MV

        
class MolecularSurfaceAreaCalculator(Descriptor):
    
    @staticmethod
    def calculate(mol: Chem.Mol) -> float:
        # Minimize energy for the molecule
        minimizedEnergyMolecule = MinimizedEnergyGenerator.generate(mol)

        # Assuming that we are using solvent-accessible surface area (SASA)
        mol_surface_area = rdMolDescriptors.CalcTPSA(minimizedEnergyMolecule)

        return mol_surface_area





