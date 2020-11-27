import numpy as np
from ase import Atoms, Atom
from ase.io import write
from ase.data import chemical_symbols
from sklearn import manifold, datasets
from sklearn.decomposition import PCA
from openbabel import pybel as pb
from openbabel import openbabel as ob
from ase.visualize import view

from rdkit.Chem.rdmolfiles import MolFromSmiles,MolToMolFile
from rdkit import Chem
from rdkit.Chem import AllChem

import argparse



def __ase2xyz__(atoms):
    """
    Prepare a XYZ string from an ASE atoms object.
    """
    # Implementation detail: If PBC should be implemented, the
    # write to xyz needs to be changed to include cell etc.
    atoms.pbc=False
    if any(atoms.get_pbc()):
        raise RuntimeError("Detected PBCs. Not supported (yet)!")
    num_atoms = len(atoms)
    types = atoms.get_chemical_symbols()
    all_atoms = zip(types, atoms.get_positions())
    a_str = str(num_atoms) + "\n" + "\n"
    for atom in all_atoms:
        a_str += atom[0] + " " + " ".join([str(x) for x in atom[1]]) + "\n"
    return a_str


def convert_ase2pybel(atoms):
    """
    Convert an ASE atoms object to pybel (openBabel) molecule.
    The ordering of the Atoms is identical.

    Parameters
    ----------
    atoms : ase.Atoms
        The ASE atoms object

    Returns
    -------
    pymol :
        The pybel molecule.
    """
    atoms.pbc=False
    a_str = __ase2xyz__(atoms)
    pymol = pb.readstring("xyz", a_str)

    return pymol

def convert_ase2rdkit(atoms, removeHs=False):
    """
    Convert an ASE atoms object to rdkit molecule.
    The ordering of the Atoms is identical.


    Important: Implemented only for clusters, not PBC!
    rdkit does not keep xyz coordinates, therefore
    a backconversion is not possible yet.

    Parameters
    ----------
    atoms : ase.Atoms
        The ASE atoms object
    removeHs : Bool
        If True, remove all H atoms from molecule.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The rdkit molecule object.
    """
    a_str = __ase2xyz__(atoms)
    pymol = pb.readstring("xyz", a_str)
    mol = pymol.write("mol")
    mol = Chem.MolFromMolBlock(mol, removeHs=removeHs)
    return mol

def pybel2ase(mol):  
    asemol = Atoms()
    species=[chemical_symbols[atm.atomicnum] for atm in mol.atoms]
    pos=np.asarray([atm.coords for atm in mol.atoms])
    pca = PCA(n_components=3)
    posnew=pca.fit_transform(pos)
    #posnew[:,2]=0.0
    atoms = Atoms(species, positions=posnew)
    sys_size = np.ptp(atoms.positions,axis=0)
    atoms.pbc=True
    atoms.cell = sys_size + 10
    atoms.center()
    
    return atoms

def rdkit2ase(m):
    pos = m.GetConformer().GetPositions()
    natoms = m.GetNumAtoms()
    species = [m.GetAtomWithIdx(j).GetSymbol() for j in range(natoms)]                
#Get the principal axes and realign the molecule
    pca = PCA(n_components=3)
    pca.fit(pos)
    posnew=pca.transform(pos)        
#Set the z to 0.0       
    #posnew[:,2]=0.0
    atoms = Atoms(species, positions=posnew)  
    sys_size = np.ptp(atoms.positions,axis=0)
    atoms.pbc=True
    atoms.cell = sys_size + 10
    atoms.center()
        
    return atoms


def pybel_opt(smile,steps):
    obconversion = ob.OBConversion()
    obconversion.SetInFormat('smi')
    obmol = ob.OBMol()
    obconversion.ReadString(obmol,smile)
    
    pbmol = pb.Molecule(obmol)
    pbmol.make3D(forcefield="uff", steps=50)
    
    pbmol.localopt(forcefield="gaff", steps=200)
    pbmol.localopt(forcefield="mmff94", steps=100)
    
    f_f = pb._forcefields["uff"]
    f_f.Setup(pbmol.OBMol)
    f_f.ConjugateGradients(steps, 1.0e-9)
    f_f.GetCoordinates(pbmol.OBMol)
    #print(f_f.Energy(), f_f.GetUnit())
    
    return pybel2ase(pbmol)

def try_rdkit(smile):
    test = Chem.MolFromSmiles(smile)
    test = Chem.AddHs(test)
    
    #test=convert_ase2rdkit(pybel_opt(smile,10))
    
    return AllChem.EmbedMolecule(test, maxAttempts=10, randomSeed=42)
    
def rdkit_opt(smile,steps):
    m = Chem.MolFromSmiles(smile)
    m = Chem.AddHs(m)
    
    AllChem.EmbedMolecule(m, maxAttempts=20, randomSeed=42)
    AllChem.UFFOptimizeMolecule(m,maxIters=steps)
    
    return rdkit2ase(m)

def mol_from_smiles(smile,steps=10000):
    if try_rdkit(smile)==0:
        #print("Going for rdkit")
        return rdkit_opt(smile,steps)
    else:
        #print("Going for pybel")
        return pybel_opt(smile,steps)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("SMILE", help="SMILE string")
    parser.add_argument("-o", "--output",type=str,default="output",help="Name for the output PDB file")
    args = parser.parse_args()
    
    # Let's convert the SMILE to an ASE object
    mol=mol_from_smiles(args.SMILE)
    write(args.output+".pdb",mol)

