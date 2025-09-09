from rdkit import Chem
from rdkit.Chem import AllChem

# Base class for Fingerprint
class Fingerprint:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)

    def generate(self):
        raise NotImplementedError("Subclasses must implement this method")

# Subclass for Morgan Fingerprint
class MorganFingerprint(Fingerprint):
    def __init__(self, smiles, radius, nBits):
        super().__init__(smiles)
        self.radius = radius
        self.nBits = nBits

    def generate(self):
        # Generate the Morgan fingerprint with specified radius and bit length
        return AllChem.GetMorganFingerprintAsBitVect(self.mol, self.radius, nBits=self.nBits)

