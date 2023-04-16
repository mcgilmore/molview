from rdkit.Chem import rdchem

class Molecule(rdchem.Mol):
    def __init__(self, mol):
        rdchem.Mol.__init__(self, mol)
        self.m1_pve = None
        self.m2_pve = None
        self.m3_pve = None
        self.m1_nve = None
        self.m2_nve = None
        self.m3_nve = None
