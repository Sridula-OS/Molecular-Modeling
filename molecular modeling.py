from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors

# Define your molecule in SMILES format
smiles = "NC1=CC=CC=C1"
# Create an RDKit molecule object
mol = Chem.MolFromSmiles(smiles)
# Generate 2D coordinates for the molecule
AllChem.Compute2DCoords(mol)
# Draw the molecule
img = Draw.MolToImage(mol)
# Save or display the image
img.save("molecule_2d.png")
img.show()

#No. of atoms
print('Number of atoms in Morphine is: ',mol.GetNumAtoms())
#Mol wt
print('Molecular weight of Morphine: ',Descriptors.MolWt(mol))
#logp
print('LogP of Morphine: ',Descriptors.MolLogP(mol))
#solubility
print('Molar refractivity of Morphine: ',Descriptors.MolMR(mol))
#bond details
i=0
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol(), 
    atom.GetAtomicNum(),atom.GetHybridization()
    , mol.GetBondWithIdx(i).GetBondType())
    i+=1
#rotatble bonds
print('Number of rotatable bonds of Morphine: ',Descriptors.NumRotatableBonds(mol))
#Calculate hydrogen bonds involving oxygen, hydrogen, and nitrogen
hbonds = AllChem.CalcNumHBD(mol)
#Print the result
print(f"Number of hydrogen bonds involving oxygen, hydrogen, and nitrogen(Hydrogen Bond Donors): {hbonds}")
# Count the number of nitrogen and oxygen atoms
num_nitrogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)  # Atomic number for nitrogen is 7
num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)  # Atomic number for oxygen is 8
print("Number of oxygen and nitrogen atoms(Hydrogen Bond Acceptors):", num_oxygen_atoms,' and ',num_nitrogen_atoms)

print('''\nLipinski's rule states that an orally active drug follows the given criteria:
    No more than 5 hydrogen bond donors (the total number of nitrogen–hydrogen and oxygen–hydrogen bonds)
    No more than 10 hydrogen bond acceptors (all nitrogen or oxygen atoms)
    A molecular mass less than 500 daltons
    A calculated octanol-water partition coefficient (Clog P) that does not exceed 5

Hence morphine is an orally active drug''')

