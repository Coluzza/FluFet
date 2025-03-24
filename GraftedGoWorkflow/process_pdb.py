from Bio import PDB
import numpy as np

# Input and Output PDB filenames
pdb_file = "8xln_ordered.pdb"
output_pdb_file = "revised_8xln.pdb"

# Initialize PDB parser
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

# Extract protein chains only (removing water, ligands, and other non-protein atoms)
protein_chains = [
    chain for model in structure for chain in model if any(res.id[0] == " " for res in chain)
]

# Sort chains alphabetically by original ID
protein_chains.sort(key=lambda x: x.get_id())

# Rename chains sequentially from 'A'
chain_names = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
for i, chain in enumerate(protein_chains):
    chain.id = chain_names[i]

# Extract atomic coordinates
atoms = [atom.get_coord() for chain in protein_chains for res in chain for atom in res]
atoms = np.array(atoms)

# Compute principal axes using PCA
center = np.mean(atoms, axis=0)
atoms_centered = atoms - center
cov_matrix = np.cov(atoms_centered.T)
eigvals, eigvecs = np.linalg.eig(cov_matrix)

# Ensure the largest eigenvector aligns with Z-axis
sorted_indices = np.argsort(eigvals)[::-1]
principal_axis = eigvecs[:, sorted_indices[0]]

# Create rotation matrix to align principal axis with Z-axis
rotation_matrix = np.eye(3)
rotation_matrix[:, 2] = -principal_axis  # Align the principal axis with Z
rotation_matrix[:, 0] = np.cross([0, 1, 0], principal_axis)  # Create orthogonal X-axis
rotation_matrix[:, 1] = np.cross(principal_axis, rotation_matrix[:, 0])  # Create Y-axis

# Normalize rotation matrix
rotation_matrix[:, 0] /= np.linalg.norm(rotation_matrix[:, 0])
rotation_matrix[:, 1] /= np.linalg.norm(rotation_matrix[:, 1])
rotation_matrix[:, 2] /= np.linalg.norm(rotation_matrix[:, 2])

# Rotate all atoms
atoms_rotated = np.dot(atoms_centered, rotation_matrix)

# Create new structure with sorted chains
new_structure = PDB.Structure.Structure("protein_reordered")
new_model = PDB.Model.Model(0)
new_structure.add(new_model)

# Add sorted chains in order
atom_index = 1  # Atom serial number
residue_index = 1  # Residue sequential numbering

for i, chain in enumerate(protein_chains):
    new_chain = PDB.Chain.Chain(chain_names[i])
    new_model.add(new_chain)

    for res in chain:
        # Assign a new residue number
        new_res_id = (' ', residue_index, ' ')
        new_res = PDB.Residue.Residue(new_res_id, res.resname, res.segid)
        new_chain.add(new_res)

        for atom in res:
            # Assign a new atom number
            new_atom = PDB.Atom.Atom(
                atom.name,
                atoms_rotated[atom_index - 1] + center,  # Adjust coordinates back
                atom.bfactor,
                atom.occupancy,
                atom.altloc,
                atom.fullname,
                atom_index,  # Sequential atom numbering
                atom.element
            )
            new_res.add(new_atom)
            atom_index += 1  # Increment atom counter

        residue_index += 1  # Increment residue counter

# Save the modified PDB file with forced order and renumbering
io = PDB.PDBIO()
io.set_structure(new_structure)
io.save(output_pdb_file)

print(f"✅ Processed PDB file saved as {output_pdb_file}")
