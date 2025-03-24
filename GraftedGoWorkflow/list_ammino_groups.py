import sys

# List of residues with amino (-NH2) groups in side chains
AMINO_GROUP_RESIDUES = {"LYS", "ARG", "ASN", "GLN", "HIS"}

def extract_residues_with_amino_group(pdb_file, selected_chains):
    """Extracts renumbered residue IDs for residues with an amino (-NH2) group in specific chains."""
    found_residues = []
    residue_map = {}  # Maps (chain, original residue ID) to a global sequential ID
    current_id = 1  # Global residue numbering

    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("ATOM"):
                chain_id = line[21]  # Extract chain ID
                res_name = line[17:20].strip()  # Extract residue name
                original_res_id = line[22:26].strip()  # Extract original residue ID

                res_key = (chain_id, original_res_id)  # Unique key with chain

                # Assign sequential numbering if not already mapped
                if res_key not in residue_map:
                    residue_map[res_key] = current_id
                    current_id += 1
                
                # Check if residue is in selected chains and has an amino group
                if chain_id in selected_chains and res_name in AMINO_GROUP_RESIDUES:
                    found_residues.append(residue_map[res_key])

    found_residues = sorted(set(found_residues))  # Remove duplicates and sort
    return found_residues

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python list_amino_residues_sequential.py <pdb_file> <chain1> <chain2> ...")
        sys.exit(1)

    pdb_filename = sys.argv[1]
    selected_chains = sys.argv[2:]  # Take chain IDs from command-line arguments
    residues = extract_residues_with_amino_group(pdb_filename, selected_chains)

    if residues:
        print(len(residues))  # Print total number of identified residues
        for res_id in residues:
            print(res_id)  # Print each relative residue ID from global count
    else:
        print("0")  # If no residues found, print 0
