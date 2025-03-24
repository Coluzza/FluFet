import MDAnalysis as mda
import numpy as np
import glob
import sys
import re

def extract_brush_atom_types(topology_file):
    """Extracts atom types of brush molecules from the LAMMPS topology file."""
    with open(topology_file, "r") as f:
        lines = f.readlines()

    atom_section = False
    brush_atom_types = set()

    for line in lines:
        if not atom_section and line.startswith("Atoms"):  
            print("✅ Found 'Atoms' section.")  # Debugging print
            atom_section = True
            continue
        
        
        if atom_section:
            data = line.split()
            if len(data) < 7:
                continue  # Skip incomplete lines
            
            mol_id = int(data[1])  # Molecule ID (Column 2)
            atom_type = int(data[2])  # Atom Type (Column 3)
            if mol_id != 1:  # Exclude protein (Mol-ID 1)
                brush_atom_types.add(atom_type)
        if atom_section and (line.strip() == "" or re.match(r"^\D", line)):  
            break  # Stop at next non-numeric line

    if not brush_atom_types:
        print("❌ No brush atoms identified. Check your topology.")
        sys.exit(1)

    brush_atom_types = sorted(brush_atom_types)
    print(f"🟢 Identified brush atom types: {brush_atom_types}")
    return brush_atom_types

def compute_brush_profile(topology_file, xtc_files, atom_types, bins=100):
    """Computes the brush profile along the Z-axis averaged over the trajectory."""

    try:
        u = mda.Universe(topology_file, xtc_files, topology_format="DATA")
    except ValueError as e:
        print(f"❌ Error loading topology: {e}")
        sys.exit(1)

    # Select brush atoms
    brush_atoms = u.select_atoms(f"type {' '.join(map(str, atom_types))}")

    if len(brush_atoms) == 0:
        print("❌ No brush atoms found. Check the topology.")
        sys.exit(1)

    # Get min/max Z to define histogram range
    z_min, z_max = brush_atoms.positions[:, 2].min(), brush_atoms.positions[:, 2].max()
    z_bins = np.linspace(z_min, z_max, bins + 1)
    z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

    # Initialize histogram accumulator
    density_profile = np.zeros(bins)

    # Process each frame
    total_frames = 0
    for ts in u.trajectory:
        z_positions = brush_atoms.positions[:, 2]  
        hist, _ = np.histogram(z_positions, bins=z_bins)
        density_profile += hist
        total_frames += 1

    # Normalize by total frames
    density_profile /= total_frames

    return z_centers, density_profile

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python Brush_profile.py <lammps_topology.data> [bins]")
        sys.exit(1)

    topology_file = sys.argv[1]  
    bins = int(sys.argv[2]) if len(sys.argv) > 2 else 100  

    xtc_files = sorted(glob.glob("7_flow.lammpstrj_*.xtc"))  

    # Extract brush atom types
    atom_types = extract_brush_atom_types(topology_file)

    # Compute brush profile
    z_values, profile = compute_brush_profile(topology_file, xtc_files, atom_types, bins)

    # Save results
    np.savetxt("brush_profile.dat", np.column_stack((z_values/10, profile)), header="Z (Å)  Density", fmt="%.5f")

    print("✅ Brush profile calculated and saved as 'brush_profile.dat'.")
