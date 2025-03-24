import numpy as np
import sys

def read_velocity_file(filename):
    """
    Reads the velocity file and processes it to compute:
    1. The velocity profile along the Z-axis (averaged over X chunks and time)
    2. The velocity profile in the Z-Y plane (averaged over X chunks and time)
    """

    z_data = {}  # Dictionary to accumulate velocity data along Z-axis
    zy_data = {}  # Dictionary to accumulate velocity data along Z-Y plane
    z_counts = {}  # Counter for averaging Z data
    zy_counts = {}  # Counter for averaging Z-Y data

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or len(line.split()) < 8:
                continue  # Skip headers and malformed lines

            parts = line.split()
            z = float(parts[3])  # Z coordinate
            y = float(parts[4])  # Y coordinate
            vx, vy, vz = map(float, (parts[5], parts[6], parts[7]))

            # Average along the Z-axis
            if z not in z_data:
                z_data[z] = np.array([vx, vy, vz])
                z_counts[z] = 1
            else:
                z_data[z] += np.array([vx, vy, vz])
                z_counts[z] += 1

            # Average over the Z-Y plane
            zy_key = (z, y)
            if zy_key not in zy_data:
                zy_data[zy_key] = np.array([vx, vy, vz])
                zy_counts[zy_key] = 1
            else:
                zy_data[zy_key] += np.array([vx, vy, vz])
                zy_counts[zy_key] += 1

    # Compute final averages
    for z in z_data:
        z_data[z] /= z_counts[z]

    for zy_key in zy_data:
        zy_data[zy_key] /= zy_counts[zy_key]

    return z_data, zy_data

def save_velocity_profile_z(z_data, output_file="velocity_profile_z.dat"):
    """ Saves velocity profile along Z in gnuplot format. """
    with open(output_file, "w") as f:
        f.write("# Z   Vx   Vy   Vz\n")
        for z in sorted(z_data.keys()):
            vx, vy, vz = z_data[z]
            f.write(f"{z:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
    print(f"✅ Velocity profile along Z saved as '{output_file}'")

def save_velocity_profile_zy(zy_data, output_file="velocity_profile_zy.dat"):
    """ Saves velocity surface plot over the Z-Y plane in gnuplot format. """
    with open(output_file, "w") as f:
        f.write("# Z   Y   Vx   Vy   Vz\n")
        for (z, y) in sorted(zy_data.keys()):
            vx, vy, vz = zy_data[(z, y)]
            f.write(f"{z:.6f} {y:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
    print(f"✅ Velocity profile over Z-Y plane saved as '{output_file}'")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_velocity.py <velocity_file>")
        sys.exit(1)

    velocity_file = sys.argv[1]

    print(f"⏳ Processing velocity file: {velocity_file}")
    z_data, zy_data = read_velocity_file(velocity_file)

    if not z_data:
        print("❌ No valid velocity data found in the file.")
        sys.exit(1)

    save_velocity_profile_z(z_data)
    save_velocity_profile_zy(zy_data)

    print("✅ Processing complete. Data ready for gnuplot visualization.")
