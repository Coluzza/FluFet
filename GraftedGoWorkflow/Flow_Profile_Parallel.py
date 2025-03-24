import numpy as np
import sys
import multiprocessing as mp

def extract_selected_chunks(input_file, output_file, selected_chunks):
    """
    Extracts only the specified X-chunks from the velocity file.
    """
    selected_chunks = set(map(int, selected_chunks))  # Convert to set for fast lookup
    extracted_lines = []
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#") or len(line.split()) < 8:
                continue  # Skip headers and malformed lines
            parts = line.split()
            chunk_id = int(parts[1])  # Extract chunk index
            if chunk_id in selected_chunks:
                extracted_lines.append(line)
    
    if not extracted_lines:
        print("❌ No matching chunks found.")
        sys.exit(1)
    
    with open(output_file, 'w') as f:
        f.writelines(extracted_lines)
    
    print(f"✅ Extracted {len(extracted_lines)} lines matching chunks {selected_chunks}.")

def process_chunk(chunk_lines):
    """Processes a subset of velocity data to compute average velocity profiles."""
    z_data = {}
    zy_data = {}
    z_counts = {}
    zy_counts = {}

    for line in chunk_lines:
        parts = line.split()
        z = float(parts[3])  # Z coordinate
        y = float(parts[4])  # Y coordinate
        vx, vy, vz = map(float, (parts[5], parts[6], parts[7]))

        if z not in z_data:
            z_data[z] = np.array([vx, vy, vz])
            z_counts[z] = 1
        else:
            z_data[z] += np.array([vx, vy, vz])
            z_counts[z] += 1

        zy_key = (z, y)
        if zy_key not in zy_data:
            zy_data[zy_key] = np.array([vx, vy, vz])
            zy_counts[zy_key] = 1
        else:
            zy_data[zy_key] += np.array([vx, vy, vz])
            zy_counts[zy_key] += 1

    return z_data, z_counts, zy_data, zy_counts

def parallel_process_velocity(file, num_threads):
    """Reads velocity data, splits it into chunks, and processes in parallel."""
    with open(file, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith("#") and len(line.split()) >= 8]
    
    num_lines = len(lines)
    chunk_size = num_lines // num_threads
    
    pool = mp.Pool(num_threads)
    chunks = [lines[i * chunk_size:(i + 1) * chunk_size] for i in range(num_threads)]
    results = pool.map(process_chunk, chunks)
    pool.close()
    pool.join()
    
    # Combine results from all threads
    z_data, z_counts, zy_data, zy_counts = {}, {}, {}, {}
    for z_chunk, zc_chunk, zy_chunk, zyc_chunk in results:
        for z in z_chunk:
            if z not in z_data:
                z_data[z] = z_chunk[z]
                z_counts[z] = zc_chunk[z]
            else:
                z_data[z] += z_chunk[z]
                z_counts[z] += zc_chunk[z]
        
        for zy in zy_chunk:
            if zy not in zy_data:
                zy_data[zy] = zy_chunk[zy]
                zy_counts[zy] = zyc_chunk[zy]
            else:
                zy_data[zy] += zy_chunk[zy]
                zy_counts[zy] += zyc_chunk[zy]
    
    # Normalize results
    for z in z_data:
        z_data[z] /= z_counts[z]
    for zy in zy_data:
        zy_data[zy] /= zy_counts[zy]
    
    return z_data, zy_data

def save_velocity_profile(data, output_file, is_zy=False):
    """Saves velocity profile to a file."""
    with open(output_file, "w") as f:
        if is_zy:
            f.write("# Z   Y   Vx   Vy   Vz\n")
            for (z, y), (vx, vy, vz) in sorted(data.items()):
                f.write(f"{z:.6f} {y:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
        else:
            f.write("# Z   Vx   Vy   Vz\n")
            for z, (vx, vy, vz) in sorted(data.items()):
                f.write(f"{z:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
    print(f"✅ Saved velocity profile to '{output_file}'")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python Flow_Profile_Parallel.py <velocity_file> <num_threads> <chunk1> <chunk2> ...")
        sys.exit(1)

    velocity_file = sys.argv[1]
    num_threads = int(sys.argv[2])
    selected_chunks = sys.argv[3:]
    
    filtered_file = "filtered_velocity_data.txt"
    extract_selected_chunks(velocity_file, filtered_file, selected_chunks)
    
    print(f"⏳ Processing velocity file: {filtered_file} with {num_threads} threads...")
    z_data, zy_data = parallel_process_velocity(filtered_file, num_threads)

    save_velocity_profile(z_data, "velocity_profile_z.dat")
    save_velocity_profile(zy_data, "velocity_profile_zy.dat", is_zy=True)

    print("✅ Processing complete. Data ready for Gnuplot visualization.")
