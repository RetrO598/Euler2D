import sys

def convert_ugr_to_su2(input_file, output_file):
    """
    Converts a UGR grid file to a SU2 grid file.
    """
    with open(input_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # 1. Find the line with number of nodes, cells, and boundaries
    idx = 0
    while idx < len(lines):
        parts = lines[idx].split()
        if len(parts) == 3 and all(part.isdigit() for part in parts):
            num_nodes, num_cells, num_boundaries = map(int, parts)
            break
        idx += 1
    else:
        raise ValueError("Could not find the line with node/cell counts.")

    idx += 1  # Move to the boundary definitions

    # 2. Read boundary definitions
    boundary_info = []
    for _ in range(num_boundaries):
        while lines[idx].startswith("#"):
            idx += 1
        parts = lines[idx].split()
        btype, last_face, last_node = int(parts[0]), int(parts[1]), int(parts[2])
        idx += 1
        name = lines[idx]
        boundary_info.append({"btype": btype, "last_face": last_face, "name": name})
        idx += 1

    # 3. Read boundary line segments
    while not lines[idx].startswith("# boundary faces"):
        idx += 1
    idx += 1  # skip header

    boundary_lines = []
    while not lines[idx].startswith("# coordinates"):
        if not lines[idx].startswith("#"):
            boundary_lines.append(tuple(map(int, lines[idx].split())))
        idx += 1
    idx += 1

    # 4. Read node coordinates
    coords = []
    while not lines[idx].startswith("# triangles"):
        if not lines[idx].startswith("#"):
            x, y = map(float, lines[idx].split())
            coords.append((x, y))
        idx += 1
    idx += 1

    # 5. Read triangle elements
    triangles = []
    while idx < len(lines):
        if not lines[idx].startswith("#"):
            triangles.append(tuple(map(int, lines[idx].split())))
        idx += 1
    
    if len(coords) != num_nodes:
        print(f"Warning: Number of nodes in header ({num_nodes}) does not match number of coordinates found ({len(coords)}).")
    
    if len(triangles) != num_cells:
        print(f"Warning: Number of cells in header ({num_cells}) does not match number of triangles found ({len(triangles)}).")


    # 6. Write SU2 file
    with open(output_file, "w") as f:
        # Dimension
        f.write("NDIME= 2\n")

        # Nodes
        f.write(f"NPOIN= {len(coords)}\n")
        for i, (x, y) in enumerate(coords):
            f.write(f"{x} {y} {i}\n")

        # Elements
        f.write(f"NELEM= {len(triangles)}\n")
        for i, tri in enumerate(triangles):
            # SU2 uses 0-based indexing for nodes
            n1, n2, n3 = tri
            f.write(f"5 {n1-1} {n2-1} {n3-1} {i}\n")

        # Boundary Markers
        f.write(f"NMARK= {len(boundary_info)}\n")
        start_face_idx = 0
        for info in boundary_info:
            f.write(f"MARKER_TAG= {info['name']}\n")
            end_face_idx = info['last_face']
            num_marker_elems = end_face_idx - start_face_idx
            f.write(f"MARKER_ELEMS= {num_marker_elems}\n")
            for i in range(start_face_idx, end_face_idx):
                n1, n2 = boundary_lines[i]
                # SU2 uses 0-based indexing
                f.write(f"3 {n1-1} {n2-1}\n")
            start_face_idx = end_face_idx

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ugr_to_su2.py <input.ugr> <output.su2>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        convert_ugr_to_su2(input_file, output_file)
        print(f"Successfully converted '{input_file}' to '{output_file}'.")
    except (ValueError, IndexError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
