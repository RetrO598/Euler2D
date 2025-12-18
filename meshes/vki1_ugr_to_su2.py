import sys
import re
from collections import defaultdict
import numpy as np

def convert_vki1_ugr_to_su2(input_file, output_file):
    """
    Converts the VKI-1 UGR grid file to a SU2 grid file,
    merging periodic boundaries into a single master/slave pair.
    """
    with open(input_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    idx = 0
    while idx < len(lines):
        parts = lines[idx].split()
        if len(parts) == 3 and all(part.isdigit() for part in parts):
            num_nodes, num_cells, num_boundaries = map(int, parts)
            break
        idx += 1
    else:
        raise ValueError("Could not find the line with node/cell counts.")
    idx += 1

    boundary_definitions = []
    bnode_start_offset = 0
    bface_start_offset = 0
    for _ in range(num_boundaries):
        while lines[idx].startswith("#"):
            idx += 1
        parts = lines[idx].split()
        btype, last_face_ugr, last_node_ugr = map(int, parts)
        idx += 1
        name = lines[idx].strip()
        idx += 1
        
        is_periodic = "periodic" in name.lower() or btype == 700
        
        info = {
            "name": name,
            "is_periodic": is_periodic,
            "num_faces": last_face_ugr - bface_start_offset if not is_periodic else 0,
            "num_nodes_or_pairs": last_node_ugr - bnode_start_offset,
        }
        boundary_definitions.append(info)

        if not is_periodic:
            bface_start_offset = last_face_ugr
        bnode_start_offset = last_node_ugr
    
    while not lines[idx].startswith("# boundary faces"):
        idx += 1
    idx += 1

    boundary_entities_raw = []
    while not lines[idx].startswith("# coordinates"):
        if not lines[idx].startswith("#"):
            boundary_entities_raw.append(tuple(map(int, lines[idx].split())))
        idx += 1
    idx += 1
    
    coords = []
    while not lines[idx].startswith("# triangles") and not lines[idx].startswith("# tetrahedra"):
        if not lines[idx].startswith("#"):
            parts = re.findall(r"-?\d+\.?\d*e?[+-]?\d*", lines[idx])
            if len(parts) >= 2:
                x, y = map(float, parts[:2])
                coords.append((x, y))
        idx += 1
    idx += 1

    elements = []
    elem_type_code = 5 # Triangles
    while idx < len(lines) and not lines[idx].startswith("#"):
        elements.append(tuple(map(int, lines[idx].split())))
        idx += 1

    # Process boundaries
    periodic_data = {}
    non_periodic_data = {}
    
    entity_idx = 0
    for info in boundary_definitions:
        name = info["name"]
        if info["is_periodic"]:
            num_pairs = info["num_nodes_or_pairs"]
            # The first column in the UGR periodic section is the node on the boundary
            pairs = boundary_entities_raw[entity_idx : entity_idx + num_pairs]
            periodic_data[name] = [p[0] for p in pairs]
            entity_idx += num_pairs
        else:
            num_faces = info["num_faces"]
            faces = boundary_entities_raw[entity_idx : entity_idx + num_faces]
            non_periodic_data[name] = faces
            entity_idx += num_faces

    # Build adjacency list for line reconstruction
    adj = defaultdict(set)
    for elem in elements:
        for i in range(len(elem)):
            n1 = elem[i]
            n2 = elem[(i + 1) % len(elem)]
            adj[n1].add(n2)
            adj[n2].add(n1)
    
    # --- Merge periodic boundaries into one master/slave pair ---
    # Assumption: Odd-numbered periodic markers are masters, even-numbered are slaves.
    master_markers = ["Periodic_1", "Periodic_3", "Periodic_5"]
    slave_markers = ["Periodic_2", "Periodic_4", "Periodic_6"]

    master_nodes_all = set()
    for name in master_markers:
        if name in periodic_data:
            master_nodes_all.update(periodic_data[name])

    slave_nodes_all = set()
    for name in slave_markers:
        if name in periodic_data:
            slave_nodes_all.update(periodic_data[name])

    def reconstruct_lines(node_set):
        lines_recon = []
        for node in node_set:
            for neighbor in adj[node]:
                if neighbor in node_set and node < neighbor:
                    lines_recon.append((node, neighbor))
        return lines_recon

    master_lines = reconstruct_lines(master_nodes_all)
    slave_lines = reconstruct_lines(slave_nodes_all)

    # --- Write SU2 file ---
    with open(output_file, "w") as f:
        f.write("NDIME= 2\n")
        f.write(f"NPOIN= {len(coords)}\n")
        for i, (x, y) in enumerate(coords):
            f.write(f"{x:.16f} {y:.16f} {i}\n")

        f.write(f"NELEM= {len(elements)}\n")
        for i, elem in enumerate(elements):
            nodes_str = " ".join(map(str, [n - 1 for n in elem]))
            f.write(f"{elem_type_code} {nodes_str} {i}\n")
        
        num_markers = len(non_periodic_data) + 2 # +2 for periodic_master and periodic_slave
        f.write(f"NMARK= {num_markers}\n")

        for name, faces in non_periodic_data.items():
            f.write(f"MARKER_TAG= {name}\n")
            f.write(f"MARKER_ELEMS= {len(faces)}\n")
            for n1, n2 in faces:
                f.write(f"3 {n1-1} {n2-1}\n")
        
        # Write merged periodic master marker
        f.write(f"MARKER_TAG= periodic_master\n")
        f.write(f"MARKER_ELEMS= {len(master_lines)}\n")
        for n1, n2 in master_lines:
            f.write(f"3 {n1-1} {n2-1}\n")

        # Write merged periodic slave marker
        f.write(f"MARKER_TAG= periodic_slave\n")
        f.write(f"MARKER_ELEMS= {len(slave_lines)}\n")
        for n1, n2 in slave_lines:
            f.write(f"3 {n1-1} {n2-1}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vki1_ugr_to_su2.py <input.ugr> <output.su2>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        convert_vki1_ugr_to_su2(input_file, output_file)
        print(f"Successfully converted '{input_file}' to '{output_file}'.")
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
