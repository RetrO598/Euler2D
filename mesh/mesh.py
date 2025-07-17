def convert_to_msh(input_file, output_file):
    with open(input_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # -------------------------------
    # 1. 找到节点/单元/边界数量行
    # -------------------------------
    idx = 0
    while idx < len(lines):
        parts = lines[idx].split()
        if len(parts) == 3 and all(part.isdigit() for part in parts):
            num_nodes, num_cells, num_boundaries = map(int, parts)
            break
        idx += 1
    else:
        raise ValueError("无法找到节点/单元数量行")

    idx += 1  # 进入边界定义部分

    # -------------------------------
    # 2. 读取边界定义
    # -------------------------------
    boundary_info = []
    for _ in range(num_boundaries):
        while lines[idx].startswith("#"):
            idx += 1
        parts = lines[idx].split()
        btype, last_face, last_node = int(parts[0]), int(parts[1]), int(parts[2])
        idx += 1
        name = lines[idx]
        boundary_info.append((btype, last_face, last_node, name))
        idx += 1

    # -------------------------------
    # 3. 读取边界线段
    # -------------------------------
    while not lines[idx].startswith("# boundary faces"):
        idx += 1
    idx += 1  # skip header

    boundary_lines = []
    while not lines[idx].startswith("# coordinates"):
        if not lines[idx].startswith("#"):
            boundary_lines.append(tuple(map(int, lines[idx].split())))
        idx += 1
    idx += 1

    # -------------------------------
    # 4. 读取节点坐标
    # -------------------------------
    coords = []
    while not lines[idx].startswith("# triangles"):
        if not lines[idx].startswith("#"):
            x, y = map(float, lines[idx].split())
            coords.append((x, y))
        idx += 1
    idx += 1

    # -------------------------------
    # 5. 读取三角形单元
    # -------------------------------
    triangles = []
    while idx < len(lines):
        if not lines[idx].startswith("#"):
            triangles.append(tuple(map(int, lines[idx].split())))
        idx += 1

    # -------------------------------
    # 6. 写入 Gmsh 2.2 msh 文件
    # -------------------------------
    with open(output_file, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        num_physical = num_boundaries + 1
        f.write("$PhysicalNames\n")
        f.write(f"{num_physical}\n")
        for i in range(num_boundaries):
            boundary_name = boundary_info[i][3]
            f.write(f'1 {boundary_info[i][0]} "{boundary_name}"\n')
            i += 1
        f.write('2 1 "fluid"\n')
        f.write("$EndPhysicalNames\n")

        # 节点部分
        f.write("$Nodes\n")
        f.write(f"{len(coords)}\n")
        for i, (x, y) in enumerate(coords, 1):
            f.write(f"{i} {x} {y} 0.0\n")
        f.write("$EndNodes\n")

        # 元素部分
        f.write("$Elements\n")
        total_elements = len(boundary_lines) + len(triangles)
        f.write(f"{total_elements}\n")

        eid = 1
        boundary_idx = 0
        bound_index = 0
        for btype, last_face, _, name in boundary_info:
            bound_index += 1
            while boundary_idx < last_face:
                n1, n2 = boundary_lines[boundary_idx]
                f.write(f"{eid} 1 2 {btype} {bound_index} {n1} {n2}\n")
                eid += 1
                boundary_idx += 1

        for tri in triangles:
            f.write(f"{eid} 2 2 1 1 {tri[0]} {tri[1]} {tri[2]}\n")
            eid += 1

        f.write("$EndElements\n")


convert_to_msh("channel.ugr", "channel.msh")
