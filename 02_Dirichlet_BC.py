import numpy as np
import gmsh
import os

# Build the path to the .msh file in the Input subfolder
input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Input')
mesh_file = os.path.join(input_dir, 'MESH_FILE_4x.msh')  # Replace with the name of the .msh file

# Initialize Gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output for debugging
gmsh.open(mesh_file)

print(f"Mesh loaded from: {mesh_file}")

# Check the physical groups defined in the file
physical_groups = gmsh.model.getPhysicalGroups()

# Definition of Dirichlet boundary conditions
face_DC = [2]  # Replace with the index of the desired face
DC_BC = np.array([
    [1],
    [1],
    [1]
])  # Each column refers to a face, and each row to the node coordinate (0 = not constrained, 1 = constrained)

DC_def = np.array([
    [0],
    [0],
    [0]
])  # Imposed displacements along x, y, z for each face

# Initialize empty arrays
D_nodes = []
D_BC_x = []
D_BC_y = []
D_BC_z = []
D_coord_x = []
D_coord_y = []
D_coord_z = []

# Iterate over the specified faces
for i in range(len(face_DC)):
    # Find the nodes belonging to the face using Gmsh
    node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, face_DC[i])

    # Check if nodes were found
    if len(node_tags) == 0:
        print(f"Error: no node found for the physical group: {face_DC[i]}.")
        continue

    # Add nodes and boundary conditions
    D_nodes.extend(node_tags)
    D_BC_x.extend([DC_BC[0, i]] * len(node_tags))
    D_BC_y.extend([DC_BC[1, i]] * len(node_tags))
    D_BC_z.extend([DC_BC[2, i]] * len(node_tags))
    D_coord_x.extend([DC_def[0, i]] * len(node_tags))
    D_coord_y.extend([DC_def[1, i]] * len(node_tags))
    D_coord_z.extend([DC_def[2, i]] * len(node_tags))

# Combine the data into a Dirichlet_BC matrix
Dirichlet_BC = np.column_stack([D_nodes, D_BC_x, D_BC_y, D_BC_z, D_coord_x, D_coord_y, D_coord_z])

# Build the output file path for the boundary conditions
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Input')
os.makedirs(output_dir, exist_ok=True)
Dirichlet_BC_path = os.path.join(output_dir, 'Dirichlet_BC.txt')

# Save the Dirichlet_BC data to a file
np.savetxt(Dirichlet_BC_path, Dirichlet_BC, fmt='%d', delimiter=',')
print(f"Dirichlet BC saved in: {Dirichlet_BC_path}")

# Finalize Gmsh
gmsh.finalize()
