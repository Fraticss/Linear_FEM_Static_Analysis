import gmsh
import numpy as np
import os

# Build the path to the .msh file in the Input subfolder
input_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Input')
mesh_file = os.path.join(input_dir, 'MESH_FILE_4x.msh')  # Replace with the name of the .msh file

# Load the mesh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output for debugging
gmsh.open(mesh_file)

print(f"Mesh loaded from: {mesh_file}")

# Check the physical groups defined in the file
physical_groups = gmsh.model.getPhysicalGroups()

# Define Neumann BC parameters
face_N = [6]  # Replace with the desired face index
N_BC = np.array([[1],
                 [0], 
                 [0]])  # Boundary condition flags
N_stress = 1e5 * np.array([[1], 
                           [0], 
                           [0]])  # Stress values

# Initialize arrays for Neumann BC
N_nodes = []
N_BC_x = []
N_BC_y = []
N_BC_z = []
N_coord_x = []
N_coord_y = []
N_coord_z = []

# Loop through the specified faces
for i in range(len(face_N)):
    # Find the nodes belonging to the face using Gmsh
    node_tags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(2, face_N[i])

    # Check if nodes were found
    if len(node_tags) == 0:
        print(f"Error: No nodes found for the physical group with index {face_N[i]}.")
        continue

    # Append data for Neumann BC
    N_nodes.extend(node_tags)
    N_BC_x.extend([N_BC[0, i]] * len(node_tags))
    N_BC_y.extend([N_BC[1, i]] * len(node_tags))
    N_BC_z.extend([N_BC[2, i]] * len(node_tags))
    N_coord_x.extend([N_stress[0, i]] * len(node_tags))
    N_coord_y.extend([N_stress[1, i]] * len(node_tags))
    N_coord_z.extend([N_stress[2, i]] * len(node_tags))

# Combine Neumann BC data into a single array
Neumann_BC = np.column_stack([N_nodes, N_BC_x, N_BC_y, N_BC_z, N_coord_x, N_coord_y, N_coord_z])

# Save Neumann BC to a file
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Input')
os.makedirs(output_dir, exist_ok=True)
Neumann_BC_path = os.path.join(output_dir, 'Neumann_BC.txt')

# Save the Neumann_BC data to a file
np.savetxt(Neumann_BC_path, Neumann_BC, fmt='%d', delimiter=',')

print(f"Neumann BC saved to: {Neumann_BC_path}")

# Finalize Gmsh
gmsh.finalize()
