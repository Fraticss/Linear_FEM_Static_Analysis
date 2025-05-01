import pyvista as pv
import numpy as np
import os

# Get the path of the folder where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
output_folder = 'Output'
# Build the relative path to the files in the 'Output' subfolder
P_new_path = os.path.join(script_dir, output_folder, 'P_new.txt')
T_path = os.path.join(script_dir, output_folder, 'mesh_T.txt')
stress_vm_path = os.path.join(script_dir, output_folder, 'stress_vm.txt')
stress_vm_nodi_path = os.path.join(script_dir, output_folder, 'stress_vm_nodes.txt')
u_nodes_path = os.path.join(script_dir, output_folder, 'u_nodes.txt')
mesh_file_path = os.path.join(script_dir, output_folder, 'MESH_FILE.msh')

# Load the files
P_new = np.loadtxt(P_new_path, delimiter=',')  # Transpose to get the correct matrix
T = np.loadtxt(T_path, delimiter=',').astype(int)
stress_vm = np.loadtxt(stress_vm_path, delimiter=',')
stress_vm_nodi = np.loadtxt(stress_vm_nodi_path, delimiter=',')
u_nodes = np.loadtxt(u_nodes_path, delimiter=',')

# Transpose P_new → (n_nodes, 3)
points = P_new.T

n_elem = T.shape[1]

print("Shape of matrix T:", T.shape)
print("Shape of matrix P:", points.shape)

cells = np.hstack([np.full((n_elem, 1), 4), T.T-1])  # each row: [4, i1, i2, i3, i4]
celltypes = np.full(n_elem, 10)  # VTK code for tetrahedron

# Create the tetrahedral mesh
grid = pv.UnstructuredGrid(cells, celltypes, points)

# Add scalar field to the mesh for the elements
grid["stress_magnitude"] = stress_vm_nodi
'''
grid['disp_x'] = u_nodes[0::3]  # Displacement in x direction
grid['disp_y'] = u_nodes[1::3]  # Displacement in y direction
grid['disp_z'] = u_nodes[2::3]  # Displacement in z direction
'''

# Visualize coloring based on the magnitude
plotter = pv.Plotter()
plotter.add_mesh(
    grid, 
    scalars="stress_magnitude", 
    cmap="jet", 
    show_edges=False, 
    clim=[0, np.max(stress_vm_nodi)]  # Set the limits of the color scale
)

# Add a reference system in the bottom left corner
plotter.add_axes()  

plotter.show()
