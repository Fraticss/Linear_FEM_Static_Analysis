import gmsh
import numpy as np
import os
import meshio
import pyvista as pv


script_dir = os.path.dirname(os.path.abspath(__file__))

# Initialize gmsh
gmsh.initialize()
#gmsh.option.setNumber("General.Terminal", 1)

# Load STEP file
step_file = os.path.join(script_dir, 'Input',"Cantilivier_beam.STEP") 
gmsh.model.occ.importShapes(step_file)
gmsh.model.occ.synchronize()

# Define mesh element size
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0)  # Minimum element size
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 15)  # Maximum element size

# Mesh quality settings
gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.5)  # Quality threshold
gmsh.option.setNumber("Mesh.Optimize", 1)  # Enable optimization
gmsh.option.setNumber("Mesh.QualityType", 2)  # Quality Criterion : aspect ratio

# Generate the mesh (thetraedrical elements)
gmsh.model.mesh.generate(3)

surfaces = gmsh.model.getEntities(2)  # Obtain surfaces
for  i,surface in enumerate(surfaces):
    gmsh.model.addPhysicalGroup(2, [surface[1]], tag=surface[1])  # Create a physical group for each surface
    gmsh.model.setPhysicalName(2, surface[1], str(surface[1]))  # Use only the surface index as the name


# Save mesh in .msh 4.1

gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.SaveAll", 1)  # Save everything

# Saving mesh as this output file
mesh4_file = os.path.join(script_dir, 'Input', 'MESH_FILE_4x.msh')
gmsh.write(mesh4_file)  # Salva il file in formato MSH 4.1

# Save mesh as .msh 2.2

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 1)  # Salva tutto

# Saving mesh as this output file
mesh22_file = os.path.join(script_dir, 'Input', 'MESH_FILE_22.msh')
gmsh.write(mesh22_file)  # Salva il file in formato MSH 2.2

# Check on element types present in the mesh 
element_types, _, _ = gmsh.model.mesh.getElements()
print("Types of elmenents generated:", element_types)


# MESH VISUALIZATION

# Load .msh 2.2 file
mesh = meshio.read(mesh22_file)

# Extract the connectivity matrix
tetra = mesh.cells_dict.get("tetra")

if tetra is None:
    raise ValueError("The file doesn't contain any thetraedrical element.")

# 'tetra' is the connectivity matrix (n_elem, 4)
T = tetra+1  # Add 1 to go from 0-based a 1-based elements

# Extract the points matrix
P = mesh.points

# Print the shape of T and P
print("Shape matrix T:", T.shape)
print("Shape matrix P:", P.shape)

# Save matricies in the Input Folder as txt files
P_path = os.path.join(script_dir,'Input', 'mesh_P.txt')
T_path = os.path.join(script_dir,'Input', 'mesh_T.txt')

np.savetxt(T_path, T, delimiter=',')
np.savetxt(P_path, P, delimiter=',')

# Set mesh visualization options
gmsh.option.setNumber("Geometry.PointNumbers", 0)
gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
gmsh.option.setNumber("Geometry.LabelType", 2)
gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
gmsh.option.setNumber("Mesh.SurfaceFaces", 0)
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.option.setNumber("Mesh.VolumeFaces", 0)

# Start the gui for the mesh visualization
gmsh.fltk.run()

# Close gmsh gui
gmsh.finalize()

# Mesh Visualization with PyVista

# Extract the points
points = P[:, :]  # Columns X, Y, Z

# Extract the connectivity matrix
connectivity = T[:, :]-1  # Columns Node1, Node2, Node3, Node4 (PyVista uses indicies 0-based)

# Create a connectivity array in PyVista format
cells = np.hstack([np.full((connectivity.shape[0], 1), 4), connectivity]).flatten()

# Define element type (10 = thetraedron in VTK)
cell_type = np.full(connectivity.shape[0], 10)

# Create an object UnstructuredGrid
grid = pv.UnstructuredGrid(cells, cell_type, points)

# Plot the mesh
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color="lightblue", opacity=1.0)
plotter.show()