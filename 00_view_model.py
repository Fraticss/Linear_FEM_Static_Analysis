import gmsh
import numpy as np
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# Initialize gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

# Load STEP file
step_file = os.path.join(script_dir, 'Input',"example_model.STEP") 
gmsh.model.occ.importShapes(step_file)
gmsh.model.occ.synchronize()

# Generate the mesh (thetraedrical elements)
gmsh.model.mesh.generate(3)

surfaces = gmsh.model.getEntities(2)  # Obtain surfaces
for  i,surface in enumerate(surfaces):
    gmsh.model.addPhysicalGroup(2, [surface[1]], tag=surface[1])  # Create a physical group for each surface
    gmsh.model.setPhysicalName(2, surface[1], str(surface[1]))  # Use only the surface index as the name

# Set mesh visualization options
gmsh.option.setNumber("Geometry.PointNumbers", 0)
gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
gmsh.option.setNumber("Geometry.LabelType", 2)
gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
gmsh.option.setNumber("Mesh.SurfaceFaces", 0)
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.option.setNumber("Mesh.VolumeFaces", 0)

# Set GUI background to black and text to white
gmsh.option.setNumber("General.ColorScheme", 3)  # Set gradient color to black
gmsh.option.setNumber("General.Light0", 0)  # Disable light

# Start the gui for the mesh visualization
gmsh.fltk.run()

# Close gmsh gui
gmsh.finalize()
