# Linear_FEM_Static_Analysis

This is the result of a project carried out by me to have a better understanding of how FEA works and how it is employed in the structural analysis. For this reason I personally wrote every function present in this repository by using only free python libraries easy downloadable by anyone. I employed gmsh library to create the mesh and to understand where to impose buondary conditions, pyvista to visualize results and the mesh itself, and other libraries like numpy, scipy and multiprocessing to create and solve the linear system. In conclusion, this was a really helpful project and lead me to create a very simple but still effective and free way to perform static analysis that I even compared with analytical results (likely the method converges with analytical results by agumenting the number of elements employed!).

# WORKFLOW
This repository contains a series of Python scripts for performing finite element analysis (FEA) on a 3D model by using thetraedral elements and linear interpolating functions. The scripts are numbered from `01` to `05`, and each script corresponds to a specific step in the workflow.

---

## Prerequisites

1. Install the required Python libraries:
   - `gmsh`
   - `numpy`
   - `os`
   - `meshio`
   - `pyvista`
   - `time`
   - `multiprocessing`
   - `scipy`
   
---

## Workflow

### Step 1: Prepare the Input Folder
- Create a folder named `Input` in the same directory as the scripts.
- Place the `.stp` or `.STEP` file of the 3D model you want to analyze into the `Input` folder.

---
### Script 00: Model Visualization (`00_view_model.py`)
This script plot the 3D model provided in the `Input` folder by using gmsh to see the index corresponding to each face (usefull for scripts 02 and 03).

1. **Modify the script:**
  - Set the name of the `.STEP` file in the `step_file` variable (line 12).

### Script 01: Mesh Generation (`01_gmsh_mesh.py`)
This script generates a tetrahedral mesh from the 3D model provided in the `Input` folder.

1. **Modify the script:**
   - Set the name of the `.STEP` file in the `step_file` variable (line 15).
   - Adjust the minimum and maximum characteristic lengths of the tetrahedral elements:
     - `Mesh.CharacteristicLengthMin` (line 20) (usually set to zero)
     - `Mesh.CharacteristicLengthMax` (line 21)
     NOTE: Given the dimensions of the 3d model, a lower h_max corresponds to more elements in the mesh

2. **Run the script:**
   - The script will generate a tetrahedral mesh from the 3D model.
   - The mesh will be displayed in the Gmsh GUI and visualized using PyVista.
   - All necessary files for subsequent steps will be saved in the `Input` folder.

---

### Script 02: Apply Dirichlet Boundary Conditions (`02_Dirichlet_BC.py`)
This script allows you to impose essential (Dirichlet) boundary conditions on specific faces of the mesh.

1. **Modify the script:**
   - Edit lines 20–32:
     - `face_DC`: A row vector containing the indices of the faces to constrain.
         example: face_DC = [6,3,7]
     - `DC_BC`: A matrix where:
       - Columns correspond to the indices in `face_DC`.
       - Rows correspond to the x, y, and z directions. Use `1` to constrain a direction and `0` to leave it free.
         example: DC_BC = np.array([
                                    [1,0,1],
                                    [1,1,0],
                                    [1,0,0]
                                ])
         In this case the face 6 is constrained on x,y,z directions, face 3 only on y and face 7 on x.
     - `DC_def`: A matrix structured like `DC_BC`, but instead of `0` or `1`, specify the imposed displacement values.
         example: DC_def = np.array([
                                    [0,0,20],
                                    [0,100,0],
                                    [0,0,0]
                                ])
         In this case, the first face will stay still since its movement is imposed to be zero on every direction.
         The second column corresponds to the face 3 that is constrained only along y, so only the value of the second row of the            third column will be considered as a constrain. For the third column will be considered only the value in the first row.

2. **Run the script:**
   - The script will save the Dirichlet boundary conditions in the `Input` folder.

---

### Script 03: Apply Neumann Boundary Conditions (`03_Neumann_BC.py`)
This script imposes constant pressure (Neumann boundary conditions) on specific faces of the mesh.

1. **Modify the script:**
   - Edit lines 20–32:
     - `face_N`: A row vector containing the indices of the faces to apply pressure.
     - `N_BC`: A matrix where:
       - Columns correspond to the indices in `face_N`.
       - Rows correspond to the x, y, and z directions. Use `1` to apply pressure in a direction and `0` to leave it free.
     - `N_stress`: A matrix structured like `N_BC`, but instead of `0` or `1`, specify the stress values.
     NOTE: this script works exactly as 02_Dirichlet_BC.py so check its examples.

2. **Run the script:**
   - The script will save the Neumann boundary conditions in the `Input` folder.

---

### Script 04: Solver (`04_solver.py`)
This script solves the finite element problem using the mesh and boundary conditions defined in the previous steps.

1. **Modify the script:**
   - Set the material properties:
     - Young's Modulus (`E`) (line 35).
     - Poisson's Ratio (`nu`) (line 36).
   - Number of cores to use (`num_cores`) (line 37).
     - Optionally, set the number of CPU cores to use for parallel processing (default is the maximum available).

2. **Run the script:**
   - The script will compute the displacements, stresses, and other results.
   - The results will be saved in the `Output` folder.

---

### Script 05: Post-Processing (`05_postprocess.py`)
This script reads the results from the `Output` folder and visualizes the deformed model with Von Mises stresses or displacements along the x, y, or z directions.

1. **Modify the script:**
   - Modify line 38 by coping one of the lines in the following section:
        grid["stress_magnitude"] = stress_vm_nodi  # Stress magnitude at nodes
        grid["stress_element"] = stress_vm  # Stress magnitude at elements
        grid['disp_x'] = u_nodes[0::3]  # Displacement in x direction
        grid['disp_y'] = u_nodes[1::3]  # Displacement in y direction
        grid['disp_z'] = u_nodes[2::3]  # Displacement in z direction

2. **Run the script:**
   - The script will plot the deformed model and display the results using PyVista.

---

## Notes
- Ensure that the `Input` folder contains all necessary files before running each script.
- The scripts are designed to work sequentially, so follow the order from `01` to `05`.
- If you encounter any errors, check the input files and ensure that the required libraries are installed.

Enjoy your finite element analysis!
