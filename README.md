# Linear_FEM_Static_Analysis

This project was developed as a personal initiative to gain a deeper understanding of the Finite Element Method (FEM) and its application in structural analysis. All functions in this repository were written from scratch using only free and widely available Python libraries.

The meshing process is handled via the Gmsh library, which also assists in identifying boundary faces. PyVista is used for mesh and results visualization. Numerical operations are carried out using NumPy and SciPy, while parallel computation is handled through the `multiprocessing` module.

Although the implementation is intentionally simple, it effectively performs static structural analysis and yields results consistent with analytical solutions as the mesh is refined.

If you find any errors or have suggestions for improvements, feel free to reach out. I paused further development to focus on upcoming exams, though there are still several parts I would like to optimize and extend with more advanced features.

A standalone GUI-based application is planned for future release on GitHub. It will allow anyone to perform static analysis of structures made of isotropic materials—without needing to modify the source code.

---

## 🧭 Workflow Overview

This repository provides a set of Python scripts to perform linear static FEM analysis on 3D models using tetrahedral elements and linear interpolation. The workflow is divided into six numbered scripts (`00` to `05`), each corresponding to a specific step in the process.

---

## ⚙️ Prerequisites

Ensure the following Python libraries are installed:

- `gmsh`
- `numpy`
- `os`
- `meshio`
- `pyvista`
- `time`
- `multiprocessing`
- `scipy`

---

## 🔧 Workflow Steps

### Step 1: Prepare the Input Directory

- Create a folder named `Input` in the root directory.
- Place the `.stp` or `.STEP` file of your 3D model inside this folder.

---

### Script 00: Model Viewer (`00_view_model.py`)

Visualizes the 3D model and displays the face indices using Gmsh—essential for defining boundary conditions in later steps.

- **Edit line 12** to specify your STEP file:
  ```python
  step_file = "your_model.step"
  ```

---

### Script 01: Mesh Generation (`01_gmsh_mesh.py`)

Generates a tetrahedral mesh from the 3D model and visualizes it using both Gmsh and PyVista.

- **Edit the following lines:**
  ```python
  step_file = "your_model.step"             # line 15
  Mesh.CharacteristicLengthMin = 0          # line 20
  Mesh.CharacteristicLengthMax = 5          # line 21
  ```
  > Smaller `CharacteristicLengthMax` results in a finer mesh.

- **Run the script** to:
  - Generate and visualize the mesh.
  - Save all relevant mesh files in the `Input` folder.

---

### Script 02: Apply Dirichlet Boundary Conditions (`02_Dirichlet_BC.py`)

Defines essential (Dirichlet) boundary conditions on selected mesh faces.

- **Edit lines 20–32:**
  ```python
  face_DC = [6, 3, 7]

  DC_BC = np.array([
      [1, 0, 1],
      [1, 1, 0],
      [1, 0, 0]
  ])

  DC_def = np.array([
      [0, 0, 20],
      [0, 100, 0],
      [0, 0, 0]
  ])
  ```

- **Run the script** to save the Dirichlet conditions in the `Input` folder.

---

### Script 03: Apply Neumann Boundary Conditions (`03_Neumann_BC.py`)

Applies pressure loads (Neumann conditions) to specified faces.

- **Edit lines 20–32:**
  ```python
  face_N = [1, 2]

  N_BC = np.array([
      [0, 0],
      [0, 1],
      [1, 1]
  ])

  N_stress = np.array([
      [0, 0],
      [0, 0],
      [50, 100]
  ])
  ```

> This script works similarly to `02_Dirichlet_BC.py`.

- **Run the script** to save Neumann conditions in the `Input` folder.

---

### Script 04: FEM Solver (`04_solver.py`)

Solves the FEM system based on the mesh and boundary conditions.

- **Edit lines 35–37:**
  ```python
  E = 210e9          # Young’s modulus
  nu = 0.3           # Poisson’s ratio
  num_cores = 4      # Number of CPU cores for parallel computation
  ```

- **Run the script** to:
  - Compute nodal displacements, stress, and more.
  - Save the results in the `Output` folder.

---

### Script 05: Post-Processing (`05_postprocess.py`)

Visualizes simulation results using PyVista.

- **Edit line 38** to select the desired output (use only one of the following lines):
  ```python
  grid["stress_magnitude"] = stress_vm_nodi  # Nodal Von Mises stress
  grid["stress_element"] = stress_vm         # Elemental Von Mises stress
  grid['disp_x'] = u_nodes[0::3]             # X displacement
  grid['disp_y'] = u_nodes[1::3]             # Y displacement
  grid['disp_z'] = u_nodes[2::3]             # Z displacement
  ```

- **Run the script** to visualize the deformed structure with the selected result type.

---

## 💡 Notes

- Ensure all necessary files are present in the `Input` folder before executing each script.
- Scripts are intended to be run in order (`00` to `05`) for consistency.
- If errors occur, check your inputs and confirm all required libraries are installed.

---
## License
This project is licensed under the MIT License – see the [LICENSE](./LICENSE) file for details.


Enjoy exploring the world of finite element analysis!
