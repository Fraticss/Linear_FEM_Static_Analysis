import os
import time
import numpy as np
import FEM_functions_par as fn_par
import multiprocessing


def main():

    # IMPORT MESH

    # Get the path of the folder where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Build the relative path to the file
    mesh_T_path = os.path.join(script_dir, 'Input', 'mesh_T.txt')
    mesh_P_path = os.path.join(script_dir, 'Input', 'mesh_P.txt')

    try:
        # Specify the delimiter as a comma
        T_matrix = np.loadtxt(mesh_T_path, delimiter=',').T
        P_matrix = np.loadtxt(mesh_P_path, delimiter=',').T
    except Exception as e:
        print(f"Error while loading the files: {e}")
        return

    elapsed_time = np.zeros(6)  # Initialize an array for execution time
    # Start the timer
    start_time = time.time()

    n_elems = np.size(T_matrix, 1)
    n_nodes = np.size(P_matrix, 1)

    # Material properties definitions
    E = 210 * 1e9
    nu = 0.33
    num_cores = multiprocessing.cpu_count()

    # Assembly of Stiffness Matrix K
    try:
        # Pass the number of processors to the function
        K = fn_par.K_global_assembly(P_matrix, T_matrix, E, nu, num_procs=num_cores)
    except Exception as e:
        print(f"Error during the assembly of the stiffness matrix: {e}")
        return
    
    # COMPUTATION TIME
    end_time = time.time()
    elapsed_time[0] = end_time - start_time
    print("Time taken for the assembly of the stiffness matrix: ", elapsed_time[0])


    start_time = time.time()

    # Imposition of Neumann BC
    f = fn_par.f_assembly(P_matrix, T_matrix, num_procs=num_cores)

    # COMPUTATION TIME
    end_time = time.time()
    elapsed_time[1] = end_time - start_time
    print("Time taken for the imposition of Neumann BC: ", elapsed_time[1])

    '''I could try to calculate in F a connectivity matrix node-->matrix 
       so that I can directly sum the contributions of the forces only on the 
       affected nodes. This way I would avoid having to check in the for loop 
       all the elements. In this case, it might not be necessary to parallelize 
       the process or at least modify it to do so in parallel for each forced 
       node and not for all the elements.'''

    start_time = time.time()

    # Application of Dirichlet BC
    K, f, constrained_nodes = fn_par.vincolo_x(K, f, P_matrix)

    # COMPUTATION TIME
    end_time = time.time()
    elapsed_time[2] = end_time - start_time
    print("Time taken for the imposition of Dirichlet BC: ", elapsed_time[2])

    start_time = time.time()

    # Calculation of displacements for unconstrained points
    toll = 1e-6
    u_new = fn_par.lin_syst_solver_par(K, f)
    print('Linear system resolution completed')

    end_time = time.time()
    elapsed_time[3] = end_time - start_time
    print("Time taken for the resolution of the linear system: ", elapsed_time[3])
    
    start_time = time.time()

    # Recompose the u vector by adding the displacements (null) of the constrained points
    u = fn_par.u_composition(u_new, constrained_nodes, n_nodes)

    end_time = time.time()
    elapsed_time[4] = end_time - start_time
    print("Time taken for the recomposition of the u vector: ", elapsed_time[4])

    start_time = time.time()

    # Calculation of element stresses
    sigma, stress_vm, stress_vm_nodes = fn_par.sigma_assembly(P_matrix, T_matrix, u, E, nu, num_procs=num_cores)

    end_time = time.time()
    elapsed_time[5] = end_time - start_time
    print("Time taken for stress calculation: ", elapsed_time[5])

    # Recalculate the position of the nodes after deformation is computed
    P_new = np.zeros((3, n_nodes))
    P_new[0, :] = P_matrix[0, :] + u[0::3]
    P_new[1, :] = P_matrix[1, :] + u[1::3]
    P_new[2, :] = P_matrix[2, :] + u[2::3]

    # COMPUTATION TIME
    print("Total time: ", sum(elapsed_time))

    # Build the path for the Output folder
    output_dir = os.path.join(script_dir, 'Output')
    # Create the Output folder if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Build the full path for the output files
    T_path = os.path.join(output_dir, 'mesh_T.txt')
    P_new_path = os.path.join(output_dir, 'P_new.txt')
    stress_vm_path = os.path.join(output_dir, 'stress_vm.txt')
    stress_vm_nodes_path = os.path.join(output_dir, 'stress_vm_nodes.txt')
    u_nodes_path = os.path.join(output_dir, 'u_nodes.txt')

    # Save the output files in the Output folder
    np.savetxt(T_path, T_matrix, delimiter=',')
    np.savetxt(P_new_path, P_new, delimiter=',')
    np.savetxt(stress_vm_path, stress_vm, delimiter=',')
    np.savetxt(stress_vm_nodes_path, stress_vm_nodes, delimiter=',')
    np.savetxt(u_nodes_path, u, delimiter=',')

    print(f"Files saved in the folder: {output_dir}")


if __name__ == '__main__':
    main()