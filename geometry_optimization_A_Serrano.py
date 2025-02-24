import os
import numpy as np
import math


def parse_molecule_file(file_path):
    atoms = []
    coordinates = []
    bonds = []

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

            for line in lines[1:]:
                parts = line.split()

                # Parse atomic data
                if len(parts) >= 10:  # Assuming atoms have at least 10 columns
                    try:
                        x, y, z = map(float, parts[:3])
                        element = parts[3]
                        atoms.append({"element": element, "coordinates": (x, y, z)})
                        coordinates.append([x, y, z])
                    except ValueError:
                        print("Error in the input structure, you should check it and then run again")
                        break

                # Parse bond data
                elif len(parts) >= 2:  # Assuming bonds have at least two integers
                    try:
                        bond = list(map(int, parts[:2]))
                        if len(bond) == 2:  # Ensure exactly two indices
                            bonds.append(bond)
                    except ValueError:
                        print("Error in the input structure, you should check it and then run again")
                        break

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return [], np.array([]), np.array([])
    except Exception as e:
        print(f"An error occurred: {e}")
        return [], np.array([]), np.array([])

    # Convert lists to numpy arrays
    coordinates_array = np.array(coordinates, dtype=float) if coordinates else np.array([])
    bonds_array = np.array(bonds, dtype=int) if bonds else np.array([])

    return atoms, coordinates_array, bonds_array


def calculate_bond_distances(coordinates, bonds):

    distances = []
    for bond in bonds:
        try:
            atom1_idx, atom2_idx = bond
            x1, y1, z1 = coordinates[atom1_idx - 1]
            x2, y2, z2 = coordinates[atom2_idx - 1]
            distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            distances.append((atom1_idx, atom2_idx, distance))
        except IndexError:
            print(f"Invalid bond indices: {bond}")
            break

    distances_array = np.array(distances, dtype=float)

    return distances_array



def calculate_bond_energies(distances, atoms, bonds):

    bond_energies = []
    for bond in bonds:
        try:
            atom1_idx, atom2_idx = bond
            atom1 = atoms[atom1_idx - 1]
            atom2 = atoms[atom2_idx - 1]

            element1 = atom1["element"]
            element2 = atom2["element"]

            # Find the corresponding distance for this bond
            rb = next(
                distance[2] for distance in distances if distance[0] == atom1_idx and distance[1] == atom2_idx
            )


            if element1 == 'C' and element2 == 'C':
                energy = 300 * (rb - 1.53) ** 2
                bond_energies.append((atom1_idx, atom2_idx, energy))
            elif element1 == 'C' and element2 == 'H':
                energy = 350 * (rb - 1.11) ** 2
                bond_energies.append((atom1_idx, atom2_idx, energy))


        except StopIteration:
            print(f"Distance not found for bond: {bond}")
            break

    bond_energies_array = np.array(bond_energies, dtype=float)

    return bond_energies_array



def calculate_angles(coordinates, bonds):
    angles = []
    seen_angles = set()

    for bond1 in bonds:
        for bond2 in bonds:
            # Ensure the bonds share a common atom
            common_atom = set(bond1) & set(bond2)
            if len(common_atom) == 1:
                common_atom_idx = list(common_atom)[0]
                other_atom1 = (set(bond1) - common_atom).pop()
                other_atom2 = (set(bond2) - common_atom).pop()

                # Ensure unique angles by sorting atoms in a consistent order
                angle_key = tuple(sorted([other_atom1, common_atom_idx, other_atom2]))
                if angle_key in seen_angles:
                    continue
                seen_angles.add(angle_key)

                try:
                    # Get coordinates of the three atoms
                    c1 = coordinates[common_atom_idx - 1]
                    c2 = coordinates[other_atom1 - 1]
                    c3 = coordinates[other_atom2 - 1]

                    # Calculate vectors
                    v1 = np.array(c2) - np.array(c1)
                    v2 = np.array(c3) - np.array(c1)

                    # Calculate angle
                    dot_product = np.dot(v1, v2)
                    norm_v1 = np.linalg.norm(v1)
                    norm_v2 = np.linalg.norm(v2)
                    angle = math.acos(dot_product / (norm_v1 * norm_v2))

                    angles.append((other_atom1, common_atom_idx, other_atom2, angle))
                except IndexError:
                    print(f"Invalid bond indices for angle calculation: {bond1}, {bond2}")
                    continue

    angles_array = np.array(angles, dtype=float)
    return angles_array



def calculate_angle_energies(angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other):

    angle_energies = []

    for angle in angles:
        atom1_idx, common_atom_idx, atom2_idx, angle_value = angle
        atom1 = atoms[int(atom1_idx) - 1]
        atom2 = atoms[int(atom2_idx) - 1]
        common_atom = atoms[int(common_atom_idx) - 1]

        element1 = atom1["element"]
        element2 = atom2["element"]
        common_element = common_atom["element"]

        # Check if the angle involves C-C-C
        if element1 == 'C' and common_element == 'C' and element2 == 'C':
            k = k_ccc
            angle0 = angle0_ccc
        else:
            k = k_other
            angle0 = angle0_other

        # Calculate energy
        energy = k * (angle_value - angle0) ** 2
        angle_energies.append((atom1_idx, common_atom_idx, atom2_idx, energy))

    angle_energies_array = np.array(angle_energies, dtype=float)
    return angle_energies_array


def calculate_dihedral_angles(coordinates, bonds):
    dihedrals = []
    seen_dihedrals = set()

    for bond1 in bonds:
        for bond2 in bonds:
            for bond3 in bonds:
                # Ensure bonds are sequential
                common1 = set(bond1) & set(bond2)
                common2 = set(bond2) & set(bond3)

                if len(common1) == 1 and len(common2) == 1:
                    common1_atom = list(common1)[0]
                    common2_atom = list(common2)[0]

                    if common1_atom == common2_atom:
                        continue

                    atom1 = (set(bond1) - common1).pop()
                    atom4 = (set(bond3) - common2).pop()

                    # Ensure unique dihedrals
                    dihedral_key = tuple(sorted([atom1, common1_atom, common2_atom, atom4]))
                    if dihedral_key in seen_dihedrals:
                        continue
                    seen_dihedrals.add(dihedral_key)

                    try:
                        # Get coordinates of four atoms
                        p1 = np.array(coordinates[atom1 - 1])
                        p2 = np.array(coordinates[common1_atom - 1])
                        p3 = np.array(coordinates[common2_atom - 1])
                        p4 = np.array(coordinates[atom4 - 1])

                        # Calculate vectors and normal vectors
                        b1 = p2 - p1
                        b2 = p3 - p2
                        b3 = p4 - p3

                        n1 = np.cross(b1, b2)
                        n2 = np.cross(b2, b3)
                        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

                        # Calculate dihedral angle
                        x = np.dot(n1, n2)
                        y = np.dot(m1, n2)
                        angle = (-1) * math.atan2(y, x)

                        dihedrals.append((atom1, common1_atom, common2_atom, atom4, angle))
                    except IndexError:
                        print(f"Invalid bond indices for dihedral calculation: {bond1}, {bond2}, {bond3}")
                        continue

    dihedrals_array = np.array(dihedrals, dtype=float)
    return dihedrals_array


def calculate_dihedral_energies(dihedrals, a, n):

    dihedral_energies = []

    for dihedral in dihedrals:
        try:
            atom1, atom2, atom3, atom4, angle = dihedral
            # Calculate energy based on the given formula
            energy = a * (1 + math.cos(n * angle))
            dihedral_energies.append((atom1, atom2, atom3, atom4, energy))
        except Exception as e:
            print(f"Error calculating energy for dihedral {dihedral}: {e}")

    dihedral_energies_array = np.array(dihedral_energies, dtype=float)
    return dihedral_energies_array




def calculate_all_unique_distances(coordinates, bonds):
    bond_set = set(tuple(sorted(bond)) for bond in bonds)  # Direct bonds

    # Determine bonded pairs via a common atom (A-C-B)
    indirect_bonds = set()
    for bond1 in bonds:
        for bond2 in bonds:
            common_atom = set(bond1) & set(bond2)
            if len(common_atom) == 1:
                common_atom = list(common_atom)[0]
                atom1 = (set(bond1) - {common_atom}).pop()
                atom2 = (set(bond2) - {common_atom}).pop()
                indirect_bonds.add(tuple(sorted([atom1, atom2])))

    excluded_pairs = bond_set | indirect_bonds

    # Calculate distances for all unique atom pairs
    n_atoms = len(coordinates)
    unique_distances = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if (i + 1, j + 1) not in excluded_pairs:  # Exclude bonded pairs
                x1, y1, z1 = coordinates[i]
                x2, y2, z2 = coordinates[j]
                distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
                unique_distances.append((i + 1, j + 1, distance))

    unique_distances_array = np.array(unique_distances, dtype=float)
    return unique_distances_array


def calculate_lennard_jones_energy(unique_distances, atoms, parameters):
    lj_energies = []

    for pair in unique_distances:
        atom1_idx, atom2_idx, r = pair
        atom1 = atoms[int(atom1_idx) - 1]
        atom2 = atoms[int(atom2_idx) - 1]

        element1 = atom1["element"]
        element2 = atom2["element"]

        # Determine the A and B parameters based on atom pair type
        pair_type = tuple(sorted([element1, element2]))
        if pair_type in parameters:
            lena_a, lena_b = parameters[pair_type]
            # Calculate Lennard-Jones potential
            energy = (lena_a / r**12) - (lena_b / r**6)
            lj_energies.append((atom1_idx, atom2_idx, energy))
        else:
            print(f"No parameters defined for pair {pair_type}")

    lj_energies_array = np.array(lj_energies, dtype=float)
    return lj_energies_array



def calculate_b_matrix(coordinates, bonds, angles, dihedrals):
    n_atoms = len(coordinates)
    n_bonds = len(bonds)
    n_angles = len(angles)
    n_dihedrals = len(dihedrals)

    # Total rows in the B matrix: bonds + angles + dihedrals
    total_rows = n_bonds + n_angles + n_dihedrals
    b_matrix = np.zeros((total_rows, n_atoms * 3))

    # Stretching contributions (bond derivatives)
    for bond_idx, bond in enumerate(bonds):
        atom1_idx, atom2_idx = bond

        atom1_coords = coordinates[atom1_idx - 1]
        atom2_coords = coordinates[atom2_idx - 1]

        bond_vector = atom2_coords - atom1_coords
        bond_distance = np.linalg.norm(bond_vector)

        if bond_distance > 0:
            v2 = bond_vector / bond_distance
        else:
            print(f"Warning: Bond distance is zero for bond {bond}")
            v2 = np.zeros(3)

        for i in range(3):
            b_matrix[bond_idx, (atom1_idx - 1) * 3 + i] = -v2[i]
            b_matrix[bond_idx, (atom2_idx - 1) * 3 + i] = v2[i]

    # Bending contributions (angle derivatives)
    for angle_idx, angle in enumerate(angles):
        row_idx = n_bonds + angle_idx
        atom1_idx, common_atom_idx, atom2_idx = map(int, angle[:3])


        atom1_coords = coordinates[atom1_idx - 1]
        common_coords = coordinates[common_atom_idx - 1]
        atom2_coords = coordinates[atom2_idx - 1]

        rba = common_coords - atom1_coords
        rbc = common_coords - atom2_coords
        p = np.cross(rba, rbc)
        norm_rba = np.linalg.norm(rba)
        norm_rbc = np.linalg.norm(rbc)
        norm_p = np.linalg.norm(p)

        if norm_rba == 0 or norm_rbc == 0 or norm_p == 0:
            print(f"Zero vector encountered for atoms {atom1_idx}, {common_atom_idx}, {atom2_idx}")
            continue

        grad_angle_common = ((-np.cross(rba, p) / (norm_rba**2 * norm_p)) +
                             (np.cross(rbc, p) / (norm_rbc**2 * norm_p)))
        grad_angle_atom1 = np.cross(rba, p) / (norm_rba**2 * norm_p)
        grad_angle_atom2 = -np.cross(rbc, p) / (norm_rbc**2 * norm_p)

        for i in range(3):
            b_matrix[row_idx, (common_atom_idx - 1) * 3 + i] = -grad_angle_common[i]
            b_matrix[row_idx, (atom1_idx - 1) * 3 + i] = -grad_angle_atom1[i]
            b_matrix[row_idx, (atom2_idx - 1) * 3 + i] = -grad_angle_atom2[i]

    # Torsional contributions (dihedral derivatives)
    for dihedral_idx, dihedral in enumerate(dihedrals):
        row_idx = n_bonds + n_angles + dihedral_idx
        atom_a, atom_b, atom_c, atom_d = map(int, dihedral[:4])


        r_a = coordinates[atom_a - 1]
        r_b = coordinates[atom_b - 1]
        r_c = coordinates[atom_c - 1]
        r_d = coordinates[atom_d - 1]

        r_ab = r_b - r_a
        r_bc = r_c - r_b
        r_cd = r_d - r_c
        r_ac = r_c - r_a
        r_bd = r_d - r_b

        t = np.cross(r_ab, r_bc)
        u = np.cross(r_bc, r_cd)

        t_norm = np.linalg.norm(t)
        u_norm = np.linalg.norm(u)
        r_bc_norm = np.linalg.norm(r_bc)

        if t_norm == 0 or u_norm == 0 or r_bc_norm == 0:
            print(f"Warning: Zero norm encountered in dihedral {dihedral}")
            continue



        term_a = np.cross(np.cross(t, r_bc) / ((t_norm ** 2) * r_bc_norm), r_bc)
        term_d = np.cross(np.cross(-u, r_bc) / ((u_norm ** 2) * r_bc_norm), r_bc)
        term_b = (np.cross(r_ac, (np.cross(t, r_bc)/((t_norm ** 2) * r_bc_norm)))) + (np.cross((np.cross(-u, r_bc)/((u_norm ** 2) * r_bc_norm)), r_cd))
        term_c = (np.cross((np.cross(t, r_bc)/((t_norm ** 2) * r_bc_norm)), r_ab)) + (np.cross(r_bd, (np.cross(-u, r_bc)/((u_norm ** 2) * r_bc_norm))))

        for i in range(3):
            b_matrix[row_idx, (atom_a - 1) * 3 + i] = term_a[i]
            b_matrix[row_idx, (atom_d - 1) * 3 + i] = term_d[i]
            b_matrix[row_idx, (atom_b - 1) * 3 + i] = term_b[i]
            b_matrix[row_idx, (atom_c - 1) * 3 + i] = term_c[i]

    return b_matrix


def calculate_stretching_gradient(coordinates, bonds, distances, k_values, r0_values):
    gradient = np.zeros_like(coordinates)

    for i, bond in enumerate(bonds):
        atom1_idx, atom2_idx = bond
        try:
            coord1 = coordinates[atom1_idx - 1]
            coord2 = coordinates[atom2_idx - 1]

            bond_vector = np.array(coord2) - np.array(coord1)
            bond_length = distances[i, 2]

            if bond_length == 0:
                raise ValueError(f"Zero-length bond found between atoms {atom1_idx} and {atom2_idx}.")

            v2 = bond_vector / bond_length
            r0 = r0_values[i]
            k = k_values[i]

            delta_r = bond_length - r0
            force = 2 * k * delta_r * v2

            gradient[atom1_idx - 1] -= force  # Negative contribution for atom1
            gradient[atom2_idx - 1] += force  # Positive contribution for atom2
        except IndexError:
            print(f"Invalid bond indices: {bond}")
        except ValueError as e:
            print(e)

    return gradient



def calculate_bending_gradient(coordinates, atoms, angles, k_ccc, angle0_ccc, k_other, angle0_other):

    bending_gradient = np.zeros_like(coordinates)

    for angle in angles:
        atom1_idx, common_atom_idx, atom2_idx = map(int, angle[:3])
        angle_value = angle[3]
        atom1_coords = coordinates[atom1_idx - 1]
        common_coords = coordinates[common_atom_idx - 1]
        atom2_coords = coordinates[atom2_idx - 1]

        # Determine constants for the angle type
        common_element = atoms[common_atom_idx - 1]["element"]
        element1 = atoms[atom1_idx - 1]["element"]
        element2 = atoms[atom2_idx - 1]["element"]

        if element1 == 'C' and common_element == 'C' and element2 == 'C':
            k = k_ccc
            angle0 = angle0_ccc
        else:
            k = k_other
            angle0 = angle0_other

        # Calculate vectors and distances
        rba = common_coords - atom1_coords
        rbc = common_coords - atom2_coords
        p = np.cross(rba, rbc)
        norm_rba = np.linalg.norm(rba)
        norm_rbc = np.linalg.norm(rbc)
        norm_p = np.linalg.norm(p)

        if norm_rba == 0 or norm_rbc == 0 or norm_p == 0:
            print(f"Zero vector encountered for atoms {atom1_idx}, {common_atom_idx}, {atom2_idx}")
            continue

        # Calculate gradient components for each atom
        grad_angle_common = ((-np.cross(rba, p) / (norm_rba**2 * norm_p)) +
                             (np.cross(rbc, p) / (norm_rbc**2 * norm_p)))
        grad_angle_atom1 = np.cross(rba, p) / (norm_rba**2 * norm_p)
        grad_angle_atom2 = -np.cross(rbc, p) / (norm_rbc**2 * norm_p)

        # Energy gradient contribution
        delta_angle = angle_value - angle0
        prefactor = 2 * k * delta_angle

        bending_gradient[common_atom_idx - 1] -= prefactor * grad_angle_common
        bending_gradient[atom1_idx - 1] -= prefactor * grad_angle_atom1
        bending_gradient[atom2_idx - 1] -= prefactor * grad_angle_atom2

    return bending_gradient



def calculate_torsional_energy_gradient(coordinates, dihedrals, a, n):

    dihedral_gradient = np.zeros_like(coordinates)

    for dihedral in dihedrals:
        atom_a, atom_b, atom_c, atom_d = map(int, dihedral[:4])
        angle = dihedral[4]

        # Get coordinates
        r_a = coordinates[atom_a - 1]
        r_b = coordinates[atom_b - 1]
        r_c = coordinates[atom_c - 1]
        r_d = coordinates[atom_d - 1]

        # Calculate vectors
        r_ab = r_b - r_a
        r_bc = r_c - r_b
        r_cd = r_d - r_c
        r_ac = r_c - r_a
        r_bd = r_d - r_b

        # Calculate cross products
        t = np.cross(r_ab, r_bc)
        u = np.cross(r_bc, r_cd)

        # Calculate norms and distances
        t_norm = np.linalg.norm(t)
        u_norm = np.linalg.norm(u)
        r_bc_norm = np.linalg.norm(r_bc)

        if t_norm == 0 or u_norm == 0 or r_bc_norm == 0:
            print(f"Warning: Zero norm encountered in dihedral {dihedral}")
            continue

        # Calculate terms
        term_a = np.cross(np.cross(t, r_bc) / ((t_norm ** 2) * r_bc_norm), r_bc)
        term_d = np.cross(np.cross(-u, r_bc) / ((u_norm ** 2) * r_bc_norm), r_bc)

        term_b = (np.cross(r_ac, (np.cross(t, r_bc)/((t_norm ** 2) * r_bc_norm)))) + (np.cross((np.cross(-u, r_bc)/((u_norm ** 2) * r_bc_norm)), r_cd))

        term_c = (np.cross((np.cross(t, r_bc)/((t_norm ** 2) * r_bc_norm)), r_ab)) + (np.cross(r_bd, (np.cross(-u, r_bc)/((u_norm ** 2) * r_bc_norm))))



        # Calculate gradient contribution
        energy_prefactor = -n * a * math.sin(n * angle)
        dihedral_gradient[atom_a - 1] += energy_prefactor * term_a
        dihedral_gradient[atom_d - 1] += energy_prefactor * term_d
        dihedral_gradient[atom_b - 1] += energy_prefactor * term_b
        dihedral_gradient[atom_c - 1] += energy_prefactor * term_c

    return dihedral_gradient


def calculate_vdw_energy_gradient(coordinates, unique_distances, atoms, parameters):

    gradient = np.zeros_like(coordinates)

    for pair in unique_distances:
        atom1_idx, atom2_idx, r = pair
        atom1 = atoms[int(atom1_idx) - 1]
        atom2 = atoms[int(atom2_idx) - 1]

        element1 = atom1["element"]
        element2 = atom2["element"]

        # Determine the A and B parameters based on the atom pair type
        pair_type = tuple(sorted([element1, element2]))
        if pair_type in parameters:
            lena_a, lena_b = parameters[pair_type]

            # Calculate the gradient contribution
            if r > 0:  # Avoid division by zero
                factor = (-12 * lena_a / r ** 14) + (6 * lena_b / r ** 8)
                delta_r = coordinates[int(atom1_idx) - 1] - coordinates[int(atom2_idx) - 1]
                grad_contribution = delta_r * factor

                gradient[int(atom1_idx) - 1] += grad_contribution
                gradient[int(atom2_idx) - 1] -= grad_contribution
        else:
            print(f"No parameters defined for pair {pair_type}")

    return gradient


def calculate_total_potential_energy(bond_energies, angle_energies, dihedral_energies, lj_energies):
    # Sum the bond stretching energies
    total_stretching_energy = np.sum(bond_energies[:, 2]) if bond_energies.size else 0

    # Sum the angle bending energies
    total_bond_energy = np.sum(angle_energies[:, 3]) if angle_energies.size else 0

    # Sum the dihedral torsional energies
    total_dihedral_energy = np.sum(dihedral_energies[:, 4]) if dihedral_energies.size else 0

    # Sum the Lennard-Jones energies
    total_lj_energy = np.sum(lj_energies[:, 2]) if lj_energies.size else 0

    # Calculate total potential energy
    total_potential_energy = total_stretching_energy + total_bond_energy + total_dihedral_energy + total_lj_energy

    return total_potential_energy, total_stretching_energy, total_bond_energy, total_dihedral_energy, total_lj_energy


def calculate_total_gradient(coordinates, bonds, distances, k_values, r0_values, angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other, dihedrals, a, n, unique_distances, parameters):

    # Calculate individual gradients
    stretching_gradient = calculate_stretching_gradient(coordinates, bonds, distances, k_values, r0_values)
    bending_gradient = calculate_bending_gradient(coordinates, atoms, angles, k_ccc, angle0_ccc, k_other, angle0_other)
    torsional_gradient = calculate_torsional_energy_gradient(coordinates, dihedrals, a, n)
    vdw_gradient = calculate_vdw_energy_gradient(coordinates, unique_distances, atoms, parameters)

    # Sum all gradients to get the total gradient
    total_gradient = stretching_gradient + bending_gradient + torsional_gradient + vdw_gradient

    return total_gradient



def initialize_inverse_hessian(num_atoms):
    dim = num_atoms * 3  # 3 dimensions (x, y, z) per atom
    m = np.zeros((dim, dim))
    np.fill_diagonal(m, 1 / 300)

    np.set_printoptions(threshold=np.inf)

    return m


def bfgs_geometry_optimization(coordinates, atoms, bonds, k_ccc, angle0_ccc, k_other, angle0_other, a, n, parameters, k_values, r0_values, max_iterations=100, rms_gradient_threshold=0.001):
    num_atoms = len(coordinates)
    m = initialize_inverse_hessian(num_atoms)

    for iteration in range(max_iterations):
        # Step 1: Compute total potential energy and gradient
        distances = calculate_bond_distances(coordinates, bonds)
        bond_energies = calculate_bond_energies(distances, atoms, bonds)
        angles = calculate_angles(coordinates, bonds)
        angle_energies = calculate_angle_energies(angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other)
        dihedrals = calculate_dihedral_angles(coordinates, bonds)
        dihedral_energies = calculate_dihedral_energies(dihedrals, a, n)
        unique_distances = calculate_all_unique_distances(coordinates, bonds)
        lj_energies = calculate_lennard_jones_energy(unique_distances, atoms, parameters)

        total_gradient = calculate_total_gradient(coordinates, bonds, distances, k_values, r0_values, angles, atoms,
                                                  k_ccc, angle0_ccc, k_other, angle0_other, dihedrals, a, n,
                                                  unique_distances, parameters)
        rms_gradient = np.sqrt(np.mean(total_gradient ** 2))

        print(f"Iteration {iteration}: RMS Gradient = {rms_gradient:.6f}")

        # Check convergence
        if rms_gradient < rms_gradient_threshold:
            print("Converged: RMS gradient below threshold.")
            break

        # Step 2: Set search direction p(k)
        p = -np.dot(m, total_gradient.flatten()).reshape(coordinates.shape)

        # Step 3: Line search with Wolfe condition
        alpha = 0.8
        while True:
            new_coordinates = coordinates + alpha * p
            new_distances = calculate_bond_distances(new_coordinates, bonds)
            new_bond_energies = calculate_bond_energies(new_distances, atoms, bonds)
            new_angles = calculate_angles(new_coordinates, bonds)
            new_angle_energies = calculate_angle_energies(new_angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other)
            new_dihedrals = calculate_dihedral_angles(new_coordinates, bonds)
            new_dihedral_energies = calculate_dihedral_energies(new_dihedrals, a, n)
            new_unique_distances = calculate_all_unique_distances(new_coordinates, bonds)
            new_lj_energies = calculate_lennard_jones_energy(new_unique_distances, atoms, parameters)

            new_total_potential_energy, *_ = calculate_total_potential_energy(new_bond_energies, new_angle_energies,
                                                                              new_dihedral_energies, new_lj_energies)
            old_total_potential_energy, *_ = calculate_total_potential_energy(bond_energies, angle_energies,
                                                                              dihedral_energies, lj_energies)

            wolfe_lhs = new_total_potential_energy
            wolfe_rhs = old_total_potential_energy + 0.1 * alpha * np.dot(p.flatten(), total_gradient.flatten())

            if wolfe_lhs <= wolfe_rhs:
                break

            alpha *= 0.8

        # Step 4: Update coordinates
        s = alpha * p
        coordinates += s

        # Step 5: Update inverse Hessian
        new_total_gradient = calculate_total_gradient(new_coordinates, bonds, new_distances, k_values, r0_values,
                                                      new_angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other,
                                                      new_dihedrals, a, n, new_unique_distances, parameters)
        y = new_total_gradient - total_gradient
        v = np.dot(m, y.flatten()).reshape(y.shape)

        sy = np.dot(s.flatten(), y.flatten())
        if sy == 0:
            print("Warning: Division by zero in Hessian update.")
            continue

        m += ((sy + np.dot(y.flatten(), v.flatten())) * np.outer(s, s) / (sy ** 2)
              - (np.outer(v, s) + np.outer(s, v)) / sy)

    else:
        print("Warning: Maximum iterations reached without convergence.")

    return coordinates


def calculate_g_matrix_and_inverse(b_matrix):
    # Calculate G = B * B^T
    g_matrix = np.dot(b_matrix, b_matrix.T)

    # Diagonalize G
    eigenvalues, eigenvectors = np.linalg.eigh(g_matrix)

    # Check for near-zero eigenvalues and handle them to avoid numerical issues
    tol = 1e-10
    inv_eigenvalues = np.array([1 / ev if ev > tol else 0 for ev in eigenvalues])

    # Calculate G^-1 using the eigen decomposition
    g_inv = np.dot(eigenvectors, np.dot(np.diag(inv_eigenvalues), eigenvectors.T))

    return g_matrix, g_inv


def calculate_internal_gradient(b_matrix, g_inv, total_gradient):
    # Flatten the gradient in Cartesian coordinates
    total_gradient_flat = total_gradient.flatten()

    # Calculate g(internal coordinates) = G^-1 * B * g(xyz)
    internal_gradient = np.dot(g_inv, np.dot(b_matrix, total_gradient_flat))

    return internal_gradient



def initialize_inverse_hessian_internal(num_bonds, num_angles, num_torsions):

    # Total number of internal coordinates
    dim = num_bonds + num_angles + num_torsions

    # Initialize an empty matrix
    inverse_hessian = np.zeros((dim, dim))

    # Set diagonal elements based on type
    for i in range(num_bonds):
        inverse_hessian[i, i] = 1 / 600  # Bonds

    for i in range(num_bonds, num_bonds + num_angles):
        inverse_hessian[i, i] = 1 / 150  # Angles

    for i in range(num_bonds + num_angles, dim):
        inverse_hessian[i, i] = 1 / 80  # Torsions

    return inverse_hessian




def bfgs_geometry_optimization_internals(coordinates, atoms, bonds, k_ccc, angle0_ccc, k_other, angle0_other, a, n, parameters, k_values, r0_values, rms_gradient_threshold=1e-3, max_iterations=100):
    b_angles = calculate_angles(coordinates, bonds)
    dihedrals = calculate_dihedral_angles(coordinates, bonds)

    num_bonds = len(bonds)
    num_angles = len(b_angles)
    num_torsions = len(dihedrals)

    mq = initialize_inverse_hessian_internal(num_bonds, num_angles, num_torsions)

    for iteration in range(max_iterations):
        distances = calculate_bond_distances(coordinates, bonds)
        angles = calculate_angles(coordinates, bonds)
        dihedrals = calculate_dihedral_angles(coordinates, bonds)
        unique_distances = calculate_all_unique_distances(coordinates, bonds)

        total_gradient = calculate_total_gradient(coordinates, bonds, distances, k_values, r0_values, angles, atoms,
                                                  k_ccc, angle0_ccc, k_other, angle0_other, dihedrals, a, n,
                                                  unique_distances, parameters)


        rms_gradient = np.sqrt(np.mean(total_gradient ** 2))

        print(f"Iteration {iteration}: RMS Gradient = {rms_gradient:.6f}")

        # Check convergence
        if rms_gradient < rms_gradient_threshold:
            print("Converged: Gradient norm below threshold.")
            break

        b_matrix = calculate_b_matrix(coordinates, bonds, b_angles, dihedrals)
        _, g_inv = calculate_g_matrix_and_inverse(b_matrix)
        internal_gradient = calculate_internal_gradient(b_matrix, g_inv, total_gradient)

        pq = -np.dot(mq, internal_gradient)
        lq = np.sqrt(np.sum(pq**2) / len(pq))
        if lq > 0.020:
            pq *= 0.020 / lq

        sq = pq
        q = np.concatenate((distances[:, 2], angles[:, 3], dihedrals[:, 4])) if num_torsions > 0 else np.concatenate((distances[:, 2], angles[:, 3]))
        q_new = q + sq

        coord = coordinates.flatten()
        cartesian_thresh = 1e-5
        pi = np.pi

        dx = np.dot(b_matrix.T, np.dot(g_inv, pq)).reshape(coord.shape)
        new_coord = coord + dx
        max_cartesian_diff = np.max(np.abs(new_coord - coord))

        while max_cartesian_diff > cartesian_thresh:
            coord_mat = new_coord.reshape(-1, 3)
            update_q = np.concatenate((
                calculate_bond_distances(coord_mat, bonds)[:, 2],
                calculate_angles(coord_mat, bonds)[:, 3],
                calculate_dihedral_angles(coord_mat, bonds)[:, 4]
            )) if num_torsions > 0 else np.concatenate((
                calculate_bond_distances(coord_mat, bonds)[:, 2],
                calculate_angles(coord_mat, bonds)[:, 3]
            ))

            sq = q_new - update_q
            for i in range(num_bonds + num_angles, num_bonds + num_angles + num_torsions):
                if sq[i] > pi:
                    sq[i] -= 2 * pi
                elif sq[i] < -pi:
                    sq[i] += 2 * pi

            coord = new_coord
            dx = np.dot(b_matrix.T, np.dot(g_inv, sq)).reshape(coord.shape)
            new_coord = coord + dx
            max_cartesian_diff = np.max(np.abs(new_coord - coord))

        coordinates = new_coord.reshape(-1, 3)  # Ensure coordinates update every step

        b_angles_new = calculate_angles(coordinates, bonds)
        dihedrals_new = calculate_dihedral_angles(coordinates, bonds)
        distances_new = calculate_bond_distances(coordinates, bonds)
        unique_distances_new = calculate_all_unique_distances(coordinates, bonds)

        total_gradient_new = calculate_total_gradient(coordinates, bonds, distances_new, k_values, r0_values, b_angles_new, atoms,
                                                      k_ccc, angle0_ccc, k_other, angle0_other, dihedrals_new, a, n,
                                                      unique_distances_new, parameters)

        b_matrix_new = calculate_b_matrix(coordinates, bonds, b_angles_new, dihedrals_new)
        _, g_inv_new = calculate_g_matrix_and_inverse(b_matrix_new)
        internal_gradient_new = calculate_internal_gradient(b_matrix_new, g_inv_new, total_gradient_new)

        sq_new = update_q - q
        for i in range(num_bonds + num_angles, num_bonds + num_angles + num_torsions):
            if sq_new[i] > pi:
                sq_new[i] -= 2 * pi
            elif sq_new[i] < -pi:
                sq_new[i] += 2 * pi

        yq = internal_gradient_new - internal_gradient
        vq = np.dot(mq, yq.flatten()).reshape(yq.shape)

        sy = np.dot(sq_new.flatten(), yq.flatten())
        if sy == 0:
            print("Warning: Division by zero in Hessian update.")
            continue

        mq += ((sy + np.dot(yq.flatten(), vq.flatten())) * np.outer(sq_new, sq_new) / (sy ** 2) - (np.outer(vq, sq_new) + np.outer(sq_new, vq)) / sy)

    return coordinates



def main():
    input_folder = "inputs"
    if not os.path.exists(input_folder):
        print(f"The folder '{input_folder}' does not exist. Please create it and add input files.")
        return

    while True:
        file_name = input("Enter the name of the file (inside the 'inputs' folder): ")
        file_path = os.path.join(input_folder, file_name)

        if not os.path.isfile(file_path):
            print(f"Error: The file '{file_name}' does not exist. Please enter a valid file name.")
            continue

        break  # Exit loop when a valid file is provided

    atoms, coordinates, bonds = parse_molecule_file(file_path)
    distances = calculate_bond_distances(coordinates, bonds)
    bond_energies = calculate_bond_energies(distances, atoms, bonds)
    b_angles = calculate_angles(coordinates, bonds)

    # Define constants for angle energy calculations
    k_ccc = 60         # Force constant for C-C-C angles
    angle0_ccc = 109.5 * (math.pi/180) # Equilibrium angle for C-C-C
    k_other = 35       # Force constant for other angles
    angle0_other = 109.5 * (math.pi/180) # Equilibrium angle for other angles

    angle_energies = calculate_angle_energies(b_angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other)
    dihedrals = calculate_dihedral_angles(coordinates, bonds)

    # Dihedral lines
    a = 0.3
    n = 3
    dihedral_energies = calculate_dihedral_energies(dihedrals, a, n)
    unique_distances = calculate_all_unique_distances(coordinates, bonds)


    # Define Lennard-Jones parameters
    parameters = {
        ("H", "H"): (4382.4417, 22.932357),
        ("H", "C"): (64393.99285, 108.6444378),
        ("C", "H"): (64393.99285, 108.6444378),
        ("C", "C"): (946181.7423, 514.714375)
    }

    lj_energies = calculate_lennard_jones_energy(unique_distances, atoms, parameters)

    # Calculate the B-matrix
    b_matrix = calculate_b_matrix(coordinates, bonds, b_angles, dihedrals)

    # Define force constants (k) and equilibrium bond lengths (r0)
    k_values = [300 if atoms[bond[0] - 1]["element"] == 'C' and atoms[bond[1] - 1]["element"] == 'C' else 350 for bond in bonds]
    r0_values = [1.53 if atoms[bond[0] - 1]["element"] == 'C' and atoms[bond[1] - 1]["element"] == 'C' else 1.11 for bond in bonds]

    stretching_gradient = calculate_stretching_gradient(coordinates, bonds, distances, k_values, r0_values)

    bending_gradient = calculate_bending_gradient(coordinates, atoms, b_angles, k_ccc, angle0_ccc, k_other, angle0_other)
    dihedral_gradient = calculate_torsional_energy_gradient(coordinates, dihedrals, a, n)
    vdw_gradient = calculate_vdw_energy_gradient(coordinates, unique_distances, atoms, parameters)

    # Calculate total potential energy and its components
    total_potential_energy, total_stretching_energy, total_bond_energy, total_dihedral_energy, total_lj_energy = calculate_total_potential_energy(
        bond_energies, angle_energies, dihedral_energies, lj_energies
    )


    total_gradient = calculate_total_gradient(coordinates, bonds, distances, k_values, r0_values, b_angles, atoms, k_ccc, angle0_ccc, k_other, angle0_other, dihedrals, a, n, unique_distances, parameters)

    num_atoms = len(coordinates)
    m = initialize_inverse_hessian(num_atoms)

#    final_coordinates = bfgs_geometry_optimization(coordinates, atoms, bonds, k_ccc, angle0_ccc, k_other, angle0_other, a, n, parameters, k_values, r0_values, max_iterations=100, rms_gradient_threshold=0.001)

    g_matrix, g_inv = calculate_g_matrix_and_inverse(b_matrix)
    internal_gradient = calculate_internal_gradient(b_matrix, g_inv, total_gradient)

    num_bonds = len(bonds)
    num_angles = len(b_angles)
    num_torsions = len(dihedrals)

    inverse_hessian_internal = initialize_inverse_hessian_internal(num_bonds, num_angles, num_torsions)

    final_coordinates_internal = bfgs_geometry_optimization_internals(coordinates, atoms, bonds, k_ccc, angle0_ccc, k_other, angle0_other, a, n, parameters, k_values, r0_values, rms_gradient_threshold=1e-3, max_iterations=100)

    if atoms:
        print("Atomic Coordinates:")
        for atom in atoms:
            print(f"Atom: {atom['element']}, Coordinates: {atom['coordinates']}")
        print("\nCoordinates once again:")
        print(coordinates)
        print("\nBonds:")
        print(bonds)
        print("\nDistances:")
        print(distances)
        print("\nBond energies:")
        print(bond_energies)
        print("\nAngles:")
        print(b_angles)
        print("\nAngle energies:")
        print(angle_energies)
        print("\nDihedrals:")
        print(dihedrals)
        print("\nDihedral energies:")
        print(dihedral_energies)
        print("\nUnique Non-Bonded Distances:")
        print(unique_distances)
        print("\nLennard-Jones Energies:")
        print(lj_energies)
        print("\nB-matrix:")
        print(b_matrix)
        print("\nStretching Gradient:")
        print(stretching_gradient)
        print("\nBending Gradient:")
        print(bending_gradient)
        print("\nTorsional Gradient:")
        print(dihedral_gradient)
        print("\nVDW Gradient:")
        print(vdw_gradient)

        print("\nTotal Potential Energy and its Components (kcal/mol):")
        print(f" - Total Stretching Energy: {total_stretching_energy}")
        print(f" - Total Bond Energy: {total_bond_energy}")
        print(f" - Total Dihedral Energy: {total_dihedral_energy}")
        print(f" - Total Lennard-Jones Energy: {total_lj_energy}")
        print(f" - Total Potential Energy: {total_potential_energy}")

        print("\nTotal Gradient:")
        print(total_gradient)
        print("\nInverse Hessian:")
        print(m)

#        print("Optimized Coordinates:")
#        print(final_coordinates)

        print("\nG-matrix:")
        print(g_matrix)
        print("\nInverse of G-matrix:")
        print(g_inv)
        print("\nGradient in internal coordinates:")
        print(internal_gradient)

        print("Inverse Hessian (Internal Coordinates):")
        print(inverse_hessian_internal)

        print("Optimized Coordinates (Internal):")
        print(final_coordinates_internal)

    else:
        print("No valid atomic data found in the file.")


if __name__ == "__main__":
    main()
