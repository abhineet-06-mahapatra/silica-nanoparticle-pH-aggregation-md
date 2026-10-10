"""
Salt concentration study, pH 10 fixed.
Reuses the existing pH 10 charged nanoparticle (np_pH10.npz) and box geometry
(system_geom_pH10.npz) from the pH sweep, and adds background NaCl or CaCl2
at two concentrations each, on top of the 234 Na+ ions already neutralizing
the nanoparticle surface charge.

4 systems: NaCl_0.05M, NaCl_0.3M, CaCl2_0.05M, CaCl2_0.3M
"""

import numpy as np
import subprocess
from scipy.spatial import cKDTree
from ase.io import read

AVOGADRO = 6.02214076e23

# --- extended atom types (added to the existing TYPE_MAP/MASSES from the pH sweep) ---
TYPE_MAP_EXT = {'Si':1, 'Ob':2, 'Oh':3, 'Ho':4, 'Ow':5, 'Hw':6, 'Na':7, 'Od':8, 'Cl':9, 'Ca':10}
MASSES_EXT = {1:28.0855, 2:15.9994, 3:15.9994, 4:1.008, 5:15.9994, 6:1.008,
              7:22.9898, 8:15.9994, 9:35.453, 10:40.078}
CHARGE_OW, CHARGE_HW, CHARGE_NA, CHARGE_CL, CHARGE_CA = -0.82, 0.41, 1.0, -1.0, 2.0
LJ_EXTRA = {9: (0.1000, 4.4000), 10: (0.1000, 2.8720)}   # (eps kcal/mol, sigma A)


def n_ions_for_concentration(molarity, box_volume_A3):
    volume_L = box_volume_A3 * 1e-27
    n_moles = molarity * volume_L
    return int(round(n_moles * AVOGADRO))


def write_ion_pdb(symbol, path):
    with open(path, 'w') as f:
        f.write(f"HETATM{1:>5} {symbol:<3}  {symbol.upper():<3} A{1:>4}    "
                f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {symbol:>2}\n")
        f.write("END\n")


def run_packmol_with_salt(label, box_size, np1_pos, np2_pos, n_water,
                            n_na_neutralize, n_cation_salt, n_anion_salt, cation_symbol):
    np1_center = np1_pos.mean(axis=0)
    np2_center = np2_pos.mean(axis=0)
    exclude_radius = 16.0

    structures = f"""
structure spc.pdb
  number {n_water}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure

structure na.pdb
  number {n_na_neutralize}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure

structure {cation_symbol.lower()}.pdb
  number {n_cation_salt}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure

structure cl.pdb
  number {n_anion_salt}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure
"""
    inp = f"tolerance 2.0\nfiletype pdb\noutput water_salt_{label}.pdb\n" + structures
    inp_path = f'packmol_{label}.inp'
    with open(inp_path, 'w') as f:
        f.write(inp)
    print(f"{label}: running Packmol (water={n_water}, Na+={n_na_neutralize}, "
          f"{cation_symbol}+={n_cation_salt}, Cl-={n_anion_salt})...")
    result = subprocess.run(['packmol'],
                              stdin=open(inp_path), capture_output=True, text=True)
    print(result.stdout[-300:])


def build_salt_data_file(label, box_size, np1_pos, np2_pos, np_symbols, np_charges, np_types,
                           n_na_neutralize, cation_symbol, cation_charge, cation_type_id):
    n_np_atoms_each = len(np_symbols)
    wat_ions = read(f'water_salt_{label}.pdb')
    wi_pos = wat_ions.get_positions()
    wi_sym = np.array(wat_ions.get_chemical_symbols())

    # Species written in Packmol order: water(O,H,H) -> Na -> cation(Na or Ca) -> Cl
    # Safer: identify by element symbol directly
    is_o = wi_sym == 'O'
    is_h = wi_sym == 'H'
    is_na = wi_sym == 'Na'
    is_cl = wi_sym == 'Cl'
    is_ca = wi_sym == 'Ca'

    # water: consecutive O,H,H triplets as written by packmol for the spc structure block
    water_pos = wi_pos[is_o | is_h]
    n_water_atoms = len(water_pos)

    if cation_symbol == 'Na':
        # all Na atoms = neutralizing Na+ (234) + salt Na+ combined, since same element
        na_all_pos = wi_pos[is_na]
        na_neutralize_pos = na_all_pos[:n_na_neutralize]
        cation_salt_pos = na_all_pos[n_na_neutralize:]
    else:
        na_neutralize_pos = wi_pos[is_na]
        cation_salt_pos = wi_pos[is_ca]

    cl_pos = wi_pos[is_cl]

    all_pos, all_type, all_charge, all_mol = [], [], [], []
    mol_id = 1
    for i in range(n_np_atoms_each):
        all_pos.append(np1_pos[i]); all_type.append(TYPE_MAP_EXT[np_types[i]])
        all_charge.append(np_charges[i]); all_mol.append(mol_id)
    mol_id += 1
    for i in range(n_np_atoms_each):
        all_pos.append(np2_pos[i]); all_type.append(TYPE_MAP_EXT[np_types[i]])
        all_charge.append(np_charges[i]); all_mol.append(mol_id)
    mol_id += 1

    np_total_atoms = 2 * n_np_atoms_each
    water_bonds, water_angles = [], []
    for w in range(n_water_atoms // 3):
        o_i, h1_i, h2_i = 3*w, 3*w+1, 3*w+2
        all_pos.append(water_pos[o_i]);  all_type.append(TYPE_MAP_EXT['Ow']); all_charge.append(CHARGE_OW); all_mol.append(mol_id)
        all_pos.append(water_pos[h1_i]); all_type.append(TYPE_MAP_EXT['Hw']); all_charge.append(CHARGE_HW); all_mol.append(mol_id)
        all_pos.append(water_pos[h2_i]); all_type.append(TYPE_MAP_EXT['Hw']); all_charge.append(CHARGE_HW); all_mol.append(mol_id)
        base = np_total_atoms + 3*w + 1
        water_bonds.append((base, base+1)); water_bonds.append((base, base+2))
        water_angles.append((base+1, base, base+2))
        mol_id += 1

    for p in na_neutralize_pos:
        all_pos.append(p); all_type.append(TYPE_MAP_EXT['Na']); all_charge.append(CHARGE_NA); all_mol.append(mol_id); mol_id += 1
    for p in cation_salt_pos:
        all_pos.append(p); all_type.append(cation_type_id); all_charge.append(cation_charge); all_mol.append(mol_id); mol_id += 1
    for p in cl_pos:
        all_pos.append(p); all_type.append(TYPE_MAP_EXT['Cl']); all_charge.append(CHARGE_CL); all_mol.append(mol_id); mol_id += 1

    all_pos = np.array(all_pos); all_type = np.array(all_type)
    all_charge = np.array(all_charge); all_mol = np.array(all_mol)
    n_total_atoms = len(all_pos)

    # pH 10 has no silanol H left (99.9% deprotonated), so no NP bonds are needed

    fname = f'system_{label}.data'
    with open(fname, 'w') as f:
        f.write(f"LAMMPS data file - Silica NP salt study, {label}\n\n")
        f.write(f"{n_total_atoms} atoms\n{len(water_bonds)} bonds\n{len(water_angles)} angles\n\n")
        f.write(f"{len(MASSES_EXT)} atom types\n1 bond types\n1 angle types\n\n")
        f.write(f"0.0 {box_size[0]:.4f} xlo xhi\n0.0 {box_size[1]:.4f} ylo yhi\n0.0 {box_size[2]:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t, m in MASSES_EXT.items():
            f.write(f"{t} {m}\n")
        f.write("\nBond Coeffs\n\n1 554.13 1.0\n\nAngle Coeffs\n\n1 45.77 109.47\n\n")
        f.write("Atoms\n\n")
        for i in range(n_total_atoms):
            x, y, z = all_pos[i]
            f.write(f"{i+1} {all_mol[i]} {all_type[i]} {all_charge[i]:.4f} {x:.4f} {y:.4f} {z:.4f}\n")
        f.write("\nBonds\n\n")
        for i, (a, b) in enumerate(water_bonds):
            f.write(f"{i+1} 1 {a} {b}\n")
        f.write("\nAngles\n\n")
        for i, (a, b, c) in enumerate(water_angles):
            f.write(f"{i+1} 1 {a} {b} {c}\n")

    print(f"{label}: total atoms={n_total_atoms}, total charge={all_charge.sum():.4f}, saved {fname}")


def write_salt_lammps_input(label, n_np_atoms_each, equil_steps=10000, prod_steps=25000):
    """LAMMPS input for a salt-study system. Same protocol as the pH sweep, with shorter runs.
    NOTE: the SHAKE line is essential - without it SPC water (rigid) blows up in minimization."""
    n1_end, n2_start, n2_end = n_np_atoms_each, n_np_atoms_each + 1, 2 * n_np_atoms_each
    script = f"""
units           real
atom_style      full
boundary        p p p
pair_style      lj/cut/coul/long 10.0
bond_style      harmonic
angle_style     harmonic
kspace_style    pppm 1.0e-4
read_data       system_{label}.data

pair_coeff 1 1  1.8405e-6 3.7064
pair_coeff 2 2  0.1554    3.5532
pair_coeff 3 3  0.1554    3.5532
pair_coeff 4 4  0.0       0.0
pair_coeff 5 5  0.1554    3.1660
pair_coeff 6 6  0.0       0.0
pair_coeff 7 7  0.1301    2.6378
pair_coeff 8 8  0.1554    3.5532
pair_coeff 9 9  0.1000    4.4000
pair_coeff 10 10 0.1000   2.8720
pair_modify     mix arithmetic

fix             rigidwater all shake 0.0001 20 0 b 1 a 1
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
group           np1_atoms id 1:{n1_end}
group           np2_atoms id {n2_start}:{n2_end}

thermo          100
thermo_style    custom step temp epair etotal press
min_style       cg
minimize        1.0e-4 1.0e-6 500 5000

velocity        all create 300.0 12345 mom yes rot yes dist gaussian
timestep        1.0
fix             nve_all all nve
fix             langevin_equil all langevin 10.0 300.0 100.0 12345
thermo          200
run             {equil_steps}
unfix           langevin_equil
unfix           nve_all

compute         com_np1 np1_atoms com
compute         com_np2 np2_atoms com
variable        dx equal "c_com_np2[1] - c_com_np1[1]"
variable        dy equal "c_com_np2[2] - c_com_np1[2]"
variable        dz equal "c_com_np2[3] - c_com_np1[3]"
variable        np_distance equal "sqrt(v_dx^2 + v_dy^2 + v_dz^2)"
fix             nvt_prod all nvt temp 300.0 300.0 100.0
fix             printdist all print 100 "$(step) ${{np_distance}}" file np_distance_{label}.dat screen no
dump            traj all custom 1000 traj_{label}.lammpstrj id type x y z
dump_modify     traj sort id
thermo          500
run             {prod_steps}
write_data      final_{label}.data
"""
    assert 'fix             rigidwater' in script
    with open(f'in.{label}', 'w') as f:
        f.write(script)
    print(f"Saved in.{label}")


if __name__ == '__main__':
    # --- load existing pH 10 nanoparticle + box geometry (already built) ---
    geom = np.load('system_geom_pH10.npz')
    box_size = geom['box_size']; np1_pos = geom['np1_pos']; np2_pos = geom['np2_pos']
    np_symbols = geom['symbols']; np_charges = geom['charges']
    n_na_neutralize = int(geom['n_na_total'])   # 468 total (234/NP), unchanged by salt
    np_types = np.load('np_pH10.npz')['types']

    write_ion_pdb('Na', 'na.pdb')
    write_ion_pdb('Cl', 'cl.pdb')
    write_ion_pdb('Ca', 'ca.pdb')

    box_volume = np.prod(box_size)
    np_volume = 2 * (4/3)*np.pi*(15.0**3)
    available_volume = box_volume - np_volume
    n_water = int(available_volume * 0.0334)

    concentrations = [0.05, 0.3]   # mol/L

    for conc in concentrations:
        # --- NaCl: adds equal Na+ and Cl- ---
        n_salt = n_ions_for_concentration(conc, available_volume)
        label = f'pH10_NaCl_{conc}M'
        run_packmol_with_salt(label, box_size, np1_pos, np2_pos, n_water,
                                n_na_neutralize, n_salt, n_salt, 'Na')
        build_salt_data_file(label, box_size, np1_pos, np2_pos, np_symbols, np_charges, np_types,
                               n_na_neutralize, 'Na', CHARGE_NA, TYPE_MAP_EXT['Na'])
        write_salt_lammps_input(label, len(np_symbols))

        # --- CaCl2: adds n_salt Ca2+ and 2*n_salt Cl- ---
        label = f'pH10_CaCl2_{conc}M'
        run_packmol_with_salt(label, box_size, np1_pos, np2_pos, n_water,
                                n_na_neutralize, n_salt, 2*n_salt, 'Ca')
        build_salt_data_file(label, box_size, np1_pos, np2_pos, np_symbols, np_charges, np_types,
                               n_na_neutralize, 'Ca', CHARGE_CA, TYPE_MAP_EXT['Ca'])
        write_salt_lammps_input(label, len(np_symbols))

