"""
Step 5-6: Place two copies of the charged nanoparticle in a box, fill with
SPC water (+ Na+ counterions) via Packmol, and write a LAMMPS full-style
data file for a given pH condition.

Requires: np_pH<X>.npz (from 02_assign_pH_charges.py), Packmol on PATH (the `packmol` command),
a spc.pdb and na.pdb template (written here).

Output: system_pH<X>.data
"""

import numpy as np
import subprocess
from scipy.spatial import cKDTree
from ase.io import read, write
from ase import Atoms

PH_VALUES = [2, 4, 7, 8, 10, 12]
SEPARATION = 40.0   # Angstrom, center-to-center NP separation
PADDING = 6.0        # Angstrom padding around the NP pair
CHARGE_NA = 1.0

TYPE_MAP = {'Si':1, 'Ob':2, 'Oh':3, 'Ho':4, 'Ow':5, 'Hw':6, 'Na':7, 'Od':8}
MASSES = {1:28.0855, 2:15.9994, 3:15.9994, 4:1.008, 5:15.9994, 6:1.008, 7:22.9898, 8:15.9994}
CHARGE_OW, CHARGE_HW = -0.82, 0.41


def write_templates():
    with open('spc.pdb', 'w') as f:
        f.write("HETATM    1  OW  SPC A   1       0.000   0.000   0.000  1.00  0.00           O\n")
        f.write("HETATM    2  HW1 SPC A   1       0.816   0.577   0.000  1.00  0.00           H\n")
        f.write("HETATM    3  HW2 SPC A   1      -0.816   0.577   0.000  1.00  0.00           H\n")
        f.write("END\n")
    # Note: correct fixed-column PDB formatting is essential here - the element
    # symbol must land in columns 77-78, or downstream parsers misread it (see
    # the Na/'N' parsing bug documented in the report, Section 6).
    with open('na.pdb', 'w') as f:
        f.write(f"HETATM{1:>5} NA   NA  A{1:>4}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          Na\n")
        f.write("END\n")


def place_two_nps_and_box(pH_val):
    data = np.load(f'np_pH{pH_val}.npz')
    positions = data['positions']; symbols = data['symbols']
    charges = data['charges']; n_na = int(data['n_na'])

    com = positions.mean(axis=0)
    pos_centered = positions - com
    np1_pos = pos_centered.copy()
    np2_pos = pos_centered.copy() + np.array([SEPARATION, 0, 0])

    all_pos = np.vstack([np1_pos, np2_pos])
    box_min = all_pos.min(axis=0) - PADDING
    box_max = all_pos.max(axis=0) + PADDING
    box_size = box_max - box_min
    shift = -box_min
    np1_pos += shift; np2_pos += shift

    np.savez(f'system_geom_pH{pH_val}.npz', box_size=box_size, np1_pos=np1_pos,
              np2_pos=np2_pos, symbols=symbols, charges=charges, n_na_total=2*n_na)
    return box_size, np1_pos, np2_pos, n_na


def run_packmol(pH_val, box_size, np1_pos, np2_pos, n_na_total):
    box_volume = np.prod(box_size)
    np_volume = 2 * (4/3)*np.pi*(15.0**3)
    n_water = int((box_volume - np_volume) * 0.0334)

    np1_center = np1_pos.mean(axis=0); np2_center = np2_pos.mean(axis=0)
    exclude_radius = 16.0

    inp = f"""
tolerance 2.0
filetype pdb
output water_box_pH{pH_val}.pdb

structure spc.pdb
  number {n_water}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure

"""
    if n_na_total > 0:   # pH 2 and pH 4: surface is neutral, no counterions needed
        inp += f"""
structure na.pdb
  number {n_na_total}
  inside box 0. 0. 0. {box_size[0]:.2f} {box_size[1]:.2f} {box_size[2]:.2f}
  outside sphere {np1_center[0]:.2f} {np1_center[1]:.2f} {np1_center[2]:.2f} {exclude_radius}
  outside sphere {np2_center[0]:.2f} {np2_center[1]:.2f} {np2_center[2]:.2f} {exclude_radius}
end structure
"""
    inp_path = f'packmol_pH{pH_val}.inp'
    with open(inp_path, 'w') as f:
        f.write(inp)
    subprocess.run(['packmol'], stdin=open(inp_path), capture_output=True, text=True)


def build_data_file(pH_val):
    geom = np.load(f'system_geom_pH{pH_val}.npz')
    box_size = geom['box_size']; np1_pos = geom['np1_pos']; np2_pos = geom['np2_pos']
    np_symbols = geom['symbols']; np_charges = geom['charges']
    n_na_expected = int(geom['n_na_total'])
    np_types = np.load(f'np_pH{pH_val}.npz')['types']

    wat_ions = read(f'water_box_pH{pH_val}.pdb')
    wi_pos = wat_ions.get_positions(); wi_sym = np.array(wat_ions.get_chemical_symbols())
    na_mask = wi_sym == 'Na'
    n_na = int(na_mask.sum())
    water_pos = wi_pos[~na_mask]; na_pos = wi_pos[na_mask]
    n_water_atoms = len(water_pos)

    assert n_na == n_na_expected, f"Na mismatch: found {n_na}, expected {n_na_expected}"

    n_np_atoms_each = len(np_symbols)
    all_pos, all_type, all_charge, all_mol = [], [], [], []
    mol_id = 1
    for i in range(n_np_atoms_each):
        all_pos.append(np1_pos[i]); all_type.append(TYPE_MAP[np_types[i]])
        all_charge.append(np_charges[i]); all_mol.append(mol_id)
    mol_id += 1
    for i in range(n_np_atoms_each):
        all_pos.append(np2_pos[i]); all_type.append(TYPE_MAP[np_types[i]])
        all_charge.append(np_charges[i]); all_mol.append(mol_id)
    mol_id += 1

    np_total_atoms = 2 * n_np_atoms_each
    water_bonds, water_angles = [], []
    for w in range(n_water_atoms // 3):
        o_i, h1_i, h2_i = 3*w, 3*w+1, 3*w+2
        all_pos.append(water_pos[o_i]);  all_type.append(TYPE_MAP['Ow']); all_charge.append(CHARGE_OW); all_mol.append(mol_id)
        all_pos.append(water_pos[h1_i]); all_type.append(TYPE_MAP['Hw']); all_charge.append(CHARGE_HW); all_mol.append(mol_id)
        all_pos.append(water_pos[h2_i]); all_type.append(TYPE_MAP['Hw']); all_charge.append(CHARGE_HW); all_mol.append(mol_id)
        base = np_total_atoms + 3*w + 1
        water_bonds.append((base, base+1)); water_bonds.append((base, base+2))
        water_angles.append((base+1, base, base+2))
        mol_id += 1

    for i in range(n_na):
        all_pos.append(na_pos[i]); all_type.append(TYPE_MAP['Na']); all_charge.append(CHARGE_NA); all_mol.append(mol_id)
        mol_id += 1

    all_pos = np.array(all_pos); all_type = np.array(all_type)
    all_charge = np.array(all_charge); all_mol = np.array(all_mol)
    n_total_atoms = len(all_pos)

    silanol_bonds = []
    for np_offset, pos in [(0, np1_pos), (n_np_atoms_each, np2_pos)]:
        o_indices = [i for i in range(n_np_atoms_each) if np_types[i] == 'Oh']
        h_indices = [i for i in range(n_np_atoms_each) if np_types[i] == 'Ho']
        if len(h_indices) == 0:
            continue
        o_pos = pos[o_indices]; h_pos = pos[h_indices]
        tree = cKDTree(h_pos)
        for k, oi in enumerate(o_indices):
            dist, idx = tree.query(o_pos[k], k=1)
            if dist < 1.3:
                silanol_bonds.append((np_offset+oi+1, np_offset+h_indices[idx]+1))

    all_bonds = silanol_bonds + water_bonds
    all_angles = water_angles

    fname = f'system_pH{pH_val}.data'
    with open(fname, 'w') as f:
        f.write(f"LAMMPS data file - Silica NP aggregation, pH {pH_val}\n\n")
        f.write(f"{n_total_atoms} atoms\n{len(all_bonds)} bonds\n{len(all_angles)} angles\n\n")
        f.write(f"{len(MASSES)} atom types\n1 bond types\n1 angle types\n\n")
        f.write(f"0.0 {box_size[0]:.4f} xlo xhi\n0.0 {box_size[1]:.4f} ylo yhi\n0.0 {box_size[2]:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t, m in MASSES.items():
            f.write(f"{t} {m}\n")
        f.write("\nBond Coeffs\n\n1 554.13 1.0\n\nAngle Coeffs\n\n1 45.77 109.47\n\n")
        f.write("Atoms\n\n")
        for i in range(n_total_atoms):
            x, y, z = all_pos[i]
            f.write(f"{i+1} {all_mol[i]} {all_type[i]} {all_charge[i]:.4f} {x:.4f} {y:.4f} {z:.4f}\n")
        f.write("\nBonds\n\n")
        for i, (a, b) in enumerate(all_bonds):
            f.write(f"{i+1} 1 {a} {b}\n")
        f.write("\nAngles\n\n")
        for i, (a, b, c) in enumerate(all_angles):
            f.write(f"{i+1} 1 {a} {b} {c}\n")

    print(f"pH {pH_val}: total atoms={n_total_atoms}, total charge={all_charge.sum():.4f}, "
          f"silanol bonds={len(silanol_bonds)}, saved {fname}")


if __name__ == '__main__':
    write_templates()
    for pH_val in PH_VALUES:
        box_size, np1_pos, np2_pos, n_na = place_two_nps_and_box(pH_val)
        run_packmol(pH_val, box_size, np1_pos, np2_pos, 2*n_na)
        build_data_file(pH_val)
