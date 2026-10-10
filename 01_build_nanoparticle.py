"""
Step 2-3: Build a 3 nm silica nanoparticle from an alpha-quartz CIF.
- Replicate unit cell into a supercell, carve a 3 nm diameter sphere
- Identify under-coordinated surface atoms and passivate with H (silanol formation)
- Fix any over-coordinated atoms introduced during capping

Input:  quartz_relaxed.cif  (alpha-quartz unit cell)
Output: nanoparticle_final.xyz  (fully passivated 3 nm silica nanoparticle)
"""

import numpy as np
from scipy.spatial import cKDTree
from ase.io import read, write
from ase.build import make_supercell
from ase import Atoms

# --- Parameters ---
CIF_PATH     = 'quartz_relaxed.cif'
TARGET_DIAM  = 30.0   # Angstrom (3 nm)
PAD          = 8.0
SI_O_CUTOFF  = 2.0    # Angstrom, real Si-O bond ~1.6 A, safe margin
SI_O_BOND_LEN = 1.62
O_H_BOND_LEN  = 1.00
COS_A = 1.0/np.sqrt(3)
SIN_A = np.sqrt(2.0/3.0)


def build_supercell_sphere():
    quartz = read(CIF_PATH)
    cell_lengths = quartz.cell.cellpar()[:3]
    reps = [int(np.ceil((TARGET_DIAM + PAD) / L)) for L in cell_lengths]
    P = np.diag(reps)
    bulk = make_supercell(quartz, P)

    center = bulk.get_center_of_mass()
    positions = bulk.get_positions()
    distances = np.linalg.norm(positions - center, axis=1)
    mask = distances <= TARGET_DIAM / 2.0
    return bulk[mask]


def compute_coordination(atoms, cutoff):
    positions = atoms.get_positions()
    symbols = np.array(atoms.get_chemical_symbols())
    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    coord = np.zeros(len(atoms), dtype=int)
    neighbor_map = {i: [] for i in range(len(atoms))}
    for i, j in pairs:
        if symbols[i] != symbols[j]:
            coord[i] += 1; coord[j] += 1
            neighbor_map[i].append(j); neighbor_map[j].append(i)
    return coord, neighbor_map, symbols


def missing_tetrahedral_dirs(existing_vecs):
    k = len(existing_vecs)
    if k == 3:
        v_sum = np.sum(existing_vecs, axis=0)
        return [-v_sum / np.linalg.norm(v_sum)]
    elif k == 2:
        v1, v2 = existing_vecs
        b1 = (v1 + v2); b1 = b1 / np.linalg.norm(b1)
        n_perp = np.cross(v1, v2); n_perp = n_perp / np.linalg.norm(n_perp)
        return [-COS_A*b1 + SIN_A*n_perp, -COS_A*b1 - SIN_A*n_perp]
    raise ValueError(f"Unexpected coordination {k}")


def main():
    nanoparticle = build_supercell_sphere()
    coord, neighbor_map, symbols = compute_coordination(nanoparticle, SI_O_CUTOFF)

    si_idx = np.where(symbols == 'Si')[0]
    o_idx  = np.where(symbols == 'O')[0]

    # remove debris (single-bond Si fragments, fully isolated O)
    delete_mask = np.zeros(len(nanoparticle), dtype=bool)
    delete_mask[si_idx[coord[si_idx] == 1]] = True
    delete_mask[o_idx[coord[o_idx] == 0]] = True
    nanoparticle = nanoparticle[~delete_mask]

    coord, neighbor_map, symbols = compute_coordination(nanoparticle, SI_O_CUTOFF)
    positions = nanoparticle.get_positions()
    si_idx = np.where(symbols == 'Si')[0]
    o_idx  = np.where(symbols == 'O')[0]

    new_positions = list(positions)
    new_symbols = list(symbols)
    new_O_from_Si = []

    for si in si_idx:
        if coord[si] in (2, 3):
            si_pos = positions[si]
            existing_vecs = [(positions[o]-si_pos)/np.linalg.norm(positions[o]-si_pos)
                              for o in neighbor_map[si]]
            for d in missing_tetrahedral_dirs(existing_vecs):
                new_o_pos = si_pos + d * SI_O_BOND_LEN
                new_positions.append(new_o_pos); new_symbols.append('O')
                new_O_from_Si.append((len(new_positions)-1, si))

    h_positions = []
    for o in o_idx:
        if coord[o] == 1:
            si_partner = neighbor_map[o][0]
            direction = positions[o] - positions[si_partner]
            direction /= np.linalg.norm(direction)
            h_positions.append(positions[o] + direction * O_H_BOND_LEN)

    for new_idx, si_parent in new_O_from_Si:
        o_pos = np.array(new_positions[new_idx]); si_pos = positions[si_parent]
        direction = (o_pos - si_pos); direction /= np.linalg.norm(direction)
        h_positions.append(o_pos + direction * O_H_BOND_LEN)

    all_positions = new_positions + h_positions
    all_symbols = new_symbols + ['H']*len(h_positions)
    nanoparticle_capped = Atoms(symbols=all_symbols, positions=all_positions)

    # --- fix any over-coordinated O/Si from capping geometry (surface curvature artifacts) ---
    positions2 = nanoparticle_capped.get_positions()
    symbols2 = np.array(nanoparticle_capped.get_chemical_symbols())

    def get_coord(pos, sym, only_si_o=False):
        tree = cKDTree(pos)
        pairs = tree.query_pairs(r=SI_O_CUTOFF, output_type='ndarray')
        c = np.zeros(len(sym), dtype=int); nm = {i: [] for i in range(len(sym))}
        for i, j in pairs:
            ok = ({sym[i], sym[j]} == {'Si', 'O'}) if only_si_o else (sym[i] != sym[j])
            if ok:
                c[i] += 1; c[j] += 1; nm[i].append(j); nm[j].append(i)
        return c, nm

    coord2, nm2 = get_coord(positions2, symbols2)
    o_idx2 = np.where(symbols2 == 'O')[0]
    overcoord_o = o_idx2[coord2[o_idx2] > 2]

    for o in overcoord_o:
        for _ in range(20):
            coord2, nm2 = get_coord(positions2, symbols2)
            if coord2[o] <= 2:
                break
            si_neighbors = nm2[o]
            dists = [np.linalg.norm(positions2[o]-positions2[si]) for si in si_neighbors]
            parent = si_neighbors[np.argmin(dists)]
            direction = positions2[o] - positions2[parent]; direction /= np.linalg.norm(direction)
            positions2[o] += direction * 0.15
            h_near = [i for i in range(len(symbols2)) if symbols2[i]=='H'
                      and np.linalg.norm(positions2[i]-positions2[o]) < 1.2]
            for h in h_near:
                positions2[h] += direction * 0.15

    nanoparticle_capped.set_positions(positions2)
    write('nanoparticle_final.xyz', nanoparticle_capped)
    print(f"Final nanoparticle: {len(nanoparticle_capped)} atoms, "
          f"formula {nanoparticle_capped.get_chemical_formula()}")


if __name__ == '__main__':
    main()
