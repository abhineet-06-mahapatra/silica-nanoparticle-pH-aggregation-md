"""
Step 4: Assign CLAYFF partial charges to the passivated nanoparticle,
and apply a pH-dependent Henderson-Hasselbalch deprotonation model.

Input:  nanoparticle_final.xyz
Output: np_pH<X>.npz  (positions, symbols, charges, types, n_na) for pH in [2, 7, 10]
"""

import numpy as np
from scipy.spatial import cKDTree
from ase.io import read

SI_O_CUTOFF = 2.0
PKA = 6.8
PH_VALUES = [2, 4, 7, 8, 10, 12]

CHARGE_SI, CHARGE_OB, CHARGE_OH, CHARGE_HO, CHARGE_OD = 2.1, -1.05, -0.95, 0.425, -1.525


def build_ph_system(pH, seed=42, xyz_path='nanoparticle_final.xyz'):
    nanoparticle_final = read(xyz_path)
    positions = nanoparticle_final.get_positions()
    symbols = np.array(nanoparticle_final.get_chemical_symbols())
    n_atoms = len(nanoparticle_final)

    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=SI_O_CUTOFF, output_type='ndarray')
    si_o_neighbors = {i: [] for i in range(n_atoms)}
    o_h_partner = {}
    for i, j in pairs:
        si_i, si_j = symbols[i], symbols[j]
        if {si_i, si_j} == {'Si', 'O'}:
            o, si = (i, j) if si_i == 'O' else (j, i)
            si_o_neighbors[o].append(si)
        elif {si_i, si_j} == {'O', 'H'}:
            o, h = (i, j) if si_i == 'O' else (j, i)
            o_h_partner[o] = h

    o_idx = np.where(symbols == 'O')[0]
    bridging_o = set(o for o in o_idx if len(si_o_neighbors[o]) == 2)
    silanol_o = [o for o in o_idx if len(si_o_neighbors[o]) == 1 and o in o_h_partner]

    frac_deprot = 1.0 / (1.0 + 10**(PKA - pH))
    n_deprot = int(round(frac_deprot * len(silanol_o)))
    rng = np.random.default_rng(seed)
    deprot_sites = set(rng.choice(silanol_o, size=n_deprot, replace=False)) if n_deprot > 0 else set()

    charges = np.zeros(n_atoms)
    types = []
    keep_mask = np.ones(n_atoms, dtype=bool)

    for i in range(n_atoms):
        if symbols[i] == 'Si':
            charges[i] = CHARGE_SI; types.append('Si')
        elif symbols[i] == 'H':
            types.append('Ho')
        elif symbols[i] == 'O':
            if i in bridging_o:
                charges[i] = CHARGE_OB; types.append('Ob')
            elif i in silanol_o:
                if i in deprot_sites:
                    charges[i] = CHARGE_OD; types.append('Od')
                    keep_mask[o_h_partner[i]] = False
                else:
                    charges[i] = CHARGE_OH; types.append('Oh')
    types = np.array(types)

    for o in silanol_o:
        if o not in deprot_sites:
            charges[o_h_partner[o]] = CHARGE_HO

    positions_kept = positions[keep_mask]
    symbols_kept = symbols[keep_mask]
    charges_kept = charges[keep_mask]
    types_kept = types[keep_mask]

    print(f"pH {pH}: {n_deprot}/{len(silanol_o)} deprotonated ({100*frac_deprot:.1f}%), "
          f"{len(symbols_kept)} atoms/NP, net charge/NP = {charges_kept.sum():.3f}, Na+/NP = {n_deprot}")

    return positions_kept, symbols_kept, charges_kept, types_kept, n_deprot


if __name__ == '__main__':
    for pH_val in PH_VALUES:
        pos, sym, chg, typ, n_na = build_ph_system(pH_val)
        np.savez(f'np_pH{pH_val}.npz', positions=pos, symbols=sym, charges=chg, types=typ, n_na=n_na)
