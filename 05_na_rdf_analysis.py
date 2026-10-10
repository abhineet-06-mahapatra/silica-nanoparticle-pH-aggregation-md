"""
Radial distribution function (RDF) of Na+ ions around the nanoparticle surface,
comparing pH 7 (61.3% deprotonated) vs pH 10 (99.9% deprotonated).

Tests the counterion condensation / charge renormalization hypothesis:
if pH 10 shows a taller, sharper first peak (Na+ piling up right at the
charged surface) than pH 7, that supports stronger counterion screening
at high surface charge density as the explanation for the smaller net
NP-NP drift observed at pH 10 despite higher nominal surface charge.

Atom types (from our LAMMPS data files):
  1=Si  2=Ob (bridging O)  3=Oh (silanol O)  4=Ho (silanol H)
  5=Ow  6=Hw  7=Na+  8=Od (deprotonated O)

Run this in Colab, pointing at the .lammpstrj trajectory files.
"""

import numpy as np
import matplotlib.pyplot as plt

# Surface atom types we measure distance TO (the charged/chargeable oxygens,
# since that's the physically relevant surface for counterion attraction)
SURFACE_TYPES = {3, 8}   # Oh (neutral silanol O) and Od (deprotonated O)
NA_TYPE = 7

R_MAX = 15.0     # Angstrom, max distance to consider
N_BINS = 100
FRAMES_TO_USE = 20   # use the last N frames of production for better statistics


def read_lammpstrj_frames(fname, n_frames=None):
    """Read all (or last n_frames) frames from a LAMMPS dump file.
    Returns list of (box_bounds, atom_id->(type,x,y,z)) tuples."""
    with open(fname) as f:
        lines = f.readlines()

    frame_starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    if n_frames is not None:
        frame_starts = frame_starts[-n_frames:]

    frames = []
    for start in frame_starts:
        n_atoms = int(lines[start+3])
        box_lines = lines[start+5:start+8]
        box_bounds = np.array([list(map(float, l.split())) for l in box_lines])
        atom_start = start + 9
        atom_data = {}
        for i in range(n_atoms):
            parts = lines[atom_start+i].split()
            aid, atype = int(parts[0]), int(parts[1])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            atom_data[aid] = (atype, x, y, z)
        frames.append((box_bounds, atom_data))
    return frames


def compute_rdf(frames, ref_types, target_type, r_max=R_MAX, n_bins=N_BINS):
    """Compute g(r)-like histogram of target_type atoms around ref_types atoms,
    using minimum-image convention for periodic distances."""
    bin_edges = np.linspace(0, r_max, n_bins+1)
    hist = np.zeros(n_bins)
    n_ref_total = 0

    for box_bounds, atom_data in frames:
        box_lo = box_bounds[:, 0]
        box_hi = box_bounds[:, 1]
        box_size = box_hi - box_lo

        ref_pos = np.array([(x,y,z) for (t,x,y,z) in atom_data.values() if t in ref_types])
        target_pos = np.array([(x,y,z) for (t,x,y,z) in atom_data.values() if t == target_type])

        if len(ref_pos) == 0 or len(target_pos) == 0:
            continue

        n_ref_total += len(ref_pos)

        for rp in ref_pos:
            delta = target_pos - rp
            delta -= box_size * np.round(delta / box_size)   # minimum image
            dist = np.linalg.norm(delta, axis=1)
            dist = dist[dist < r_max]
            h, _ = np.histogram(dist, bins=bin_edges)
            hist += h

    # normalize: divide by number of reference atoms and by shell volume (uniform density reference)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    shell_volumes = 4*np.pi*bin_centers**2 * (bin_edges[1]-bin_edges[0])
    n_frames_used = len(frames)

    # avoid divide by zero
    density_norm = hist / (n_ref_total * shell_volumes + 1e-12)
    return bin_centers, density_norm, hist


def main():
    # Usage: python 05_na_rdf_analysis.py   (expects traj_pH7.lammpstrj and traj_pH10.lammpstrj)
    print("Reading pH 7 trajectory (last frames)...")
    frames_ph7 = read_lammpstrj_frames('traj_pH7.lammpstrj', n_frames=FRAMES_TO_USE)
    print(f"  {len(frames_ph7)} frames loaded")

    print("Reading pH 10 trajectory (last frames)...")
    frames_ph10 = read_lammpstrj_frames('traj_pH10.lammpstrj', n_frames=FRAMES_TO_USE)
    print(f"  {len(frames_ph10)} frames loaded")

    r_ph7, g_ph7, raw_ph7 = compute_rdf(frames_ph7, SURFACE_TYPES, NA_TYPE)
    r_ph10, g_ph10, raw_ph10 = compute_rdf(frames_ph10, SURFACE_TYPES, NA_TYPE)

    plt.figure(figsize=(8,5))
    plt.plot(r_ph7, g_ph7, label='pH 7 (143 Na+/NP, 61.3% deprotonated)')
    plt.plot(r_ph10, g_ph10, label='pH 10 (234 Na+/NP, 99.9% deprotonated)')
    plt.xlabel('Distance from surface O (Å)')
    plt.ylabel('Na+ density (unnormalized g(r)-like)')
    plt.title('Na+ ion distribution around nanoparticle surface: pH 7 vs pH 10')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xlim(0, 12)
    plt.savefig('na_rdf_comparison_pH7_pH10.png', dpi=150, bbox_inches='tight')
    plt.show()

    # quantitative peak comparison
    peak_ph7_idx = np.argmax(g_ph7[:40])   # restrict search to near-surface region
    peak_ph10_idx = np.argmax(g_ph10[:40])
    print(f"\npH 7  first peak: r = {r_ph7[peak_ph7_idx]:.2f} Å, height = {g_ph7[peak_ph7_idx]:.4f}")
    print(f"pH 10 first peak: r = {r_ph10[peak_ph10_idx]:.2f} Å, height = {g_ph10[peak_ph10_idx]:.4f}")
    print(f"\nIf pH 10's peak is taller/sharper (closer to surface) than pH 7's,")
    print("this supports counterion condensation being stronger at pH 10,")
    print("consistent with the smaller net NP-NP repulsive drift observed there.")


if __name__ == '__main__':
    main()
