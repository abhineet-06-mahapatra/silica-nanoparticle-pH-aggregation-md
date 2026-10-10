"""
Plot NP-NP centre-of-mass separation vs time from the LAMMPS fix-print output
(np_distance_pH<X>.dat, columns: step distance) and report simple statistics.

Usage:  python 07_plot_results.py            # pH sweep: 2, 4, 7, 8, 10, 12
        python 07_plot_results.py salt       # salt study (4 runs)

Statistics are descriptive only: each condition is ONE short trajectory, so the
net drift is compared against the fluctuation (std) of the same trace.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

PH_VALUES = [2, 4, 7, 8, 10, 12]
SALT_LABELS = ['pH10_NaCl_0.05M', 'pH10_NaCl_0.3M', 'pH10_CaCl2_0.05M', 'pH10_CaCl2_0.3M']


def load(path):
    return np.loadtxt(path, comments='#')


def stats(step, d, n_blocks=5):
    blocks = np.array_split(d, n_blocks)
    means = [b.mean() for b in blocks]
    return dict(mean=d.mean(), std=d.std(), net_drift=means[-1] - means[0],
                block_means=means, range=d.max() - d.min())


def plot_ph():
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), sharey=False)
    print(f"{'pH':>4} {'mean(A)':>9} {'std(A)':>8} {'net drift(A)':>13}")
    for ax, pH in zip(axes.ravel(), PH_VALUES):
        try:
            a = load(f'np_distance_pH{pH}.dat')
        except OSError:
            ax.set_title(f'pH {pH}: file not found'); continue
        step, d = a[:, 0], a[:, 1]
        s = stats(step, d)
        ax.plot(step, d, lw=1)
        ax.set_title(f'pH {pH}'); ax.set_xlabel('timestep (1 fs)'); ax.set_ylabel('NP-NP distance (A)')
        ax.grid(alpha=.3)
        print(f"{pH:>4} {s['mean']:>9.3f} {s['std']:>8.3f} {s['net_drift']:>13.3f}")
    fig.tight_layout(); fig.savefig('np_distance_all_pH.png', dpi=200)


def plot_salt():
    plt.figure(figsize=(10, 5.5))
    for lab in SALT_LABELS:
        try:
            a = load(f'np_distance_{lab}.dat')
        except OSError:
            continue
        plt.plot(a[:, 0], a[:, 1], label=lab)
        s = stats(a[:, 0], a[:, 1])
        print(f"{lab:>20}: mean {s['mean']:.3f}  std {s['std']:.3f}  net drift {s['net_drift']:+.3f} A")
    plt.xlabel('timestep (1 fs)'); plt.ylabel('NP-NP distance (A)')
    plt.legend(); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig('salt_comparison.png', dpi=200)


if __name__ == '__main__':
    plot_salt() if len(sys.argv) > 1 and sys.argv[1] == 'salt' else plot_ph()
