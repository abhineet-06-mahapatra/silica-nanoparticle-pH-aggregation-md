"""
Step 7-8: Generate LAMMPS input scripts (one per pH) for the
minimization -> Langevin equilibration -> NVT production protocol used
throughout this project.

Run with (2 OpenMP threads):
    OMP_NUM_THREADS=2 lmp -sf omp -pk omp 2 -in in.pH<X> -log run_pH<X>.log
"""

def write_lammps_input(pH_val, n_np_atoms_each, equil_steps=20000, prod_steps=50000):
    np1_start, np1_end = 1, n_np_atoms_each
    np2_start, np2_end = n_np_atoms_each + 1, 2 * n_np_atoms_each

    script = f"""
units           real
atom_style      full
boundary        p p p

pair_style      lj/cut/coul/long 10.0
bond_style      harmonic
angle_style     harmonic
kspace_style    pppm 1.0e-4

read_data       system_pH{pH_val}.data

pair_coeff 1 1 1.8405e-6 3.7064
pair_coeff 2 2 0.1554    3.5532
pair_coeff 3 3 0.1554    3.5532
pair_coeff 4 4 0.0       0.0
pair_coeff 5 5 0.1554    3.1660
pair_coeff 6 6 0.0       0.0
pair_coeff 7 7 0.1301    2.6378
pair_coeff 8 8 0.1554    3.5532

pair_modify     mix arithmetic

fix             rigidwater all shake 0.0001 20 0 b 1 a 1

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

group           np1_atoms id {np1_start}:{np1_end}
group           np2_atoms id {np2_start}:{np2_end}

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
fix             printdist all print 100 "$(step) ${{np_distance}}" file np_distance_pH{pH_val}.dat screen no

dump            traj all custom 1000 traj_pH{pH_val}.lammpstrj id type x y z
dump_modify     traj sort id

thermo          500
run             {prod_steps}

write_data      final_pH{pH_val}.data
"""
    fname = f'in.pH{pH_val}'
    with open(fname, 'w') as f:
        f.write(script)
    print(f"Saved: {fname}")


if __name__ == '__main__':
    import numpy as np
    PH_VALUES = [2, 4, 7, 8, 10, 12]
    for pH_val in PH_VALUES:
        # atoms per nanoparticle are read from the file written by 02_assign_pH_charges.py
        n_atoms = len(np.load(f'np_pH{pH_val}.npz')['symbols'])
        write_lammps_input(pH_val, n_atoms)
