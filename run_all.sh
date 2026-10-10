#!/usr/bin/env bash
# Full pipeline. Requirements: python3 (ase, numpy, scipy, matplotlib), packmol on PATH, LAMMPS (lmp) with OPENMP.
# Run from a working directory that contains quartz_relaxed.cif. Use  SCRIPTS=/path/to/scripts  if needed.
set -e
S=${SCRIPTS:-$(dirname "$0")}
cp -n "$S/../data/quartz_relaxed.cif" . 2>/dev/null || true
python3 "$S/01_build_nanoparticle.py"
python3 "$S/02_assign_pH_charges.py"
python3 "$S/03_build_system_and_solvate.py"
python3 "$S/04_generate_lammps_inputs.py"
for pH in 2 4 7 8 10 12; do
  OMP_NUM_THREADS=2 lmp -sf omp -pk omp 2 -in in.pH$pH -log run_pH$pH.log
done
python3 "$S/07_plot_results.py"
python3 "$S/05_na_rdf_analysis.py"
# Salt study (pH 10 + NaCl / CaCl2): generates data files and inputs, then run each
python3 "$S/06_salt_concentration_study.py"
for l in pH10_NaCl_0.05M pH10_NaCl_0.3M pH10_CaCl2_0.05M pH10_CaCl2_0.3M; do
  OMP_NUM_THREADS=2 lmp -sf omp -pk omp 2 -in in.$l -log run_$l.log
done
python3 "$S/07_plot_results.py" salt
