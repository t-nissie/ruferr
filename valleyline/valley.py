#!/usr/bin/env python3
import numpy as np
from typing import IO
from pathlib import Path

from ase import Atoms
from ase.filters import ExpCellFilter
from ase.optimize.bfgs import BFGS, BFGSMethod
from mace.calculators import mace_mp


class ConstrainedBFGS(BFGS):
    """
    BFGS optimizer implementing the fixed-u constraint (Valley Line Method).
    Fixes the norm of atomic displacements u = ||pos - ref_pos|| while 
    allowing internal coordinates and cell vectors to relax.
    """
    def __init__(self, atoms, reference_scaled_positions, **kwargs):
        super().__init__(atoms, **kwargs)
        # remve x,y,z translation here
        self.u_x = # calculate u_x, u_y, u_z here

    def step(self, gradient=None):
        gradient = self._get_gradient(gradient)
        optimizable = self.optimizable

        pos = optimizable.get_x()
        dpos, steplengths = self.prepare_step(pos, gradient)
        dpos = self.determine_step(dpos, steplengths)

        # Proposed next coordinates
        new_pos = pos + dpos

        # Enforce fixed-u constraint on atomic position degrees of freedom
        natoms = len(self.atoms)
        atomic_dpos = new_pos[:natoms * 3].reshape(natoms, 3)
        ref_flat = self.reference_scaled_positions.reshape(natoms, 3)

        # Center displacement vector (remove translation)
        disp = atomic_dpos - ref_flat
        disp -= disp.mean(axis=0)

        # Rescale displacement magnitude to exactly match target_u
        current_u = np.linalg.norm(disp)
        if current_u > 1e-12:
            disp = disp * (self.target_u / current_u)

        # Re-apply constrained atomic positions
        atomic_dpos = ref_flat + disp
        new_pos[:natoms * 3] = atomic_dpos.ravel()

        optimizable.set_x(new_pos)

        # Dump trajectory/restart state
        if hasattr(self.atoms, 'orig_cell'):
            self.dump((self.state.hessian, self.pos0, self.forces0, self.maxstep, self.atoms.orig_cell))
        else:
            self.dump((self.state.hessian, self.pos0, self.forces0, self.maxstep))


# 1. Structure setup
reference_scaled_positions=[
    [0.0, 0.0, 0.0],  # Ba
    [0.5, 0.5, 0.5],  # Ti
    [0.0, 0.5, 0.5],  # O1
    [0.5, 0.0, 0.5],  # O2
    [0.5, 0.5, 0.0],  # O3
]
scaled_positions = reference_scaled_positions.copy()
scaled_positions[0][2]=scaled_positions[0][2]+0.2
a = 4.01
atoms = Atoms(
    symbols=['Ba', 'Ti', 'O', 'O', 'O'],
    scaled_positions=scaled_positions,
    cell=[a, a, a + 0.1],
    pbc=True
)

atoms.calc = mace_mp(model="medium")

# 2. Display initial energy
initial_energy = atoms.get_potential_energy()
print(f"Initial potential energy: {initial_energy:.4f} eV")

# 3. Apply target amplitude u (e.g., u = 0.2 Å displacement norm)
cell_relax = ExpCellFilter(atoms, reference_scaled_positions)

target_u = 0.2  # Set desired displacement amplitude target
optimizer = ConstrainedBFGS(cell_relax, target_u=target_u)

optimizer.run(fmax=0.05)

# 4. Display optimized energy and results
final_energy = atoms.get_potential_energy()
print(f"Optimized potential energy (u={target_u:.2f} Å): {final_energy:.4f} eV")
print("Optimization finished successfully!")
print(atoms)
print(atoms.get_scaled_positions(wrap=True))
