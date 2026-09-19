#!/usr/bin/env python3
import copy
import numpy as np
from typing import IO
from pathlib import Path

from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize.bfgs import BFGS, BFGSMethod
from mace.calculators import mace_mp


# class ConstrainedBFGS(BFGS):
#     """
#     BFGS optimizer implementing the fixed-u constraint (Valley Line Method).
#     Fixes the norm of atomic displacements u = ||pos - ref_pos|| while 
#     allowing internal coordinates and cell vectors to relax.
#     """
#     def __init__(self, atoms, reference_scaled_positions, **kwargs):
#         super().__init__(atoms, **kwargs)
#         # remve x,y,z translation here
#         self.u_x = # calculate u_x, u_y, u_z here

#     def step(self, gradient=None):
#         gradient = self._get_gradient(gradient)
#         optimizable = self.optimizable

#         pos = optimizable.get_x()
#         dpos, steplengths = self.prepare_step(pos, gradient)
#         dpos = self.determine_step(dpos, steplengths)

#         # Proposed next coordinates
#         new_pos = pos + dpos

#         # Enforce fixed-u constraint on atomic position degrees of freedom
#         natoms = len(self.atoms)
#         atomic_dpos = new_pos[:natoms * 3].reshape(natoms, 3)
#         ref_flat = self.reference_scaled_positions.reshape(natoms, 3)

#         # Center displacement vector (remove translation)
#         disp = atomic_dpos - ref_flat
#         disp -= disp.mean(axis=0)

#         # Rescale displacement magnitude to exactly match target_u
#         current_u = np.linalg.norm(disp)
#         if current_u > 1e-12:
#             disp = disp * (self.target_u / current_u)

#         # Re-apply constrained atomic positions
#         atomic_dpos = ref_flat + disp
#         new_pos[:natoms * 3] = atomic_dpos.ravel()

#         optimizable.set_x(new_pos)

#         # Dump trajectory/restart state
#         if hasattr(self.atoms, 'orig_cell'):
#             self.dump((self.state.hessian, self.pos0, self.forces0, self.maxstep, self.atoms.orig_cell))
#         else:
#             self.dump((self.state.hessian, self.pos0, self.forces0, self.maxstep))

def remove_translations(self, centrosymmetric):
    """
    Remove x, y, z translations from displacements from a given scaled centrosymmetric structure.
    Return a vector of [u_x^2, u_y^2, u_z^2]
    """
    sp = self.get_scaled_positions()
    if len(sp) != len(centrosymmetric):
        raise ValueError("The number of atoms of Atoms and the given centrosymmetric structure is different.")
    displacement = sp - centrosymmetric
    displacement = (displacement + 0.5) % 1.0 - 0.5   # within [-0.5, 0.5)
    mean_translation = np.mean(displacement, axis=0)
    new_sp = sp - mean_translation
    self.set_scaled_positions(new_sp)
    new_displacement = new_sp - centrosymmetric
    new_displacement = (new_displacement + 0.5) % 1.0 - 0.5   # within [-0.5, 0.5)
    u2=np.sum(new_displacement**2, axis=0)
    return u2

# add a method to Atoms class
Atoms.remove_translations = remove_translations


#===Main===================================================================
# 1. Structure setup
reference_scaled_positions=np.array([
    [0.0, 0.0, 0.0],  # Ba
    [0.5, 0.5, 0.5],  # Ti
    [0.0, 0.5, 0.5],  # O1
    [0.5, 0.0, 0.5],  # O2
    [0.5, 0.5, 0.0],  # O3
])
scaled_positions = reference_scaled_positions.copy()
scaled_positions[0][2]=scaled_positions[0][2]+0.1
scaled_positions[1][2]=scaled_positions[1][2]+0.1
a = 4.01
atoms = Atoms(
    symbols=['Ba', 'Ti', 'O', 'O', 'O'],
    scaled_positions=scaled_positions,
    cell=[a, a, a + 0.1], # tetragonal
    pbc=True
)
print(atoms.remove_translations(reference_scaled_positions))
print(atoms.get_scaled_positions())



atoms.calc = mace_mp(model="medium")

# 2. Display initial energy
initial_energy = atoms.get_potential_energy()
print(f"Initial potential energy: {initial_energy:.4f} eV")

# 3. Apply target amplitude u (e.g., u = 0.2 Å displacement norm)
cell_relax = FrechetCellFilter(atoms)

target_u = 0.2  # Set desired displacement amplitude target
#optimizer = ConstrainedBFGS(cell_relax, reference_scaled_positions)
optimizer = BFGS(cell_relax)

optimizer.run(fmax=0.05)

# 4. Display optimized energy and results
final_energy = atoms.get_potential_energy()
print(f"Optimized potential energy (u={target_u:.2f} Å): {final_energy:.4f} eV")
print("Optimization finished successfully!")
print(atoms)
print(atoms.get_scaled_positions(wrap=True))
