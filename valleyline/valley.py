#!/usr/bin/env python3
import copy
import numpy as np
from typing import IO
from pathlib import Path

from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.optimize.bfgs import BFGS
from ase.optimize.optimize import Optimizer, UnitCellFilter
from mace.calculators import mace_mp


class valleyLineBFGS(BFGS):
    """
    BFGS optimizer implementing the fixed-u constraint (Valley Line Method).
    Fixes the norm of atomic displacements u = ||pos - ref_pos|| while 
    allowing internal coordinates and cell vectors to relax.
    """
    def __init__(self, filteredAtoms, scaled_centrosymmetrics, symmetry, **kwargs):
        self.scaled_centrosymmetrics=scaled_centrosymmetrics
        self.symmetry = symmetry
        displacements=filteredAtoms.atoms.remove_translations(self.scaled_centrosymmetrics, self.symmetry)
        self.norm=np.linalg.norm(displacements)
        self.unit_vector_n=displacements.ravel()/self.norm
        self.len=self.unit_vector_n.size
        super().__init__(filteredAtoms, **kwargs)
        print('displacements:\n', displacements)
        print(self.norm)
        print('unit_vector_n:\n', self.unit_vector_n.reshape([-1,3]))
        print(self.len)
        print(filteredAtoms.atoms.get_scaled_positions())

    def step(self, gradient=None):
        gradient = self._get_gradient(gradient)
        optimizable = self.optimizable
        pos = optimizable.get_x()
        dpos, steplengths = self.prepare_step(pos, gradient)
        #print(dpos.reshape([-1,3]))
        dpos[:self.len] -= (dpos[:self.len] @ self.unit_vector_n) * self.unit_vector_n
        print(dpos.reshape([-1,3]))
        steplengths = self.optimizable.gradient_norm(dpos)
        dpos = self.determine_step(dpos, steplengths)
        optimizable.set_x(pos+dpos)
        displacements=self.atoms.atoms.remove_translations(self.scaled_centrosymmetrics, self.symmetry, self.norm)
        norm=np.linalg.norm(displacements)
        self.unit_vector_n=displacements.ravel()/norm
        #print('displacements:\n', displacements)
        #print(norm)
        #print('unit_vector_n:\n', self.unit_vector_n.reshape([-1,3]))
        if isinstance(self.atoms, UnitCellFilter):
            self.dump((self.state.hessian, self.pos0, self.forces0,
                       self.maxstep, self.atoms.orig_cell))
        else:
            self.dump((self.state.hessian, self.pos0, self.forces0,
                       self.maxstep))

    # exit the valley-line loop other than fmax


def remove_translations(self, scaled_centrosymmetrics, symmetry, set_norm=-1000.0):
    """
    Remove x, y, z translations from displacements from a given scaled centrosymmetric structure.
    Return array of scaled displacements [v_x_scaled, v_y_scaled v_z_scaled].
    """
    sp = self.get_scaled_positions()
    if len(sp) != len(scaled_centrosymmetrics):
        raise ValueError("The number of atoms of Atoms and the given scaled centrosymmetric structure is different.")
    displacements = sp - scaled_centrosymmetrics
    displacements = (displacements + 0.5) % 1.0 - 0.5   # within [-0.5, 0.5)
    mean_translation = np.mean(displacements, axis=0)
    new_displacements = (displacements - mean_translation) * symmetry
    if set_norm>0.0:
        new_displacements = new_displacements * set_norm / np.linalg.norm(new_displacements)
    self.set_scaled_positions(scaled_centrosymmetrics+new_displacements)
    return new_displacements
# add a method to Atoms class
Atoms.remove_translations = remove_translations


#===Main===================================================================
# 1. Structure setup
scaled_centrosymmetrics=np.array([
    [0.0, 0.0, 0.0],  # Ba
    [0.5, 0.5, 0.5],  # Ti
    [0.0, 0.5, 0.5],  # O1
    [0.5, 0.0, 0.5],  # O2
    [0.5, 0.5, 0.0],  # O3
])
a = 4.01
cell=[a,a,a]
symmetry=[1,0,0] # [1,0,0] tetragonal, [1,1,0] orthorhombic, [1,1,1] rhombohedral

# make slightly displaced scaled positions.
scaled_positions = scaled_centrosymmetrics.copy()
amp=0.045
for index, value in enumerate(symmetry):
    scaled_positions[0][index]=scaled_positions[0][index]+amp*value
    scaled_positions[1][index]=scaled_positions[1][index]+amp*value
    cell[index]=cell[index]+amp*value
atoms = Atoms(
    symbols=['Ba', 'Ti', 'O', 'O', 'O'],
    scaled_positions=scaled_positions,
    cell=cell,
    pbc=True
)
atoms.calc = mace_mp()

# 2. Display initial energy
initial_energy = atoms.get_potential_energy()
print(f"Initial potential energy: {initial_energy:.4f} eV")

# 3. valley line optimization
if amp==0.000:
    optimizer =           BFGS(FrechetCellFilter(atoms))
    optimizer.run(fmax=0.02,steps=40)
    displacements = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]])
else:
    optimizer = valleyLineBFGS(FrechetCellFilter(atoms), scaled_centrosymmetrics, symmetry)
    optimizer.run(fmax=0.02,steps=40)
    displacements = atoms.remove_translations(           scaled_centrosymmetrics, symmetry)

# 4. report the result
norms=np.sqrt(np.sum(displacements**2, axis=0))
a,b,c,alpha,beta,gamma = atoms.cell.cellpar(radians=True)
v=norms*[a,b,c]
u=np.sqrt( v[0]**2 + v[1]**2 + v[2]**2 + 2*(  v[1]*v[2]*np.cos(alpha)
                                            + v[2]*v[0]*np.cos(beta )
                                            + v[0]*v[1]*np.cos(gamma)) )
eta_x = displacements[:,0]/np.linalg.norm(displacements[:,0])
print(gamma)
# valley: u E eta_A_x eta_B_x eta_O_I_x eta_O_II_x eta_O_III_x a b c alpha beta gamma
print('valley:', u, atoms.get_potential_energy(), *eta_x, *atoms.cell.cellpar())
