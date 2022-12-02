""" Small extension of the ASE.Atom class, which is currently not very useful.
"""

import ase
import numpy as np


class Atom(ase.Atoms):
    """Atom object.

    Inherits from the ase-Atom object.
    See "https://wiki.fysik.dtu.dk/ase/ase/atom.html" for more information
    """

    def copy(self):
        """Creates a deepcopy of the Atom."""

        atom = copy.deepcopy(self)
        return atom
