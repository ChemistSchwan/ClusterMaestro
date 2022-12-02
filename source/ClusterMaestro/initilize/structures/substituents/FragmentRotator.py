#!/usr/bin/env python
# -*- coding=utf-8 -*

"""
Rotates chosen Fragment to match the z-Axis. It can then be added MANUALLY to the FragemntList"
The Carbon, that belongs to the framestructure must be the frist atom in the file.
The atom attached to this carbon must be the second one.
"""
import numpy as np
import sys, os
from ase.io import read
from scipy.linalg import expm


def findVector(vector1, vector2):
    """
    Returns the vector from point "vector1" to point "vector2"
    """
    vector = np.zeros(3)
    for i in range(len(vector)):
        vector[i] = vector1[i] - vector2[i]
    return vector


def rotateVector(vector1, vector2, vector3):
    """
    Rotates vector 3 to match the direction of vector1. Vector2 is the starting orientation of the fragment.
    """
    c = np.dot(vector1, vector2) / np.linalg.norm(vector1) / np.linalg.norm(vector2)
    angle = np.arccos(np.clip(c, -1, 1))
    cross = np.cross(vector1, vector2)
    M0 = expm(np.cross(np.eye(3), cross / np.linalg.norm(cross) * angle))
    return np.dot(M0, vector3)


def main(fragmentFile):
    molecule = read(fragmentFile)
    molecule.positions -= molecule.positions[0]
    CX_vec = findVector(molecule.positions[1], molecule.positions[0])
    for i in range(len(molecule.positions)):
        molecule.positions[i] = rotateVector(
            CX_vec, np.array([0, 0, 1]), molecule.positions[i]
        )
    molecule.write(fragmentFile)


main(sys.argv[1])
# main("Ad.xyz")

with open(sys.argv[1], "r") as source:
    with open("tmp", "w") as target:
        counter = 0
        for line in source:
            if counter == 0:
                target.write("%d\n" % (int(line.split()[0]) - 1))
            elif counter != 2:
                target.write(line)
            counter += 1
os.rename("tmp", sys.argv[1])
