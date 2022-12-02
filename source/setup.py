#!/usr/bin/env python3

from setuptools import setup
from Cython.Build import cythonize

setup(name='ClusterMaestro',
      version='0.1.0',
      description='Python package created the handling of molecular clustermaterials',
      author='Sebastian Schwan',
      author_email='sschwan1994@gmail.com',
      package_dir={'ClusterMaestro': 'ClusterMaestro'},
      packages=setuptools.find_packages(),
      ext_modules = cythonize([ "ClusterMaestro/util/rmsdFAST.pyx",
                                 "ClusterMaestro/lib/solidFuncs.pyx",
                                 "ClusterMaestro/lib/shiftMols.pyx"]),
      include_package_data=True,
      install_requires=[
          "ase",
          "numpy",
          "scipy",
          "matplotlib",
          "pymatgen",
          "statistics",
          "MDAnalysis",
          "rmsd",
          ]
      )

