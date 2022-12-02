Installation
============

.. attention::
   This module was developed and tested on linux-based systems exclusively. The usage under windows,
   MAC or BSD based systems has not been tested and might not work as intended.

Most important python-internal dependencies (will be installed automatically during the setup process with pip):

 - numpy
 - scipy
 - matplotlib
 - MDAnalysis
 - Cython
 - PyMatGen

External dependencies need to be installed/ downloaded manually be made executable on the local linux-machine, if are used that rely on these applications.
A note is written in the functions description, if this is the case:

 - Packmol (http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml)
 - xTB (https://github.com/grimme-lab/xtb)
 - Turbomole (https://www.turbomole.org/)

For developers:

 - sphinx and sphinx-rtd-theme (via pip) for the documentation
 - pytest and coverage for testing code

It is generally recommended to know how to use the python-package ASE (atomistic simulation environment), as big
parts of the code are using ASE-routines and all classes used to handle atomic structures of some sort are based on the ASE-Atoms class.

If you need the code/ a new version of the code, please feel free to contact me.

In order to actually install this module, download it, go to the folder with the setup.py-file and type:

.. code-block:: bash
   :linenos:

   pip3 install .

This should install the module and compile the cython-files for improved performance.
If this is not the case, please contact me and I will take a look at the problem/ patch the setup-file to work better.
