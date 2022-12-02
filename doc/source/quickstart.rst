Quickstart
==========

This module is specifically written to handle molecular cluster materials, which consist of a core and a substituent structure.
In order to give a reasonable overview over the handling of atomic structrure, some code examples are given here. In order to
perform more complicated tasks, we recommend looking up the specific task in detail.

Initilizing structures
----------------------

Initilizing structures is very straight-forward. Just load the module, call the function and set the core-
and substituent structure and for dimer and solids some other parameters as well. Note, that the dimer as well as the
amorphSolid classes have defaults that might not fit your needs.:

.. code-block:: python
   :linenos:

   import ClusterMaestro as cm

   molecule = cm.initilize.init.monomer(core="SnS", substituent="Me")
   dimer = cm.initilize.init.dimer(core="GeS", substituent="Ph", distance=7)
   solid = cm.initilize.init.amorphSolid(core="Ad", substituent="Np", number=20, density=0.5)


Writing the structure into a simple structure file is also very easy. All the *System*-classes are based on
the *ase.Atoms* class and can use all its functions, including *write*:

.. code-block:: python
   :linenos:

   molecule.write("molecule.xyz", format="xyz")
   solid.write("solid.cif", format="cif")

Also, existing structures can be easily read-in using the *ase*-package and formatted as an object of the
desired ClusterMaestro-class by using core in the following fashion. This way, the core and substituent structures
are not being analyzed automatically. Depending on the structure, the analysation might not work as intended, manual
review of the resulting structures might be necessary. If however the structure was created using this ClusterMaestro
and the order of the atoms in the system was NOT changed, the function *matchPattern* can be used, which is much more
reliable in this case:

.. code-block:: python
   :linenos:

   import ClusterMaestro as cm
   from ase.io import read

   molecule = cm.System.Molecule(read("molecule.xyz"))
   molecule.analyze(core="SnS", substituent="Me")
   molecule.matchPattern(core="SnS", substituent="Me")

Writing these three code-blocks in a script and executing it will create the *molecule*, *dimer* and *solid*
object, write the *molecule*-object in an xyz-file and read it back in. This makes it easy to create structures,
prepare systems for more advanced tasks with it in other programs (e.g. DFT calculations) and read the result
back in as an object of the desired class.

.. note::
   It is HIGHLY recommended to create the structures with this very molecule, if more complex analyses should
   be done with it. This makes it very  as many functions use the *matchPattern* function to analyze and decompose the structure
   into molecules, cores and substituents. It makes it a lot easier to do so, if the order of the atoms
   is known.

Calling the various different functions is just as easy as in any other python-module. You need to know, what
you want to do and hope that this is implemented. As said before, this whole module is specificly written for
the use of the developer and might therefore not suit all of your needs.
