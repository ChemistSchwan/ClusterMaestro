Screeing
========

The Screening module only works with the ASE Turbomole-calculator class. For some calculations, the modified Turbomole class is
needed for the actual analysis (This modified Turbomole Calculator will be included in this package soon). These functions are
written to word with Turbomole 7.5. Other versions might also work, but this is not guaranteed.

A full analysis consists of:

1. Conformer analysis for monomer and dimer
2. Optiization for best monomer- and dimer-conformer
3. *more*, including potnetially eiger-analysis for MO-diagrams and aoforce frequency analysis
4. Complete binding analysis with decomposition of the structure parts
5. Analysis of various structures created by CREST
6. Analysis of the energy barrier of the rotation of a substituent
7. TDDFT calculation with escf

A script could look like the following:

.. code-block:: python
   :linenos:

   import ClusterMaestro as cm
   from ase.io import read, write

   params = {"title": "test_job",
                     "task": "optimize",
                     "basis set name": "cc-pVDZ",
                     "use dft" : True,
                     "density functional": "b3-lyp",
                     "dispersion correction": "d4",
                     "use redundant internals": True,
                     'use resolution of identity': True,
                     "ri memory": 4000,
                     'multiplicity': 1,
                     'total charge': 0,
                     }

   dftFiles = []

   # Conformer analysis and optimization of the created base-structure
   cm.analyze.screening.crest_opt(core="Ad", substituent="Me", threads=24, dftFiles=dftFiles,
           crest=True, startOpt=True, optScript="opt.sh", params=params)

   # Additional information. In this case: eiger and aoforce
   cm.analyze.screening.moreInfo(dftFiles=dftFiles, params=params)

   # Analysis of the binding energy. Including: decomposition of the structure, single point calculation and analysis.
   cm.analyze.screening.bindingInfo(startBind=True, dftFiles=dftFiles, spScript="sp.sh", params=params)
   cm.analyze.screening.bindingAnalysis()

   # Reading in the optimized structure.
   system = cm.System.Molecule(read("monomer/opt/coord", format="turbomole"))
   system.find_core()
   system.find_substituents()
   # Performing a rotation of the first substituent. constrains will be added automatically. Followed by analysis.
   cm.analyze.screening.createRotationStructs(system, "tm", params=params, rotrange=[0, 180, 5])
   cm.analyze.screening.analyseRotation()

   # Perform and analyse a TDDFT calculation from the monomer structure.
   cm.analyze.screening.createEcxitation(params=params, exMode="rpas", exString="a 20")

This script will analyse the Tetramethly Adamantane molecule completely. The individual parameters can be looked up
at the description of the different structures.

For analysis, the functions *sumAllInfo* can be used. Further evaluation, summary and presentation is left for the user.
I, for example, am using a latex-file that just loops over all the different folders in the summary-folder and created one
concise PDF-file to share and reference.


.. automodule:: ClusterMaestro.analyze.screening
   :members:
   :undoc-members:
   :show-inheritance:
