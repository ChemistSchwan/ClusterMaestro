"""
The lib package ist a performance-enhanced library of ClusterMaestro. Some functionalities (especially concerning the
trajectory class are quite slow due to the shere quantity of information a classical MD or similar creates. To
comensate this, parts of the code were compiled using Cython.
"""

import ClusterMaestro.lib.solidFuncs

# import ClusterMaestro.lib.trajFuncs
