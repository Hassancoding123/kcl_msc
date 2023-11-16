# A system which has 3 SA-I peptides 
import MDAnalysis as mda
import numpy as np

# defining the universe of atoms using the gro and xtc files
u = mda.Universe("trimer_step4_1.gro", "trimer_step5_0.xtc")

# defining the 3 SA-I peptides - one SA-I has 850 atoms
peptide1 = u.atoms[0:850]
peptide2 = u.atoms[850:1700]
peptide3 = u.atoms[1700:2550]

# trajectory defined
trajectory = u.trajectory

# defining our norm function to find the vector norms
from numpy import linalg as LA
norm = LA.linalg.norm

# finding the size of the periodic box
from itertools import product
dimensions = u.dimensions[0]
pbc_scale = dimensions/2
