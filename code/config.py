#!/usr/bin/env python2
# -*- coding: utf-8 -*-

MIN_DIST = 0.5 # Minimum atoms distance
RES_LIGAND_MIN_DIST = 30.0 # Maximum RNA residue centroid and ligand centroid distance to calculate interactions
CUT_OFF_SIMPLE  = 4.0 # Distance cutoff for SIMPLE & PBS interactions
MAX_HB_DIST = 3.9 # Maximum Donor-Acceptor distance in Hydrogen Bond
MIN_HB_ANGLE = 100.0 # Minimum Donor-Hydrogen-Acceptor angle in Hydrogen Bond
MAX_HB_ANGLE = 260.0 # Maximum Donor-Hydrogen-Acceptor angle in Hydrogen Bond
MAX_HAL_DIST = 4.0 # Maximum Donor-Acceptor distance in Halogen Bond
HALOGEN_ACC_ANGLE = 165.0  # Preffered Donor-Halogen-Acceptor angle in Halogen Bond
HALOGEN_DON_ANGLE = 120.0 # Preffered Halogen-Acceptor-Acceptor' angle in Halogen Bond
HALOGEN_ANGLE_DEV = 30.0 # Tolerated halogen angles deviation
MAX_CA_DIST = 5.5 # Maximum cation-anion distance in electrostatic interaction
PI_ION_ANGLE = 90.0 # Preferred aromatic ring - cation/anion angle
PI_ION_ANGLE_DEV = 30.0 # Maximum angle deviation for Pi-cation/anion interaction
PI_ION_DISTANCE = 6.0 # Maximum cation/anion - aromatic ring center distance
RING_RING_MAX = 5.5 # Aromatic rings' centroids maximum distance
PISTACK_OFFSET_MAX = 2.0  # Maximum offset of the two aromatic rings (corresponds to the radius of benzene + 0.5 A)
PI_ANGLE_DISPLACED = 30.0 # Maximum angle value for Pi-stacking parallel & displaced interactions
PLANAR_ANGLE_TSHAPED = 90.0 # Preferred angle value for T-shaped Pi-stacking interaction
PLANAR_ANGLE_TSHAPED_DEV = 30.0 # Maximum angle deviation for the T-shaped Pi-stacking inetraction

GROUPS = [["P","OP1","OP2","OP3"],\
          ["C2","C4","C5","C6","C8","N1","N2","N3","N4","N6","N7","N9","O2","O4","O6","H1","H2","H21","H22","H3","H41","H42","H5","H6","H61","H62","H8","C7","H71","H72","H73"],\
          ["C1'","C2'","C3'","C4'","C5'","O2'","O3'","O4'","O5'","H1'","H2'","H3'","HO3'","H4'","H5'","H5''","HO2'","HO5'","H2''"]]
CANONICAL_RESIDUES = ['A','G','C','U','T']

OXYGEN_NUM = 8 # Oxygen atomic number
