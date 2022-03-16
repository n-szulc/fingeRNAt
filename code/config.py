### fingeRNAt configuration file ###

MIN_DIST = 0.5 # Minimum atoms distance
RES_LIGAND_MAX_DIST = 85.0 # Maximum nucleic acid residue's centroid - ligand's centroid distance to calculate interactions; decrease to speed up software but some interactions may be omitted!
CUT_OFF_SIMPLE  = 4.0 # Distance cutoff for SIMPLE & PBS interactions
MAX_HB_DIST = 3.9 # Maximum Donor-Acceptor distance in Hydrogen Bond (Torshin, Weber, & Harrison, 2002)
MIN_HB_ANGLE = 100.0 # Minimum Donor-Hydrogen-Acceptor angle in Hydrogen Bond (Adasme et al., 2021)
MAX_HB_ANGLE = 260.0 # Maximum Donor-Hydrogen-Acceptor angle in Hydrogen Bond (Adasme et al., 2021)
MAX_HAL_DIST = 4.0 # Maximum Donor-Acceptor distance in Halogen Bond (Auffinger et al., 2004)
HALOGEN_ACC_ANGLE = 165.0  # Preffered Donor-Halogen-Acceptor angle in Halogen Bond (Auffinger et al., 2004)
HALOGEN_DON_ANGLE = 120.0 # Preffered Halogen-Acceptor-Acceptor' angle in Halogen Bond (Auffinger et al., 2004)
HALOGEN_ANGLE_DEV = 30.0 # Tolerated halogen angles deviation
MAX_CA_DIST = 5.5 # Maximum cation-anion distance in electrostatic interaction (Barlow and Thornton, 1983)
PI_ION_ANGLE = 90.0 # Preferred aromatic ring - cation/anion angle
PI_ION_ANGLE_DEV = 30.0 # Maximum angle deviation for Pi-cation/anion interaction
PI_ION_DISTANCE = 6.0 # Maximum cation/anion - aromatic ring center distance (Gallivan and Dougherty, 1999)
RING_RING_MAX = 5.5 # Aromatic rings' centroids maximum distance (McGaughey, 1998)
PISTACK_OFFSET_MAX = 2.0  # Maximum offset of the two aromatic rings (corresponds to the radius of benzene + 0.5 A)
PI_ANGLE_DISPLACED = 30.0 # Maximum angle value for Pi-stacking parallel & displaced interactions
PLANAR_ANGLE_TSHAPED = 90.0 # Preferred angle value for T-shaped Pi-stacking interaction
PLANAR_ANGLE_TSHAPED_DEV = 30.0 # Maximum angle deviation for the T-shaped Pi-stacking interaction
MAX_ION_DIST = 3.9 # Maximum ligand's nitrogen/oxygen/sulphur atom - positively charged ion or residue's nitrogen/oxygen atom - positively charged ion distance
MAX_MAGNESIUM_DIST = 3.2 # Maximum ligand's nitrogen/oxygen/sulphur atom - magnesium or residue's nitrogen/oxygen atom - magnesium distance (Zheng et al., 2015)
MAX_POTASSIUM_DIST = 3.9 # Maximum ligand's nitrogen/oxygen/sulphur atom - potassium or residue's nitrogen/oxygen atom - potassium distance (Zheng et al., 2008)
MAX_SODIUM_DIST = 3.6 # Maximum ligand's nitrogen/oxygen/sulphur atom - sodium or residue's nitrogen/oxygen atom - sodium distance (Zheng et al., 2008)
MAX_OTHER_ION_DIST = 3.5 # Maximum ligand's nitrogen/oxygen/sulphur atom - other than magnesium/potassium/sodium ion or residue's nitrogen/oxygen atom - other than magnesium/potassium/sodium ion distance (Zheng et al., 2008)
MAX_WATER_DIST = 3.5 # Maximum ligand's hydrogen donors/acceptors - water molecule (oxygen) or nucleic acid's hydrogen donors/acceptors - water molecule (oxygen) distance (Poornima & Dean, 1995)
MAX_LIPOHILIC_DIST = 4.0 # Maximum ligand - residue lipophilic contact distance (Padroni et al., 2020)

GROUPS = [["P","OP1","OP2","OP3"],\
          ["C2","C4","C5","C6","C8","N1","N2","N3","N4","N6","N7","N9","O2","O4","O6","C7"],\
          ["C1'","C2'","C3'","C4'","C5'","O2'","O3'","O4'","O5'"]]
WHICH_GROUP = {0:'PHOSPHATE', 1:'BASE', 2:'SUGAR'}
CANONICAL_RESIDUES = ['A','G','C','U','T']
POS_CHARGED_IONS = ["MG", "CA", "ZN", "NA", "K", "MN", "FE", "CO", "NI", "LI", "CU", "AL",
                    "PT", "HG", "CD", "YB", "SR", "GD", "HO", "AU", "RB", "CS", "EU", "SM",
                    "PB", "CE", "BA"]

OXYGEN_NUM = 8 # Oxygen atomic number
NITROGEN_NUM = 7 # Nitrogen atomic number
SULPHUR_NUM = 16 # Sulphur atomic number
HYDROGEN_NUM = 1 # Hydrogen atomic number
PHOSPHORUS_NUM = 15 # Phosphorus atomic number
CARBON_NUM = 6 # Carbon atomic number
