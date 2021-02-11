import numpy as np

from openbabel import openbabel
from openbabel import pybel
from rdkit import Chem

from tqdm import tqdm
import collections

# Own module
import config


def vector(p1, p2):
    """Calculates vector between 2 points

    :param p1: coordinates of point 1
    :param p2: coordinates of point 2
    :type p1: tuple, numpy.ndarray
    :type p2: tuple, numpy.ndarray
    :return: vector between two points
    :rtype: numpy.ndarray
    """

    return None if len(p1) != len(p2) else np.array([p2[i] - p1[i] for i in range(len(p1))])

def measure_distance(v1, v2):
    """Calculates distance between 2 vectors

    :param v1: vector 1 coordinates
    :param v2: vector 2 coordinates
    :type v1: tuple, numpy.ndarray
    :type v2: tuple, numpy.ndarray
    :return: distance between 2 vectors
    :rtype: numpy.float64
    """

    return np.round(np.linalg.norm(v1 - v2), 5)

def check_distance(v1, v2, cutoff):
    """Checks if distance between 2 vectors is within cutoff range and return bool

    :param v1: vector 1 coordinates
    :param v2: vector 2 coordinate
    :param cutoff: cutoff value
    :type v1: numpy.ndarray
    :type v2: numpy.ndarray
    :type cutoff: float
    :return: `True` if distance between 2 vectors is within cutoff range, `False` otherwise
    :rtype: bool
    """

    return (np.all(np.linalg.norm(v1 - v2)  <= cutoff))

def calculate_angle(v1, v2, deg = True):
    """Calculates angle between 2 vectors

    :param v1: vector 1 coordinates
    :param v2: vector 2 coordinates
    :param deg: output type; default degrees
    :type v1: numpy.ndarray
    :type v2: numpy.ndarray
    :type deg: bool
    :return: calculated angle value; in degrees or radians
    :rtype: numpy.float64
    """

    if np.array_equal(v1, v2): return 0.0
    dm = np.dot(v1, v2)
    cm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(round(dm/cm, 10))

    return np.degrees([angle, ])[0] if deg else angle

def modify_HB_result_list(precision, result, dist):
    """Modifies hydrogen bond output list according to it's presence and strength

    :param precision: fingerprint type
    :param result: list to modify
    :param dist: hydrogen bond donor - hydrogen bond acceptor distance
    :type precision: str
    :type result: list
    :type dist: numpy.float64
    :return: modified hydrogen bond output list
    :rtype: NoneType
    """

    if precision == 'FULL':
        result[-1] = 1

    else:
        result[-4]+= 1

        if dist <=  2.5 and dist >= 2.2: # Strong hydrogen bond
            result[-3] += 1
        elif dist <= 3.5 and dist > 2.5: # Moderate hydrogen bond
            result[-2] += 1
        elif dist > 3.5: # Weak hydrogen bond
            result[-1] += 1
        else:
            print('Too short D-A distance: %s A' %str(dist))

def calculate_planar(atoms_list):
    """Calculates planarity normal vector

    :param atoms_list: list of 3 atoms creating planar space
    :type atoms_list: numpy.ndarray
    :return: planarity normal vector
    :rtype: numpy.ndarray
    """

    v1 = np.array(atoms_list[-1]) - np.array(atoms_list[-3])
    v2 = np.array(atoms_list[-2]) - np.array(atoms_list[-3])

    cp = np.cross(v1, v2)
    a, b,c  = cp
    pl = [a, b, c]

    return np.round(pl, 5)

def centroid(coo):
    """Finds ring center

    :param coo: list of ring's atoms
    :type coo: list
    :return: coordinates of ring center
    :rtype: numpy.ndarray
    """

    return np.round(np.array(list(map(np.mean, (([c[0] for c in coo]), ([c[1] for c in coo]), ([c[2] for c in coo]))))), 5)

def get_ligand_name_pose(dictionary, title):
    """Gets currently processed ligand's name and it's pose number from ligands input file

    :param dictionary: dictionary of all so far processed ligands
    :param title: currently processed ligand's title obtained from parsing by OpenBabel
    :type: dictionary: dict
    :type: title: str
    :return: ligand's name as ligand_name^pose_number
    :rtype: str
    """

    dictionary_keys_only_ligand_names = [el.split('^')[0] for el in dictionary.keys()] # list of ligands names (prefixes; without pose no) that are so far present in the dictionary
    title_cleaned = title.split('/')[-1]

    if title_cleaned == "":
        title_cleaned = "COMPOUND"

    if title_cleaned not in dictionary_keys_only_ligand_names: # if molecule title (prefix) is not present in dict keys prefixes
        name = title_cleaned + '^1' # prefix^pose
    else:
        previous_poses_numbers = [int(el.split('^')[-1]) for el in dictionary.keys() if el.split('^')[0] == title_cleaned] # find all dict keys with the same ligand name (prefix) and save only pose no to list
        previous_pose_number = max(previous_poses_numbers) # find biggest value of ligand pose no
        current_pose_number = str(int(previous_pose_number)+1) # to get next pose no add 1
        name = title_cleaned + '^' + current_pose_number

    return name

def find_ligands_all_atoms(mols):
    """Creates dictionary of ligand_name^pose_number with list of sublists of ligand's all atoms coords (except hydrogens)

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary of ligand_name^pose_number : [list of sublists of ligand's all atoms coords]
    :rtype: dict
    """

    dictionary = {}

    for i in range(len(mols)):
        name = get_ligand_name_pose(dictionary, mols[i].title)
        tmp=[]

        for atom in mols[i].atoms:
            if atom.atomicnum != 1:
                tmp.extend(np.array([atom.coords]))
        dictionary[name] = tmp # {'prefix^pose':[list of sublists of all atoms coords]}

    return dictionary

def projection(pnormal1, ppoint, tpoint):
    """Calculates the centroid from a 3D point cloud and returns the coordinates

    :param pnormal1: normal of plane
    :param ppoint: coordinates of point in the plane
    :param tpoint: coordinates of point to be projected
    :type pnormal1: numpy.ndarray
    :type ppoint: numpy.ndarray
    :type tpoint: numpy.ndarray
    :return: coordinates of point orthogonally projected on the plane
    :rtype: list
    """

    # Choose the plane normal pointing to the point to be projected
    pnormal2 = np.array([coo*(-1) for coo in pnormal1])
    d1 = measure_distance(tpoint, pnormal1 + ppoint)
    d2 = measure_distance(tpoint, pnormal2 + ppoint)
    pnormal = pnormal1 if d1 < d2 else pnormal2

    # Calculate the projection of tpoint to the plane
    sn = -np.dot(pnormal, vector(ppoint, tpoint))
    sd = np.dot(pnormal, pnormal)
    sb = sn / sd

    return [c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])]

def assign_interactions_results(result, RESULTS, RNA_LENGTH, RNA_seq_index, FINGERPRINT_DESCRIPTORS_NO, desc_index, multiple_results = False):
    """Adds interaction from residue - ligand pair into fingerprint held in a dictionary indexed by ligand ID

    :param result: result obtained from calling one of functions that calculate molecular interaction between residue-ligand
    :param RESULTS: final dictionary holding all SIFs for each ligand - nucleic acid complex
    :param RNA_LENGTH: total length of nucleic acid chains from the input file
    :param RNA_seq_index: index of position of particular residue in the dictionary value list, does not have to correspond to residue number due to possible shifts in sequence
    :param FINGERPRINT_DESCRIPTORS_NO: value responding to total no of calculated interactions for each fingerprint type
    :param descriptor_index: index of particular descriptor in the dictionary value list
    :param multiple_results: indicates whether we deal with result with more than one value to add to the dictionary value list
    :type result: list
    :type RESULTS: dict
    :type RNA_LENGTH: int
    :type RNA_seq_index: int
    :type FINGERPRINT_DESCRIPTORS_NO: int
    :type descriptor_index: int
    :type multiple_results: bool
    :return: modified RESULTS dictionary
    :rtype: None
    """

    ligand_id = result[0]
    descriptor_index = desc_index

    if not multiple_results: # Most results have only one value to add

        # XP fingerprint 3 additional HB strong/moderate/weak descriptors have to be taken into account
        if FINGERPRINT_DESCRIPTORS_NO == 9:
            # Change index in XP fingerprint for halogen bondings, cation-anion, Pi-interactions (XP hydrogen bonding is the first calculated interaction and it took 3 additional places in results list)
            descriptor_index += 3

        RESULTS[ligand_id][(RNA_seq_index * FINGERPRINT_DESCRIPTORS_NO) + descriptor_index] = result[-1]

    else: # In case of PBS and HB_XP results

        if FINGERPRINT_DESCRIPTORS_NO == 9: # Hydrogen bondings XP result
            for i in range(4):
                RESULTS[ligand_id][(RNA_seq_index*FINGERPRINT_DESCRIPTORS_NO) + descriptor_index+i] = result[-4+i]

        else: # PBS result
            for i in range(3):
                RESULTS[ligand_id][(RNA_seq_index*FINGERPRINT_DESCRIPTORS_NO) + descriptor_index+i] = result[-3+i]

def wrap_results(wrapper, RESULTS, RNA_nucleotides, fingerprint_length, wrapper_length):
    """Wrap results from calculating fingerprint to the desired output

    :param wrapper: wrapper type
    :rtype wrapper: str
    :param RESULTS: calculated fingerprint results
    :rtype RESULTS: dict
    :param RNA_nucleotides: list of nucleic acid nucleotides names
    :rtype RNA_nucleotides: list
    :param fingerprint_length: number of calculated molecular interactions
    :rtype fingerprint: int
    :param wrapper_length: length of wrapper
    :rtype wrapper_length: int
    :return: list of wrapped results
    :rtype: list

    """
    WRAP_RESULTS = {}

    if wrapper == 'ACUG':
        letter_order = {'A':0, 'C':1, 'U':2, 'G':3}
    if wrapper == 'PuPy':
        letter_order = {'A':0, 'C':1, 'U':1, 'G':0}

    for key, values in RESULTS.items():

        WRAP_RESULTS[key] = [0] * fingerprint_length * wrapper_length

        i=0

        while i < len(values):
            chunk = values[i:(i+fingerprint_length)]

            if wrapper == 'ACUG' or wrapper == 'PuPy':

                try:
                    nucleotide_index = letter_order[RNA_nucleotides[(int(i/fingerprint_length))]]

                except KeyError: # non-canonical nucleotide
                    i += fingerprint_length
                    continue

                for el in range(len(chunk)):
                    if fingerprint_length != 9: # if not XP, we overwrite 0 with 1
                        if WRAP_RESULTS[key][fingerprint_length*nucleotide_index+el] < chunk[el]:
                             WRAP_RESULTS[key][fingerprint_length*nucleotide_index+el] = chunk[el]
                    else:
                        WRAP_RESULTS[key][fingerprint_length*nucleotide_index+el] += chunk[el] # XP, we sum all interactions

            else: # wrapper Counter
                for el in range(len(chunk)):
                    WRAP_RESULTS[key][el] += chunk[el]

            i += fingerprint_length

    return WRAP_RESULTS

def find_ligands_HBA_HBD(mols, verbose):
    """Finds all donors/acceptors in all ligands

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary indexed by ligand name, with the coords od all ligand's hydrogen bonds acceptors & donors
    :rtype: dict
    """

    dictionary = {}

    if verbose:
        print("Looking for hydrogen bonds donors & acceptors...")

    for i in tqdm(range(len(mols)), disable=(not verbose)):

        name = get_ligand_name_pose(dictionary, mols[i].title)
        dictionary[name] = [] # {'prefix^pose':[[sublist of acceptors coords],[sublist of tuples donor-H coords]]}
        acceptors = []
        donors = []

        for atom in filter(lambda at: at.OBAtom.IsHbondAcceptor(), mols[i].atoms): # Find all acceptors
            if atom.atomicnum not in [9, 17, 35, 53]: # Exclude halogen atoms
                acceptors.append(atom.coords)

        for atom in filter(lambda at: at.OBAtom.IsHbondDonor(), mols[i].atoms): # Find all donors with their hydrogens
            for neighbour in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                if neighbour.IsHbondDonorH():
                    neighbour_coords = np.array([neighbour.GetX(),neighbour.GetY(),neighbour.GetZ()])
                    donors.append((atom.coords, neighbour_coords))

        dictionary[name].append(acceptors)
        dictionary[name].append(donors)

    return dictionary

def find_ligands_HAL_don(mols, verbose):
    """Finds all halogens donors in all ligands

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary indexed by ligand name, with the coords od all ligand's halogens & halogen bonds donors
    :rtype: dict
    """

    dictionary = {}

    if verbose:
        print("Looking for halogen bonds donors & acceptors...")

    for i in tqdm(range(len(mols)), disable=(not verbose)): # for molecule in ligand file

        name = get_ligand_name_pose(dictionary, mols[i].title)
        dictionary[name] = [] # {'prefix^pose':[list of tuples (C,halogen)]}
        halogen_donors = []

        for atom in mols[i].atoms:
            if atom.atomicnum in [9, 17, 35, 53]:  # Include only halogen atoms
                for neighbour in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                    neighbour_coords = np.array([neighbour.GetX(),neighbour.GetY(),neighbour.GetZ()]) # Get Carbon coords
                    halogen_donors.append((neighbour_coords, atom.coords))

        dictionary[name].extend(halogen_donors)

    return dictionary

def find_ligands_CA(mols, verbose):
     """Finds all cations & anions in all ligands

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary indexed by ligand name, with the coords of all ligand's cations & anions
    :rtype: dict
    """

     dictionary = {}

     if verbose:
         print("Looking for cation-anion interactions...")

     for i in tqdm(range(len(mols)), disable=(not verbose)): # For molecule in ligand file

        name = get_ligand_name_pose(dictionary, mols[i].title)
        dictionary[name]=[] # {'prefix^pose':[[list of cations coords],[list of anions coords]]}
        cations = []
        anions = []

        for atom in mols[i]:
            if atom.formalcharge >= 1:
                cations.append(atom.coords)
            if atom.formalcharge <= -1:
                anions.append(atom.coords)

        dictionary[name].append(cations)
        dictionary[name].append(anions)

     return dictionary

def find_ligands_lipophilic(mols, verbose):
    """Finds lipophilic fragments in all ligands

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary indexed by ligand name, with the coords od all ligand's lipophilic fragments
    :rtype: dict
    """

    # SMARTS pattern:
    # [CH0,CH1,CH2,#9,#17,#35,#53] - aliphatic C with 0,1 or 2 H (ie not CH3) or halogens
    # ;+0   and only neutral (charge zero)
    # ;!$(C~O);!$(C~N)  and not C=O, C=N with any bonds
    # ;!$(*~[+1]);!$(*~[-1]) and not connected to a cation or anion
    # ICM: [C&!$(C=O)&!$(C#N),S&^3,#17,#15,#35,#53]
    # smarts = pybel.Smarts("[CH0,CH1,CH2,CH3,#9,#17,#35,#53,SX2;+0;!$(C~O);!$(C~N);!$(*~[+1]);!$(*~[-1])]")
    smarts = pybel.Smarts("[C,#9,#17,#35,#53,SX2;+0;!$(C~O);!$(C~N);!$(*~[+1]);!$(*~[-1])]")

    dictionary = {}

    if verbose:
        print("Looking for lipophilic fragments")

    for i in tqdm(range(len(mols)), disable=(not verbose)): # for molecule in ligand file

        name = get_ligand_name_pose(dictionary, mols[i].title)
        dictionary[name] = [] # {'prefix^pose':[list of tuples (C,halogen)]}

        atomSets = smarts.findall(mols[i]) # list of atoms fulfilling this pattern
        atomsList = [ id[0] for id in atomSets ]


        for atom in atomsList:
            dictionary[name].append( mols[i].atoms[atom-1].coords )

    return dictionary



def findAromaticRingsWithRDKit(sdfFile):
    """Finds all aromatic rings coords in all ligands using RDKit
    :param sdfFile: path to a valid sdf file
    :type sdfFile: str
    :return: dictionary indexed by ligand name, with the coords of all ligand's aromatic rings
    :rtype: dict
    """

    suppl = Chem.SDMolSupplier(sdfFile, sanitize=False) # do not sanitize at the beginning, do it later

    aromaticRingsDict = {}

    # iterate over molecules
    for m in suppl:
        aromaticRingsAtomsCoords = []

        # Ignore Explicit Valence Error - do only a partial sanitization here
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)

        if m is not None:
            # get rings
            ssr = Chem.GetSymmSSSR(m)
            molTitleTemp = m.GetProp('_Name')
            molTitle = get_ligand_name_pose(aromaticRingsDict, molTitleTemp)

            for ring in ssr:
                is_aromatic = True
                for atom_idx in ring:
                    # if any atom of the ring is not aromatic...
                    if not m.GetAtoms()[atom_idx].GetIsAromatic():
                        is_aromatic = False
                        break

                if is_aromatic == True:
                    tempRingCoords = []
                    for ring_idx, atom_idx in enumerate(ring):
                        # only aromatic rings here
                        position = m.GetConformer().GetAtomPosition(atom_idx)
                        tempRingCoords.append((position.x, position.y, position.z))

                    aromaticRingsAtomsCoords.extend([tempRingCoords])

            aromaticRingsDict[molTitle] = aromaticRingsAtomsCoords

    return aromaticRingsDict

def find_RNA_rings(structure, extension_structure):
    """Finds all aromatic rings in whole nucleic acid

    :param structure: nucleic acid structure object
    :param extension_structure: extension of nucleic acid input file
    :type structure: pybel.Molecule
    :type extension_structure: str
    :return: list of coords of all aromatic rings in nucleic acid found by OpenBabel
    :rtype: list
    """

    all_atoms = structure.atoms
    rings_candidates = structure.OBMol.GetSSSR()
    rings = []

    for ring in rings_candidates:
            ring_atoms = [a for a in all_atoms if ring.IsMember(a.OBAtom)]
            res = structure.OBMol.GetAtom(ring_atoms[0].idx).GetResidue() # Residue according to first ring's atom

            if extension_structure == 'pdb':
                res_id = res.GetName()
            elif extension_structure == 'mol2':
                res_id = res.GetName()[:-len(str(res.GetNum()))] # Need res.GetName() without last characters of its res number, because when nucleic acid is in mol2 format, res.GetName() gives for example U22 as residue name
            else:
                raise Exception('Unknown nucleic acid structure format')

            if res_id in config.CANONICAL_RESIDUES:
                sugar = False
                if len([atom for atom in ring_atoms if (atom.atomicnum == config.OXYGEN_NUM)]) > 0:
                    sugar = True
                if not sugar:
                    rings.append(ring)
            else:
                if ring.IsAromatic():
                    rings.append(ring)

    return rings

def find_RNA_HB_HAL_acc_don(residue):
    """Finds all hydrogen/halogen bonds acceptors with all of their neighbours and all hydrogen bonds donors together with hydrogens in nucleic acid residue

    :param residue: nucleic acid residue object
    :type residue: openbabel.OBResidue
    :return: coords of residue's halogen/hydrogen bonds acceptors & their neighbours
    :return: coords of residue's hydrogen bonds donors together with hydrogens
    :rtype: list
    :rtype: list
    """

    residue_atoms_list = []
    for atom in openbabel.OBResidueAtomIter(residue):
        residue_atoms_list.append(atom) # Create list of all residue's atoms

    acceptors_RNA = [] # List of sublists [[acceptor,acceptor',acceptor'(if the 2nd one exists)]]
    donors_RNA = []

    for atom in filter(lambda at: at.IsHbondAcceptor(), residue_atoms_list): # Find all acceptors
        acceptors_RNA.append([atom])
        for neighbour in pybel.ob.OBAtomAtomIter(atom):
             if neighbour.GetAtomicNum() != 1: # If not hydrogen
                 acceptors_RNA[-1].append(neighbour) # Append acceptors' - all Y

    for atom in filter(lambda at: at.IsHbondDonor(), residue_atoms_list): # Find all donors with their hydrogens
        for neighbour in pybel.ob.OBAtomAtomIter(atom):
            if neighbour.IsHbondDonorH():
                donors_RNA.append((atom, neighbour))

    return acceptors_RNA, donors_RNA

def find_RNA_anions(residue):
    """Finds all nucleic acid's residue's anions

    :param residue: nucleic acid's residue object
    :type residue: openbabel.OBResidue
    :return: coords of residue's anions
    :rtype: list
    """

    P = None
    anions_RNA = []

    for atom in openbabel.OBResidueAtomIter(residue):
        if atom.GetAtomicNum() == config.PHOSPHORUS_NUM:
            P = atom
            break # Assuming there is only 1 Phosphate in residue

    if P: # 5'end has no phosphate group
        O = [] # List of 2 negative charged Oxygens
        for oxygen in pybel.ob.OBAtomAtomIter(P):
            # https://open-babel.readthedocs.io/en/latest/UseTheLibrary/migration.html
            if oxygen.GetExplicitValence() == 1 or oxygen.HasDoubleBond():
                O.append(oxygen)
        anions_RNA.extend(O)

    return anions_RNA

def check_if_RNA(all_nucleotides):
    """Check if input structure is RNA or DNA

    :param all_nucleotides: list of all structure nucleotides
    :type all_nucleotides: list
    :return: information if structure is RNA (True)
    :rtype: bool
    """

    if 'U' in all_nucleotides:
        if 'T' in all_nucleotides:
            raise Exception('Invalid input structure, has both uracil and thymine')
        return True
    return False

########################### FUNCTIONS FOR DEBUG/DETAIL MODE ###########################

def ligands_coords_atom_index_dict(mols):
    """For the debug/detail mode only; creates a dictionary of dictionaries - each ligand's name and the dictionary of its coords with their atom indices and number of ligand's occurrence in sdf file as values

    :param mols: list of Pybel-parsed ligands' objects
    :type mole: list
    :return: dictionary indexed by ligand name, with the coords as subdictionary keys and their atom indices and number of ligand's occurrence in sdf file as values
    :rtype: dict
    """

    dictionary = {}

    for i in range(len(mols)):

        name = get_ligand_name_pose(dictionary, mols[i].title)
        dictionary[name] = {}

        for atom in mols[i]:
            dictionary[name][atom.coords] = (atom.idx, i+1)

    return dictionary

def rna_coords_atom_index_dict(structure):
    """For the debug/detail mode only; creates a dictionary or nucleic acid's atoms' coords with their atom ids as values

    :param structure: nucleic acid structure object
    :type structure: pybel.Molecule
    :return: dictionary indexed by nucleic acid's atoms' coords, with their atom ids as values
    :rtype: dict
    """

    dictionary = {}

    for residue in openbabel.OBResidueIter(structure.OBMol):
        for atom in openbabel.OBResidueAtomIter(residue):
            dictionary[(atom.GetX(), atom.GetY(), atom.GetZ())] = atom.GetResidue().GetAtomID(atom).strip()

    return dictionary

def print_debug_info(ligands_hba_hbd, ligands_HAL, ligands_CA, arom_ring_ligands_info, debug_dict_ligand,
RNA_HB_acc_don_info, RNA_anion_info, arom_RNA_ligands_info, HB_RNA_acc_info, HB_RNA_donor_info,
HAL_info, Cation_Anion_info, Pi_Cation_info, Pi_Anion_info, Sandwich_Displaced_info, T_shaped_info, columns):
    """Prints all collected information in debug mode of SIFt type FULL/XP about ligands/nucleic acid properties and detected interactions.

    :param ligands_hba_hbd: dictionary indexed by ligand name, with the coords od all ligand's hydrogen bonds acceptors & donors
    :type ligands_hba_hbd: dict
    :param ligands_HAL: dictionary indexed by ligand name, with the coords od all ligand's halogens & halogen bonds donors
    :type ligands_HAL: dict
    :param ligands_CA: dictionary indexed by ligand name, with the coords od all ligand's cations & anions
    :type ligands_CA: dict
    :param arom_ring_ligands_info: dictionary indexed by ligand name, with set of indices of atoms building ligand's aromatic rings
    :type arom_ring_ligands_info: dict
    :param debug_dict_ligand: dictionary of dictionaries of ligand's atom's coords with their corresponding atom index - {ligand_name : {(x1,y1,z1) : 1, (x2,y2,z2) : 2}}
    :type debug_dict_ligand: dict
    :param RNA_HB_acc_don_info: dictionary of dictionaries of nucleic acid's residue's numbers and list of it's HB acceptors & HB donors {chain : {res_no : [[list of HB acceptors][list of HB donors]]}}
    :type RNA_HB_acc_don_info: dict
    :param RNA_anion_info: dictionary of dictionaries of nucleic acid's residue's numbers and list of it's Anions {chain : {res_no : [list of Anions]}}
    :type RNA_anion_info: dict
    :param arom_RNA_ligands_info: dictionary of dictionaries of nucleic acid's residue's numbers and set of it's aromatic ring's atoms IDs {chain : {res_no : {set of aromatic ring's atoms IDs}}}
    :type arom_RNA_ligands_info: dict
    :param HB_RNA_acc_info: all found Hydrogen Bonds where nucleic acid is HB acceptor
    :type HB_RNA_acc_info: str
    :param HB_RNA_donor_info: all found Hydrogen Bonds where nucleic acid is HB donor
    :type HB_RNA_donor_info: str
    :param HAL_info: all found Halogen Bonds
    :type HAL_info: str
    :param Cation_Anion_info: all found Cation-Anion interactions
    :type Cation_Anion_info: str
    :param Pi_Cation_info: all found Pi-Cation interactions
    :type Pi_Cation_info: str
    :param Pi_Anion_info: all found Pi-Anion interactions
    :type Pi_Anion_info: str
    :param Sandwich_Displaced_info: all found Pi-Stacking type Sandwich/Displaced interactions
    :type Sandwich_Displaced_info: str
    :param T_shaped_info: all found Pi-Stacking type T-shaped interactions
    :type T_shaped_info: str
    :param columns: width of terminal
    :type columns: int
    :return: prints debug info
    :rtype: None
    """

    print(('*** LIGANDS PROPERTIES ***').center(columns))
    print()
    print()

    print(('#'*len('# Acceptors & Donors atom indices #')))
    print(('# 1. HB ACCEPTORS/DONORS          #'))
    print(('# Acceptors & Donors atom indices #'))
    print(('#'*len('# Acceptors & Donors atom indices #')))

    for k in ligands_hba_hbd.keys():
        print()
        print('### {} ###'.format(k))
        print('Acceptors')
        for c0 in ligands_hba_hbd[k][0]:
            print(str(debug_dict_ligand[k][c0]), end=" ")
        print()
        print('Donors')
        # prints atom inedx for each Donor-Hydrogen pair therefore if donor has 2 hydrogens it will be printed twice
        for c1 in ligands_hba_hbd[k][1]:
            print(str(debug_dict_ligand[k][c1[0]]), end=" ")
        print()

    print(('#'*len('# Carbons & Halogens atom indices #')))
    print(('#       2. HAL DONORS             #'))
    print(('# Carbons & Halogens atom indices #'))
    print(('#'*len('# Carbons & Halogens atom indices #')))

    for k in ligands_HAL.keys():
        print()
        print('### {} ###'.format(k))
        for c0 in ligands_HAL[k]:
            tuple_c = (c0[0][0], c0[0][1], c0[0][2])
            tuple_h = (c0[1][0], c0[1][1], c0[1][2])
            print('{}&{}'.format(debug_dict_ligand[k][tuple_c], debug_dict_ligand[k][tuple_h]), end=" ")
        print()

    print()
    print(('#'*len('# Cations & Anions atom indices #')))
    print(('# 3. CATIONS & ANIONS           #'))
    print(('# Cations & Anions atom indices #'))
    print(('#'*len('# Cations & Anions atom indices #')))

    for k in ligands_CA.keys():
        print()
        print('### {} ###'.format(k))
        print('Cations')
        for c0 in ligands_CA[k][0]:
            print(str(debug_dict_ligand[k][c0]), end=" ", sep=',')
        print('\nAnions')
        for c1 in ligands_CA[k][1]:
            print(str(debug_dict_ligand[k][c1]), end=" ", sep=',')
        print()

    print()
    print(('#'*len("# Number of ligand's aromatic rings #")))
    print(('#         4. AROMATIC RINGS          #'))
    print(("# Number of ligand's aromatic rings #"))
    print(('#'*len("# Number of ligand's aromatic rings #")))

    for key in ligands_CA.keys():
        print()
        print('### {} ###'.format(key))
        if key not in arom_ring_ligands_info.keys():
            print('0')
        else:
            print(arom_ring_ligands_info[key])
        print()

    print()
    print(('*** NUCLEIC ACID PROPERTIES ***').center(columns))
    print()
    print()
    print(('#'*len('# 1. HB ACCEPTORS/DONORS #')))
    print(('# 1. HB ACCEPTORS/DONORS #'))
    print(('#'*len('# 1. HB ACCEPTORS/DONORS #')))

    for ch in RNA_HB_acc_don_info.keys():
        print()
        print('@@@ CHAIN {} @@@'.format(ch))
        for k, v in collections.OrderedDict(sorted(RNA_HB_acc_don_info[ch].items())).items():
            print()
            print('### Res {} ###'.format(k))
            print('Acceptors')
            for el in v[0]:
                print('{}'.format(el), end=' ')
            print('\nDonors')
            for el in v[1]:
                print('{}'.format(el), end=' ')
            print()

    print()
    print(('#'*len("# IDs of nucleic acid anions #")))
    print(('#         2. ANIONS          #'))
    print(('#'*len("# IDs of nucleic acid anions #")))
    print()
    print('@@@ CHAIN {} @@@'.format(ch))
    for k, v in collections.OrderedDict(sorted(RNA_anion_info[ch].items())).items():
        print()
        print('### Res {} ###'.format(k))
        for el in v:
            print('{}'.format(el), end=' ')
        print()

    print()
    print(('#'*len("# IDs of nucleic acid aromatic rings #")))
    print(('#         3. AROMATIC RINGS          #'))
    print(('#'*len("# IDs of nucleic acid aromatic rings #")))

    for ch in arom_RNA_ligands_info.keys():
        print()
        print('@@@ CHAIN {} @@@'.format(ch))
        for k, v in collections.OrderedDict(sorted(arom_RNA_ligands_info[ch].items())).items():
            print()
            print('### Res {} ###'.format(k))
            for el in v:
                print('{}'.format(el))

    print()
    print()
    print(('*** DETECTED INTERACTIONS ***').center(columns))
    print()
    print()
    print(('#'*len("# HB - nucleic acid as acceptor #")))
    print(("# HB - nucleic acid as acceptor #"))
    print(('#'*len("# HB - nucleic acid as acceptor #")))
    print((HB_RNA_acc_info))
    print(('#'*len("# HB - nucleic acid as donor #")))
    print(("# HB - nucleic acid as donor #"))
    print(('#'*len("# HB - nucleic acid as donor #")))
    print((HB_RNA_donor_info))
    print(('#'*len("# HAL #")))
    print(("# HAL #"))
    print(('#'*len("# HAL #")))
    print(HAL_info)
    print(('#'*len("# AC #")))
    print(("# AC #"))
    print(('#'*len("# AC #")))
    print(Cation_Anion_info)
    print(('#'*len("# Pi-Cation #")))
    print(("# Pi-Cation #"))
    print(('#'*len("# Pi-Cation #")))
    print(Pi_Cation_info)
    print(('#'*len("# Pi-Anion #")))
    print(("# Pi-Anion #"))
    print(('#'*len("# Pi-Anion #")))
    print(Pi_Anion_info)
    print(('#'*len("# Pi-Stacking: Sandwich/Displaced #")))
    print(("# Pi-Stacking: Sandwich/Displaced #"))
    print(('#'*len("# Pi-Stacking: Sandwich/Displaced #")))
    print(Sandwich_Displaced_info)
    print(('#'*len("# Pi-Stacking: T-shaped #")))
    print(("# Pi-Stacking: T-shaped #"))
    print(('#'*len("# Pi-Stacking: T-shaped #")))
    print(T_shaped_info)
