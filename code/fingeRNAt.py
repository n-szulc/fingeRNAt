#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
fingeRNAt is a software tool for detecting non-covalent interactions formed within complexes of nucleic acids with ligands.

Interactions are encoded and saved i.e. in the form of bioinformatic-friendly Structural Interaction Fingerprint (SIFt) - a binary string, where the respective bit in the fingerprint is set to 1 in case of a presence of a particular interaction and to 0 otherwise.
This enables high-throughput analysis of the interaction data using data analysis techniques.

Authors:\n
Natalia A. Szulc, nszulc@iimcb.gov.pl\n
Filip Stefaniak, fstefaniak@genesilico.pl

If you use this software, please cite:\n
Natalia A. Szulc, Zuzanna Mackiewicz, Janusz M. Bujnicki, Filip Stefaniak\n
[in preparation]

fingeRNAt requires Python 3.5 - 3.8
'''

import argparse
import pandas as pd
import numpy as np
import os
import copy
import shutil
from openbabel import openbabel
from openbabel import pybel
from tqdm import tqdm

# Own modules
import config
from preprocessing import measure_distance, vector, calculate_angle, calculate_planar, centroid
from preprocessing import get_ligand_name_pose, projection, find_ligands_all_atoms, assign_interactions_results
from preprocessing import wrap_results, find_ligands_HBA_HBD, find_ligands_HAL_don, find_ligands_CA
from preprocessing import find_ligands_ions, find_ligands_water, find_ligands_lipophilic, find_RNA_rings, find_RNA_HB_HAL_acc_don
from preprocessing import find_RNA_anions, check_if_RNA, findAromaticRingsWithRDKit, parseYaml, find_atoms_from_SMARTS
from preprocessing import rna_coords_atom_index_dict, rna_coords_residue_index_dict, addHwithRDKit, ligands_coords_atom_index_dict, print_debug_info


##################################################
#  FUNCTIONS CALCULATING MOLECULAR INTERACTIONS  #
##################################################

def calculate_SIMPLE(residue, ligand_name, ligand_atoms, centroid_ligand, CUTOFF):
    """Calculates SIMPLE interaction between residue - ligand pair:\n
            1. Check nucleic acid residue - ligand distance
            2. Compare the distance to CUTOFF:
                - write down 1 if the distance <= CUTOFF
                - write down 0 if the distance > CUTOFF

        :param residue: residue as OpenBabel object
        :param ligand_name: ligand_name^pose_number
        :param ligand_atoms: coordinates of all ligand's atoms (except hydrogens)
        :param CUTOFF: declared cutoff value for interaction
        :type residue: openbabel.OBResidue
        :type ligand_name: str
        :type ligand_atoms: list
        :type CUTOFF: float
        :return: [ligand_name^pose_number, residue_number:residue_chain, binary info about interaction (0/1)]
        :rtype: list
    """

    # List of sublists of all atoms' coords of the residue
    residue_atoms = []

    for atom in openbabel.OBResidueAtomIter(residue):
        if atom.GetAtomicNum() != config.HYDROGEN_NUM: # do not consider hydrogens
            residue_atoms.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))

    result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0]

    if measure_distance(centroid(residue_atoms), centroid_ligand) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
        return result

    # Flag to iterate over residue's atoms as long as we do not find an atom within CUTOFF distance from ligand
    flag=True

    for rna_atom in residue_atoms:
        if flag:
            for ligand_atom in ligand_atoms:
                dist = measure_distance(ligand_atom, rna_atom)
                if config.MIN_DIST < dist <= CUTOFF:
                        result[-1]=1 # Condition met; write down 1

                        if debug:
                            print('### {} - {} first below cutoff dist: {} ###'.format(filename_RNA.split(sys_sep)[-1], ligand_name, np.round(dist, 4)))
                            print('    between\t{}:{}:{}\t {} atom {}'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])], debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0], ligand_name))

                        if detail:
                            global detail_list
                            if filename_ligand:
                                ligand_name_detail = ligand_name.split('^')[0]
                                ligand_pose_detail = ligand_name.split('^')[1]
                            else:
                                ligand_name_detail = ligand_name.split(':')[0]
                                ligand_pose_detail = 0
                            detail_list.append([ligand_name_detail, ligand_pose_detail, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'SIMPLE',
                            debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0],
                            ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                            residue.GetName(), residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])],
                            rna_atom[0], rna_atom[1], rna_atom[2],
                            np.round(dist, 4)])

                        if not detail:
                            flag=False
                            break
        else:
            break

    return result

def calculate_PBS(residue, ligand_name, ligand_atoms, centroid_ligand, CUTOFF):
    """Calculates PBS interaction between residue - ligand pair:\n
            1. Divide nucleic acid's residue into 3 groups Phosphate/Base/Sugar (P/B/S)
            2. Check each group - ligand distance
            3. Compare the distance to CUTOFF:
                - write down 1 if the distance <= CUTOFF
                - write down 0 if the distance > CUTOFF

        :param residue: residue as OpenBabel object
        :param ligand_name: ligand_name^pose_number
        :param ligand_atoms: coordinates of all ligand's atoms (except hydrogens)
        :param CUTOFF: declared cutoff value for interaction
        :type residue: openbabel.OBResidue
        :type ligand_name: str
        :type ligand_atoms: list
        :type CUTOFF: float
        :return: [ligand_name^pose_number, residue_number:residue_chain, binary info about interaction in P group (0/1), binary info about interaction in B group (0/1), binary info about interaction in S group (0/1)]
        :rtype: list
    """

    # List of residue's atoms as OBAtoms objects
    residue_atoms = []
    residue_atoms_coords = []

    for atom in openbabel.OBResidueAtomIter(residue):
        if atom.GetAtomicNum() != config.HYDROGEN_NUM: # do not consider hydrogens
            residue_atoms.append(atom)
            residue_atoms_coords.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))


    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0, 0, 0]

    if measure_distance(centroid(residue_atoms_coords), centroid_ligand) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
        return result

    # Flag to iterate over residue's atoms as long as we do not find an atom of the defined P/B/S group within CUTOFF distance from ligand
    flags=[True, True, True]

    for rna_atom in residue_atoms:
        # Flag to check if we deal with atom that does not belong to any of the defined P/B/S groups
        atom_group_type = False
        rna_atom_name = residue.GetAtomID(rna_atom).strip()
        rna_atom_coords = np.array([rna_atom.GetX(), rna_atom.GetY(), rna_atom.GetZ()])

        for g in range(len(config.GROUPS)):
            # If atom present in the group
            if rna_atom_name in config.GROUPS[g]:
                atom_group_type = True # If atom belongs to one of the 3 groups
                if flags[g]:
                    for ligand_atom in ligand_atoms:
                        dist = measure_distance(ligand_atom, rna_atom_coords)
                        if config.MIN_DIST < dist <= CUTOFF:
                                result[g-3]= 1 # Condition met; write down 1

                                if debug:
                                    print('### {} - {} first below cutoff {} GROUP dist: {} ###'.format(filename_RNA.split(sys_sep)[-1], ligand_name, config.WHICH_GROUP[g], np.round(dist, 4)))
                                    print('    between\t{}:{}:{}\t {} atom {}'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2])], debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0], ligand_name))

                                if detail:
                                    global detail_list
                                    if filename_ligand:
                                        ligand_name_detail = ligand_name.split('^')[0]
                                        ligand_pose_detail = ligand_name.split('^')[1]
                                    else:
                                        ligand_name_detail = ligand_name.split(':')[0]
                                        ligand_pose_detail = 0

                                    detail_list.append([ligand_name_detail, ligand_pose_detail, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'PBS',
                                    debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0],
                                    ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                                    residue.GetName(), residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2])],
                                    rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2],
                                    np.round(dist, 4)])

                                if not detail:
                                    flags[g]=False
                                    break

        if not atom_group_type:
            raise Exception('Unknown atom type %s' %rna_atom_name)

    return result

def calculate_HB(residue, acceptors_RNA, donors_RNA, ligand_name, ligand_donors_acceptors, check_dha=False):
    """Calculates hydrogen bond between residue - ligand pair.\n
        Simplified graphical representation:\n
        A \*\*\*\*\* H --- D\n
            where:\n
            - A  - hydrogen bond acceptor
            - H  - hydrogen
            - D  - hydrogen bond acceptor
            - \*\*\*\*\* - hydrogen bond\n
        Default geometric rule:
            1. D-A distance < 3.9 A\n
            - If check_dha is True:\n
              2. 100 < D-H-A angle < 260\n

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors (list of sublists [[acceptor,acceptor',acceptor'(if the 2nd one exists)]])
        :param donors_RNA: list of tuples tuples (Donor, Hydrogen) of the residue
        :param ligand_name: ligand_name^pose_number
        :param ligand_donors_acceptors: [[sublist of acceptors coords], [sublist of tuples (D, H) coords]]
        :param check_dha: flag indictacting if D-H-A angle should also be considered
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type donors_RNA: list
        :type ligand_name: str
        :type ligand_donors_acceptors: list
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0]

    # List of ligand's (D, H) tuples
    ligand_donors_coords = ligand_donors_acceptors[1]
    # List of ligand's acceptors
    ligand_acceptors_coords = ligand_donors_acceptors[0]
    # Important for 'FULL' fingerprint as we are searching only for the first hydrogen bond
    searching_flag = True

    if detail:
        global detail_list

    for RNA_acceptor_set in acceptors_RNA:

        RNA_acc = RNA_acceptor_set[0] # We do not need coords of acceptor's neighbors

        if searching_flag:
            RNA_acc_coords = np.array([RNA_acc.GetX(), RNA_acc.GetY(), RNA_acc.GetZ()])

            for donor in ligand_donors_coords:
                dist = measure_distance(donor[0], RNA_acc_coords) # Measure D-A distance

                if config.MIN_DIST < dist < config.MAX_HB_DIST:
                    interaction_found = False
                    if not check_dha:
                        result[-1] = 1
                        interaction_found = True
                    else:
                        dh = vector(donor[0], donor[1])
                        ha = vector(RNA_acc_coords, donor[1])
                        angle = calculate_angle(dh, ha)
                        if config.MIN_HB_ANGLE < angle < config.MAX_HB_ANGLE:
                            result[-1] = 1
                            interaction_found = True

                    if interaction_found and debug:
                        global HB_RNA_acc_info
                        HB_RNA_acc_info += '***\n'
                        if not check_dha:
                            HB_RNA_acc_info += ('{} acceptor - {} donor\ndist: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4)))
                        else:
                            HB_RNA_acc_info +=('{} acceptor - {} donor\ndist: {}; angle: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(angle, 4)))
                        HB_RNA_acc_info += ('{}:{}:{}\t{} atom of {}\n'.format(residue.GetChain(), residue.GetNum(), residue.GetAtomID(RNA_acc).strip(), debug_dict_ligand[ligand_name][donor[0]][0], ligand_name))

                    if interaction_found and detail:
                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][donor[0]][1], 'HB',
                        debug_dict_ligand[ligand_name][donor[0]][0],
                        donor[0][0], donor[0][1],  donor[0][2],
                        residue.GetName(), residue.GetNum(), residue.GetChain(), residue.GetAtomID(RNA_acc).strip(),
                        RNA_acc_coords[0], RNA_acc_coords[1], RNA_acc_coords[2],
                        dist])

                    if interaction_found and not detail:
                        searching_flag = False
                        break
        else: break

    if searching_flag:

        for RNA_don in donors_RNA:

            if searching_flag:

                RNA_don_coords = np.array([RNA_don[0].GetX(), RNA_don[0].GetY(), RNA_don[0].GetZ()])
                if check_dha:
                    RNA_donH_coords = np.array([RNA_don[1].GetX(), RNA_don[1].GetY(), RNA_don[1].GetZ()])

                for acceptor in ligand_acceptors_coords:
                    dist = measure_distance(RNA_don_coords,acceptor) # Measure D-A distance

                    if config.MIN_DIST < dist < config.MAX_HB_DIST:
                        interaction_found = False
                        if not check_dha:
                            result[-1] = 1
                            interaction_found = True
                        else:
                            dh = vector(RNA_don_coords, RNA_donH_coords)
                            ha = vector(acceptor, RNA_donH_coords)
                            angle = calculate_angle(dh, ha)
                            if config.MIN_HB_ANGLE < angle < config.MAX_HB_ANGLE:
                                result[-1] = 1
                                interaction_found = True

                        if interaction_found and debug:
                            global HB_RNA_donor_info
                            HB_RNA_donor_info += '***\n'
                            if not check_dha:
                                HB_RNA_donor_info +=('{} donor - {} acceptor\ndist: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4)))
                            else:
                                HB_RNA_donor_info += ('{} donor - {} acceptor\ndist: {}; angle: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(angle, 4)))
                            HB_RNA_donor_info +=('{}:{}:{}\t{} atom of {}\n'.format(residue.GetNum(), residue.GetChain(), residue.GetAtomID(RNA_don[0]), debug_dict_ligand[ligand_name][acceptor][0], ligand_name))

                        if interaction_found and detail:
                            detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][acceptor][1], 'HB',
                            debug_dict_ligand[ligand_name][acceptor][0],
                            acceptor[0], acceptor[1],  acceptor[2],
                            residue.GetName(), residue.GetNum(), residue.GetChain(), residue.GetAtomID(RNA_don[0]),
                            RNA_don_coords[0], RNA_don_coords[1], RNA_don_coords[2],
                            dist])

                        if interaction_found and not detail:
                            searching_flag = False
                            break

            else: break

    return result

def calculate_HAL(residue, acceptors_RNA, ligand_name, ligand_donors_coords):
    """Calculates halogen bond between residue - ligand pair.\n
        Simplified graphical representation:\n
        Y --- O \*\*\*\*\* X --- C\n
            where:\n
            - Y  - acceptor'; atom covalently bond to the acceptor
            - O  - halogen bond acceptor
            - X  - halogen [F, Br, I, Cl]
            - C  - halogen donor: Carbon
            - \*\*\*\*\* - halogen bond\n
        .. note::
            There may be two Ys, if O is part of the ring. If so, we need to apply above rules to both Ys.\n
        Default geometric rules:
            - X-O distance < 4.0 A
            - C-X-O angle ~ 165 +/- 30
            - X-O-Y angle ~ 120 +/- 30\n

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors (list of sublists [[acceptor,acceptor',acceptor'(if the 2nd one exists)]])
        :param ligand_name: ligand_name^pose_number
        :param ligand_donors_coords: [list of tuples (C, halogen)]
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type ligand_name: str
        :type ligand_donors_coords: list
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0]

    searching_flag = True

    for RNA_acceptor_set in acceptors_RNA:

        if searching_flag:

            RNA_acc_coords = np.array([RNA_acceptor_set[0].GetX(), RNA_acceptor_set[0].GetY(), RNA_acceptor_set[0].GetZ()])
            # Flag to mark if we have found halogen bonding for one particular acceptor
            acc_bond_found = False

            for donor in ligand_donors_coords:

                if not acc_bond_found:

                    for y in range(len(RNA_acceptor_set[1:])): # For all Y's neighbours (max 2)

                        RNA_acc_y_coords = np.array([RNA_acceptor_set[y+1].GetX(), RNA_acceptor_set[y+1].GetY(), RNA_acceptor_set[y+1].GetZ()]) # coords of nucleic acid acceptor' - Y
                        dist = measure_distance(donor[1], RNA_acc_coords) # Measure X-O distance

                        if config.MIN_DIST < dist < config.MAX_HAL_DIST:

                            dh = vector(donor[0], donor[1])
                            ha = vector(RNA_acc_coords, donor[1])
                            ah = vector(donor[1], RNA_acc_coords)
                            aa = vector(RNA_acc_y_coords, RNA_acc_coords)
                            angle_acc = calculate_angle(dh, ha) # Calculate C-X-O angle
                            angle_don = calculate_angle(ah, aa) # Calculate X-O-Y angle

                            if (abs(angle_acc - config.HALOGEN_ACC_ANGLE) < config.HALOGEN_ANGLE_DEV) and (abs(angle_don - config.HALOGEN_DON_ANGLE) < config.HALOGEN_ANGLE_DEV):
                                result[-1] = 1

                                if not detail:
                                    searching_flag = False # Just found first halogen bond, no need to search further

                                if debug:
                                    global HAL_info
                                    HAL_info += '***\n'
                                    HAL_info += ('{} acceptor - {} donor\ndist: {}; C-X-O angle: {}; X-O-Y angle: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(angle_acc, 4), round(angle_don, 4)))
                                    HAL_info += ('{}:{}:{}\t{} atom of {}\n'.format(residue.GetChain(), residue.GetNum(), residue.GetAtomID(RNA_acceptor_set[0]).strip(), debug_dict_ligand[ligand_name][(donor[0][0], donor[0][1], donor[0][2])][0], ligand_name))

                                if detail:
                                    global detail_list
                                    detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][(donor[0][0], donor[0][1], donor[0][2])][1], 'HAL',
                                    debug_dict_ligand[ligand_name][(donor[0][0], donor[0][1], donor[0][2])][0],
                                    donor[0][0], donor[0][1],  donor[0][2],
                                    residue.GetName(), residue.GetNum(), residue.GetChain(), residue.GetAtomID(RNA_acceptor_set[0]).strip(),
                                    RNA_acc_coords[0], RNA_acc_coords[1], RNA_acc_coords[2],
                                    dist])

                                # If we found halogen bond for one O-Y pair, there is no need to check angles for another Y of the same O (if O has 2 neighbours)
                                acc_bond_found = True
                                break
        else:
           break

    return result

def calculate_CATION_ANION(residue, RNA_anions, ligand_name, ligand_cation_coords):
    """ Calculates cation-anion interaction between residue - ligand pair.\n
        Simplified graphical representation:\n
        C \*\*\*\*\* A\n
            where:\n
            - C  - cation
            - A  - anion\n
        Default geometric rule:
            - 0.5 A < cation-anion distance < 5.5 A\n

        :param residue: residue as OpenBabel object
        :param RNA_anions: residue's anions coordinates [OP1, OP2]
        :param ligand_name: ligand_name^pose_number
        :param ligand_cation_coords: list of ligand's cations coords
        :type residue: openbabel.OBResidue
        :type RNA_anions: list
        :type ligand_name: str
        :type ligand_cation_coords: list
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0]

    searching_flag = True

    for anion in anions_RNA:

        if searching_flag:
            RNA_anion_coords = np.array([anion.GetX(), anion.GetY(), anion.GetZ()])

            for cation in ligand_cation_coords:

                dist = measure_distance(cation, RNA_anion_coords) # Measure cation-anion distance

                if config.MIN_DIST < dist < config.MAX_CA_DIST:
                    result[-1] = 1

                    if debug:
                        global Cation_Anion_info
                        Cation_Anion_info += '***\n'
                        Cation_Anion_info += ('{} - {}\ndist: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4)))
                        Cation_Anion_info += ('{}:{}:{}\t{} atom of {}\n'.format(residue.GetChain(), residue.GetNum(), residue.GetAtomID(anion).strip(), debug_dict_ligand[ligand_name][cation][0], ligand_name))

                    if detail:
                        global detail_list
                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][cation][1], 'CA',
                        debug_dict_ligand[ligand_name][cation][0],
                        cation[0], cation[1],  cation[2],
                        residue.GetName(), residue.GetNum(), residue.GetChain(), residue.GetAtomID(anion).strip(),
                        RNA_anion_coords[0], RNA_anion_coords[1], RNA_anion_coords[2],
                        dist])

                    if not detail:
                        searching_flag = False # Just found first cation-anion interaction, no need to search further
                        break
        else:
            break

    return result

def calculate_ANION_PI(residue, RNA_anions, ligand_name, ligand_rings_coords):
    """ Calculates anion-Pi interaction between residue - ligand pair.\n
        Simplified graphical representation:\n
        A \*\*\*\*\* R\n
            where:\n
            - A  - nucleic acids's anion
            - R  - ligand's aromatic ring\n
        Default geometric rule:
            - Anion - aromatic ring's center distance < 6.0 A
            - Angle between anion and aromatic ring's planar ~ 90 +/- 30\n

        :param residue: residue as OpenBabel object
        :param RNA_anions: residue's anions coordinates [OP1, OP2]
        :param ligand_name: ligand_name^pose_number
        :param ligand_rings_coords: list of ligand's aromatic rings coords
        :type residue: openbabel.OBResidue
        :type RNA_anions: list
        :type ligand_name: str
        :type ligand_rings_coords: list
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0]

    searching_flag = True

    for anion in anions_RNA:

        if searching_flag:
            RNA_anion_coords = np.array([anion.GetX(), anion.GetY(), anion.GetZ()])

            for aring in ligand_rings_coords:
                ring_atoms_ligand = aring
                ring_center_ligand = centroid([ra for ra in ring_atoms_ligand])
                dist = measure_distance(RNA_anion_coords, ring_center_ligand)

                if config.MIN_DIST < dist < config.PI_ION_DISTANCE: # Measure anion - ligand's ring center distance
                    atoms_creating_planar_space_ligand = np.array([ring_atoms_ligand[0], ring_atoms_ligand[1], ring_atoms_ligand[2]], dtype=np.longdouble) # Add 3 atoms (we do not need more) from ring to calculate planar
                    planar_ligand = calculate_planar(atoms_creating_planar_space_ligand)
                    ligand_ring_center = vector(RNA_anion_coords, ring_center_ligand)
                    angle = calculate_angle(ligand_ring_center, planar_ligand) # Calculate angle between cation/anion-ring center and aromatic ring's planar

                    if angle > 90: angle = 180 - angle # We are 'on the other side'
                    anion_pi_angle = 90 - angle # The angle is equal to the complementary acute angle

                    if abs(config.PI_ION_ANGLE - anion_pi_angle) < config.PI_ION_ANGLE_DEV:
                        result[-1] = 1

                        if debug or detail:
                            debug_ligand = ''
                            for i in range(len(ring_atoms_ligand)):
                                debug_ligand += str(debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][0])
                                if i != len(ring_atoms_ligand)-1:
                                    debug_ligand += ','

                            if debug:
                                global Anion_Pi_info
                                Anion_Pi_info += '***\n'
                                Anion_Pi_info += ("{} - {}\ndist: {}; anion to ring's planar angle: {}\n".format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(anion_pi_angle, 4)))
                                Anion_Pi_info += ('{}:{}:{}\t{} atoms of {}\n'.format(residue.GetChain(), residue.GetNum(), residue.GetAtomID(anion).strip(), debug_ligand, ligand_name))

                            if detail:
                                global detail_list
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][1], 'Pi_Anion',
                                debug_ligand,
                                ring_center_ligand[0], ring_center_ligand[1],  ring_center_ligand[2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), residue.GetAtomID(anion).strip(),
                                RNA_anion_coords[0], RNA_anion_coords[1], RNA_anion_coords[2],
                                dist])

                        if not detail:
                            searching_flag = False # Just found first cation-anion interaction, no need to search further
                            break
        else:
            break

    return result


def calculate_PI_INTERACTIONS(RNA_rings, RNA_all_atoms, all_ligands_CA_dict, all_ligands_rings_dict):
    """Calculates Pi-cation, Pi-anion & Pi-stacking interactions between all residues' aromatic rings - all ligands' pairs.\n
        .. note::
            Each ring of purines is considered separately.\n
        Pi-cation & Pi-anion default geometric rules:
            - Cation/anion - aromatic ring's center distance < 6.0 A
            - Angle between aromatic ring's planar and cation/anion ~ 90 +/- 30\n
        Three types of Pi-stacking interactions:
            - Sandwich
            - Parallel Displaced
            - T-shaped\n
        Pi-stacking default geometric rules:
            - All types' common rules:
                - Aromatic ring center - aromatic ring center distance < 5.5 A
                - Aromatic rings' offset < 2.0 A\n
            - For sandwich & parallel displaced types:
                - Angle between aromatic rings' planars < 30\n
            - For T-shaped type:
                - Angle between aromatic rings' planars ~ 90 +/- 30\n\n

        :param RNA_rings: list of all RNA aromatic rings found by OpenBabel
        :param RNA_all_atoms: list of all RNA atoms
        :param all_ligands_CA_dict: dictionary containing ligand's cations & anions {'prefix^pose':[[list of cations coords],[list of anions coords]]}
        :param all_ligands_rings_dict: dictionary containing ligand's aromatic rings {'prefix^pose':[[list of atoms in aromatic ring coords]]}
        :type RNA_rings: list
        :type RNA_all_atoms: list
        :type all_ligands_CA_dict: dict
        :type all_ligands_rings_dict: dict
        :return: calculated 3 Pi-interaction for RNA - all ligands
        :rtype: list
       """

    #########################################
    #  Common part for all Pi-interactions  #
    #########################################

    # There will be 3 results of Pi-interactions: : Pi-cation, Pi-anion, Pi-stacking
    RESULTS = [[],[],[]]

    if detail:
        global detail_list

    # Looping over nucleic acid rings
    for ring in RNA_rings: # Unlike in previous functions, iteration is over all nucleic acid rings

        ring_atoms_RNA = [a for a in RNA_all_atoms if ring.IsMember(a.OBAtom)]
        atoms_creating_planar_space_RNA = np.array([ring_atoms_RNA[0].coords,ring_atoms_RNA[1].coords,ring_atoms_RNA[2].coords],dtype=np.longdouble) # add 3 atoms (we do not need more) from nucleic acid's ring to calculate planar
        planar_RNA = calculate_planar(atoms_creating_planar_space_RNA)
        residue = structure.OBMol.GetAtom(ring_atoms_RNA[0].idx).GetResidue() # Get nucleic acid's ring's residue
        ring_center_RNA = centroid([ra.coords for ra in ring_atoms_RNA])

        results = [[],[],[]] # There will 3 be results for each nucleic acid's residue from 3 Pi-interactions: : Pi-cation, Pi-anion, Pi-stacking

        if debug:
            global arom_RNA_ligands_info
            if residue.GetChain() not in arom_RNA_ligands_info.keys():
                arom_RNA_ligands_info[residue.GetChain()] = {}
            if residue.GetNum() not in arom_RNA_ligands_info[residue.GetChain()].keys():
                arom_RNA_ligands_info[residue.GetChain()][residue.GetNum()] =  set()
            arom_RNA_ligands_info[residue.GetChain()][residue.GetNum()].add(','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]))

        #################################################
        #  Calculate Pi-Cation & Pi-Anion interactions  #
        #################################################

        for ligand_name in all_ligands_CA_dict.keys():

            ligand_ion_coords = [all_ligands_CA_dict[ligand_name][0], all_ligands_CA_dict[ligand_name][1]]  # Calculate Pi-cations & Pi-anions interactions

            for j in range(2): # Pi-cation & Pi-anion interactions

                results[j].append([ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0])

                for ion in ligand_ion_coords[j]:
                    dist = measure_distance(ion, ring_center_RNA)
                    if config.MIN_DIST < dist < config.PI_ION_DISTANCE: # Measure ring center-cation/anion distance

                        ion_ring_center = vector(ring_center_RNA, ion)
                        angle = calculate_angle(ion_ring_center, planar_RNA) # Calculate angle between cation/anion-ring center and aromatic ring's planar

                        if angle > 90: angle = 180 - angle # We are 'on the other side'
                        pi_ion_angle = 90 - angle # The angle is equal to the complementary acute angle

                        if abs(config.PI_ION_ANGLE - pi_ion_angle) < config.PI_ION_ANGLE_DEV:
                            results[j][-1][-1] = 1

                            if debug:
                                if j == 0:
                                    global Pi_Cation_info
                                    Pi_Cation_info += '***\n'
                                    Pi_Cation_info +=("{} - {}\ndist: {}; ring's planar to Cation angle: {}\n".format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(pi_ion_angle, 4)))
                                    Pi_Cation_info +=('{}:{}:{}\t{} atom of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_dict_ligand[ligand_name][ion][0], ligand_name))
                                else:
                                    global Pi_Anion_info
                                    Pi_Anion_info += '***\n'
                                    Pi_Anion_info +=("{} - {}\ndist: {}; ring's planar to Anion angle: {}\n".format(filename_RNA.split(sys_sep)[-1], ligand_name, round(dist, 4), round(pi_ion_angle, 4)))
                                    Pi_Anion_info +=('{}:{}:{}\t{} atom {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_dict_ligand[ligand_name][ion][0], ligand_name))

                            if detail:
                                if j == 0:
                                    detail_inter = 'Pi_Cation'
                                else:
                                    detail_inter = 'Pi_Anion'
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][ion][1], detail_inter,
                                debug_dict_ligand[ligand_name][ion][0],
                                ion[0], ion[1],  ion[2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]),
                                ring_center_RNA[0], ring_center_RNA[1], ring_center_RNA[2],
                                measure_distance(ion, ring_center_RNA)])

        ########################################
        #  Calculate Pi-Stacking interactions  #
        ########################################

        for ligand_name in all_ligands_rings_dict.keys():

            results[2].append([ligand_name,  str(residue.GetNum()) + ':' + str(residue.GetChain()), 0])

            ligand_rings = all_ligands_rings_dict[ligand_name] # Take list of ligand's aromatic rings atoms coords

            for aring in ligand_rings:

                if debug:
                    global arom_ring_ligands_info
                    arom_ring_ligands_info[ligand_name] = []
                    debug_ligand = ''
                    for i in range(len(aring)):
                        debug_ligand += str(debug_dict_ligand[ligand_name][aring[i]][0])
                        if i != len(aring) - 1:
                            debug_ligand += ','

                    arom_ring_ligands_info[ligand_name].append(debug_ligand)

                ring_atoms_ligand = aring
                atoms_creating_planar_space_ligand = np.array([ring_atoms_ligand[0], ring_atoms_ligand[1], ring_atoms_ligand[2]], dtype=np.longdouble) # Add 3 atoms (we do not need more) from ring to calculate planar
                planar_ligand = calculate_planar(atoms_creating_planar_space_ligand)
                ring_center_ligand = centroid([ra for ra in ring_atoms_ligand])
                centroid_distance = measure_distance(ring_center_RNA, ring_center_ligand) # Measure ring - ring distance

                # Calculate aromatic rings' center offset (project each ring center into the other ring)
                proj1 = projection(planar_ligand, ring_center_ligand, ring_center_RNA)
                proj2 = projection(planar_RNA, ring_center_RNA, ring_center_ligand)
                offset = min(measure_distance(proj1, ring_center_ligand), measure_distance(proj2, ring_center_RNA))

                if config.MIN_DIST < centroid_distance < config.RING_RING_MAX and offset < config.PISTACK_OFFSET_MAX:
                    planar_angle = calculate_angle(planar_RNA, planar_ligand) # Calculate planar - planar angle, which is equal to angle between normal vectors of both planes

                    if planar_angle > 90: planar_angle = 180 - planar_angle # We are 'on the other side'

                    if planar_angle < config.PI_ANGLE_DISPLACED: # Sandwich & Displaced Pi-stacking interactions
                        results[2][-1][-1] = 1

                        if debug or detail:
                            debug_ligand = ''
                            for i in range(len(ring_atoms_ligand)):
                                debug_ligand += str(debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][0])
                                if i != len(ring_atoms_ligand) - 1:
                                    debug_ligand += ','

                            if debug:
                                global Sandwich_Displaced_info
                                Sandwich_Displaced_info += '***\n'
                                Sandwich_Displaced_info += "{} - {}\nrings center dist: {}; rings offset: {}; rings planars angle: {}\n".format(filename_RNA.split(sys_sep)[-1], ligand_name, round(centroid_distance, 4), round(offset, 4), round(planar_angle, 4))
                                Sandwich_Displaced_info += '{}:{}:{}\t{} atoms of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_ligand, ligand_name)

                            if detail:
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][1], 'Pi_Stacking',
                                debug_ligand,
                                ring_center_ligand[0], ring_center_ligand[1],  ring_center_ligand[2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]),
                                ring_center_RNA[0], ring_center_RNA[1], ring_center_RNA[2],
                                centroid_distance])

                    elif abs(config.PLANAR_ANGLE_TSHAPED - planar_angle) < config.PLANAR_ANGLE_TSHAPED_DEV: # T-shaped Pi-stacking interaction
                        results[2][-1][-1] = 1

                        if debug or detail:
                            debug_ligand = ''
                            for i in range(len(ring_atoms_ligand)):
                                debug_ligand += str(debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][0])
                                if i != len(ring_atoms_ligand)-1:
                                    debug_ligand += ','

                            if debug:
                                global T_shaped_info
                                Sandwich_Displaced_info += '***\n'
                                Sandwich_Displaced_info += "{} - {}\nrings center dist: {}; rings offset: {}; rings planars angle: {}\n".format(filename_RNA.split(sys_sep)[-1], ligand_name, round(centroid_distance, 4), round(offset, 4), round(planar_angle, 4))
                                Sandwich_Displaced_info += '{}:{}:{}\t{} atoms of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_ligand, ligand_name)

                            if detail:
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][ring_atoms_ligand[i]][1], 'Pi_Stacking',
                                debug_ligand,
                                ring_center_ligand[0], ring_center_ligand[1],  ring_center_ligand[2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]),
                                ring_center_RNA[0], ring_center_RNA[1], ring_center_RNA[2],
                                centroid_distance])

                    else:
                        continue

        ###########################################################
        #  Merge all the Pi-interactions results and return them  #
        ###########################################################

        for j in range(3): # Append results from calculated 3 Pi-interactions: Pi-cation, Pi-anion, Pi-stacking
            RESULTS[j].extend(results[j])

    return RESULTS


def calculate_ION_MEDIATED(residue, residue_atoms, ligand_name, ions, ions_dict):
    """ Calculates ion-mediated ligand-residue interaction.\n
        Simplified graphical representation:\n
        R \*\*\*\*\* I \*\*\*\*\* A\n
            where:\n
            - R  - nucleic acid residue (nitrogen or oxygen atom)
            - I  - inorganic ion
            - A  - ligand (nitrogen, oxygen or sulphur atom)\n
        Default geometric rule:
            - 0.5 A < ion-anion distance <= 4.0 A\n
            - 0.5 A < residue-ion distance <= 4.0 A\n

        :param residue: residue as OpenBabel object
        :param residue_atoms: residue's nitrogen or oxygen atoms' coordinates
        :param ligand_name: ligand_name^pose_number
        :param ions: list of ions in contact with ligand's nitrogen, oxygen or sulphur atoms
        :param ions_dict: dictionary with ions' names and coordinates as values
        :type residue: openbabel.OBResidue
        :type residue_atoms: list
        :type ligand_name: str
        :type ions: list
        :type ions_dict: dict
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0, 0, 0, 0]

    # Flag to iterate over residue's atoms as long as we do not find an atom within CUTOFF distance from ligand
    searching_flag = True

    for rna_atom in residue_atoms:
        if searching_flag:
            for ion in ions:
                dist = measure_distance(ions_dict[ion], rna_atom)
                if config.MIN_DIST < dist <= config.MAX_ION_DIST:

                        ion_name = ion.split(':')[0]
                        found_interaction = True
                        if ion_name == 'MG' and dist <= config.MAX_MAGNESIUM_DIST:
                            result[-4] = 1

                        elif ion_name == 'K' and dist <= config.MAX_POTASSIUM_DIST:
                            result[-3] = 1
                        elif ion_name == 'NA' and dist <= config.MAX_SODIUM_DIST:
                            result[-2] = 1
                        elif ion_name not in ['MG', 'K', 'NA'] and dist <= config.MAX_OTHER_ION_DIST:
                            result[-1] = 1
                        else:
                            found_interaction = False

                        if found_interaction and (debug or detail):

                            shortest_ligand_ion = 9999
                            ligand_atom = None
                            for atom in ligands_all_atoms[ligand_name]:
                                ligand_ion_dist = measure_distance(atom, ions_dict[ion])
                                if ligand_ion_dist < shortest_ligand_ion:
                                    shortest_ligand_ion = ligand_ion_dist
                                    ligand_atom = atom

                            if debug:
                                global ion_mediated_info
                                ion_mediated_info += '***\n'
                                ion_mediated_info += '{} - {} - {}\n shortest dist: {}\t{}\n'.format(filename_RNA.split(sys_sep)[-1],  ion, ligand_name, shortest_ligand_ion, np.round(dist, 4))
                                ion_mediated_info += '{}:{}:{}\t ion {}\t{} atom {}\n'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])], ion, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0], ligand_name)

                            if detail:
                                global detail_list
                                ligand_name_detail = ligand_name.split('^')[0]
                                ligand_pose_detail = ligand_name.split('^')[1]

                                detail_list.append([ligand_name_detail, ligand_pose_detail, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'Ion-mediated',
                                debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0],
                                ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                                ion.split(':')[0], ion.split(':')[1], ion.split(':')[2], ion.split(':')[0],
                                ions_dict[ion][0][0], ions_dict[ion][0][1], ions_dict[ion][0][2],
                                shortest_ligand_ion])

                                detail_list.append([ion.split(':')[0], 0, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'Ion-mediated', ion.split(':')[1] + ':' + ion.split(':')[2],
                                ions_dict[ion][0][0], ions_dict[ion][0][1],  ions_dict[ion][0][2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])],
                                rna_atom[0], rna_atom[1], rna_atom[2],
                                np.round(dist, 4)])

                        if found_interaction and not detail:
                            searching_flag = False
                            break
        else:
            break

    return result

def calculate_WATER_MEDIATED(residue, acceptors_RNA, donors_RNA, ligand_name, water_molecules, water_dict):
    """ Calculates water-mediated ligand-residue interaction.\n
        Simplified graphical representation:\n
        R \*\*\*\*\* W \*\*\*\*\* L\n
            where:\n
            * R  - nucleic acid's hydrogen bond donor/acceptor
            * W  - water molecule (oxygen)
            * A  - ligand's hydrogen bond donor/acceptor\n
        Default geometric rule:
            - 0.5 A < water - ligand distance <= 3.5 A
            - 0.5 A < residue - water distance <= 3.5 A\n

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors (list of sublists [[acceptor,acceptor',acceptor'(if the 2nd one exists)]])
        :param donors_RNA: list of tuples tuples (Donor, Hydrogen) of the residue
        :param ligand_name: ligand_name^pose_number
        :param water_molecules: list of water molecules in contact with ligand
        :param water_dict: dictionary with water names and coordinates as values
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type donors_RNA: list
        :type ligand_name: str
        :type water_molecules: list
        :type water_dict: dict
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0]
    searching_flag = True
    acceptors_donors_RNA = [item[0] for item in acceptors_RNA] + [item[0] for item in donors_RNA]

    for rna_object in acceptors_donors_RNA:
        if searching_flag:
            rna_atom = np.array([rna_object.GetX(), rna_object.GetY(), rna_object.GetZ()])
            for water in water_molecules:
                dist = measure_distance(water_dict[water], rna_atom)

                if config.MIN_DIST < dist <= config.MAX_WATER_DIST:
                        result[-1] = 1

                        if debug or detail:
                            shortest_ligand_water = 9999
                            ligand_atom = None
                            for atom in ligands_all_atoms[ligand_name]:
                                ligand_water_dist = measure_distance(atom, water_dict[water])
                                if ligand_water_dist < shortest_ligand_water:
                                    shortest_ligand_water = ligand_water_dist
                                    ligand_atom = atom

                            if debug:
                                global water_mediated_info
                                water_mediated_info += '***\n'
                                water_mediated_info += '{} - {} - {}\n shortest dist: {}\t{}\n'.format(filename_RNA.split(sys_sep)[-1],  water, ligand_name, shortest_ligand_water, np.round(dist, 4))
                                water_mediated_info += '{}:{}:{}\t ion {}\t{} atom {}\n'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])], water, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0], ligand_name)

                            if detail:
                                global detail_list
                                ligand_name_detail = ligand_name.split('^')[0]
                                ligand_pose_detail = ligand_name.split('^')[1]

                                detail_list.append([ligand_name_detail, ligand_pose_detail, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'Water-mediated',
                                debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][0],
                                ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                                water.split(':')[0], water.split(':')[1], water.split(':')[2], "O",
                                water_dict[water][0][0], water_dict[water][0][1], water_dict[water][0][2],
                                shortest_ligand_water])

                                detail_list.append([water.split(':')[0], 0, debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])][1], 'Water-mediated', water.split(':')[1] + ':' +water.split(':')[2],
                                water_dict[water][0][0], water_dict[water][0][1],  water_dict[water][0][2],
                                residue.GetName(), residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])],
                                rna_atom[0], rna_atom[1], rna_atom[2],
                                np.round(dist, 4)])

                        if not detail:
                            flag = False
                            break
        else:
            break

    return result

def calculate_lipophilic_interactions(residue, residue_atoms, ligand_name, ligands_lipophilic_coords):
    """ Calculates lipohilic ligand-residue interaction.\n
        1. Check nucleic acid residue (carbon atoms only) - ligand (detected lipohilic atoms) distance
        2. Compare the distance to CUTOFF:\n
            - write down 1 if the distance <= CUTOFF\n
            - write down 0 if the distance > CUTOFF\n

        :param residue: residue as OpenBabel object
        :param residue_atoms: residue's carbon atoms coordinates
        :param ligand_name: ligand_name^pose_number
        :param ligands_lipophilic_coords: list of ligand's detected lipohilic atoms coords
        :type residue: openbabel.OBResidue
        :type residue_atoms: list
        :type ligand_name: str
        :type ligands_lipophilic_coords: list
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0]

    # Flag to iterate over residue's atoms as long as we do not find an atom within CUTOFF distance from ligand
    flag = True

    for rna_atom in residue_atoms:
        if flag:
            for lipophilic in ligands_lipophilic_coords:
                dist = measure_distance(lipophilic, rna_atom)
                if config.MIN_DIST < dist <= config.MAX_LIPOHILIC_DIST:
                    result[-1] = 1

                    if debug:
                        global lipophilic_info
                        lipophilic_info += '***\n'
                        lipophilic_info += '{} - {} \n dist: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, np.round(dist, 4))
                        lipophilic_info += '{}:{}:{}\t{} atom {}\n'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])], debug_dict_ligand[ligand_name][lipophilic][0], ligand_name)

                    if detail:
                        global detail_list
                        ligand_name_detail = ligand_name.split('^')[0]
                        ligand_pose_detail = ligand_name.split('^')[1]

                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][lipophilic][1], 'Lipophilic',
                        debug_dict_ligand[ligand_name][lipophilic][0],
                        lipophilic[0], lipophilic[1],  lipophilic[2],
                        residue.GetName(), residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])],
                        rna_atom[0], rna_atom[1], rna_atom[2],
                        np.round(dist, 4)])

                    if not detail:
                        flag = False
                        break
        else:
            break

    return result

def detect_user_def_interaction(res_name, residue_atoms, ligand_name, ligand_coords, interaction_type, dist_min, dist_max, angle1_min=None, angle1_max=None, angle2_min=None, angle2_max=None):
    """ Calculates user-defined interaction for the SMARTS-defined atoms.\n
        1. If only 2 SMARTS (one for the receptor and one for the ligand) are given:
            - Check nucleic acid SMARTS-defined atom - ligand's SMARTS-defined atom distance
            - Compare the distance to user-defined range:\n
                - write down 1 if the distance within user-defined range\n
                - write down 0 if the distance exceeds user-defined range \n
        2. If 3 SMARTS: two for the receptor and one for the ligand are given:
            - Check nucleic acid SMARTS 2-defined atom - ligand's SMARTS-defined atom distance
            - Compare the distance to user-defined range:\n
                - if the distance within user-defined range: check nucleic acid SMARTS 1-defined atom - nucleic acid SMARTS 2-defined atom - ligand's SMARTS-defined atom angle \n
                    - write down 1 if the angle within user-defined range \n
                    - write down 0 if the angle exceeds user-defined range \n
                - write down 0 if the distance exceeds user-defined range \n
        3. If 3 SMARTS: one for the receptor and two for the ligand are given:
            - Check nucleic acid SMARTS-defined atom - ligand's SMARTS 1-defined atom distance
            - Compare the distance to user-defined range:\n
                - if the distance within user-defined range: check nucleic acid SMARTS-defined atom - ligand's SMARTS 1-defined atom - ligand's SMARTS 2-defined atom angle \n
                    - write down 1 if the angle within user-defined range \n
                    - write down 0 if the angle exceeds user-defined range \n
                - write down 0 if the distance exceeds user-defined range \n
        4. If 4 SMARTS: two for the receptor and two for the ligand are given:
            - Check nucleic acid SMARTS 2-defined atom - ligand's SMARTS 1-defined atom distance
            - Compare the distance to user-defined range:\n
                - if the distance within user-defined range: check nucleic acid SMARTS 1-defined atom - nucleic acid SMARTS 2-defined atom - ligand's SMARTS 1-defined atom angle \n
                    - if the first angle within first user-defined range: check nucleic acid SMARTS 2-defined atom - ligand's SMARTS 1-defined atom - ligand's SMARTS 2-defined atom angle \n
                        - write down 1 if the second angle within second user-defined range \n
                        - write down 0 if the second angle exceeds second user-defined range \n
                    - write down 0 if the first angle exceeds first user-defined range \n
                - write down 0 if the distance exceeds user-defined range \n

        :param res_name: residue name as chain_no:residue_no:residue_name
        :param residue_atoms: residue's SMARTS-defined atoms coordinates as list of tuples
        :param ligand_name: ligand_name^pose_number
        :param ligand_coords: ligand's SMARTS-defined atoms coordinates as list of tuples
        :param interaction_type: name of calculated interaction
        :param dist_min: minimum distance between atoms
        :param dist_max: maximum distance between atoms
        :param angle1_min: minimum first (or only) angle between atoms
        :param angle1_max: maximum first (or only) angle between atoms
        :param angle2_min: minimum second angle between atoms
        :param angle2_max: maximum second angle between atoms
        :type res_name: str
        :type residue_atoms: list
        :type ligand_name: str
        :type ligand_coords: list
        :type interaction_type: str
        :type dist_min: float
        :type dist_max: float
        :type angle1_min: float
        :type angle1_max: float
        :type angle2_min: float
        :type angle2_max: float
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, res_name, 0]
    # Flag to iterate over residue's atoms as long as we do not find an atom within user-defined distance/angle from ligand
    flag = True

    for rna_atom in residue_atoms:
        if flag:
            for ligand_atom in ligand_coords:

                interaction_found = False

                if angle1_min is None:
                    dist = measure_distance(np.array(ligand_atom[1]), np.array(rna_atom[0]))
                elif rna_atom[1] is None: # 3 SMARTS; 1 for receptor & 2 for ligands
                    dist = measure_distance(np.array(ligand_atom[0]), np.array(rna_atom[0]))
                else:
                    if ligand_atom[0] is None: # 3 SMARTS; 2 for receptor & 1 for ligands
                        dist = measure_distance(np.array(ligand_atom[1]), np.array(rna_atom[1]))
                    else: # 4 SMARTS; 2 for receptor & 2 for ligands
                        dist = measure_distance(np.array(ligand_atom[0]), np.array(rna_atom[1]))

                if dist_min < dist <= dist_max:
                    if angle1_min is None:
                        result[-1] = 1
                        interaction_found = True
                    else:

                        if rna_atom[1] is None:
                            ll = vector(np.array(ligand_atom[1]), np.array(ligand_atom[0]))
                            rl = vector(np.array(rna_atom[0]), np.array(ligand_atom[0]))
                            angle1 = calculate_angle(ll, rl)
                        else:
                            rr = vector(np.array(rna_atom[0]), np.array(rna_atom[1]))
                            if ligand_atom[0] is None:
                                rl = vector(np.array(ligand_atom[1]), np.array(rna_atom[1]))
                            else:
                                rl = vector(np.array(ligand_atom[0]), np.array(rna_atom[1]))
                            angle1 = calculate_angle(rr, rl)

                        if angle1_min < angle1 < angle1_max:
                            if angle2_min is None:
                                result[-1] = 1
                                interaction_found = True
                            else:
                                lr = vector(np.array(rna_atom[1]), np.array(ligand_atom[1]))
                                ll = vector(np.array(ligand_atom[0]), np.array(ligand_atom[1]))
                                angle2 = calculate_angle(lr, ll)
                                if angle2_min < angle2 < angle2_max:
                                    result[-1] = 1
                                    interaction_found = True

                    if interaction_found and debug:
                        global new_interactions_info
                        if interaction_type not in new_interactions_info.keys():
                            new_interactions_info[interaction_type] = []
                        if angle1_min is None:
                            new_interactions_info[interaction_type].append('{} - {} \ndist: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, np.round(dist, 4))
                            + '{}:{}:{}\t{} atom {}'.format(res_name.split(':')[1], res_name.split(':')[0], debug_dict_rna[(rna_atom[0][0], rna_atom[0][1], rna_atom[0][2])], debug_dict_ligand[ligand_name][ligand_atom[1]][0], ligand_name))
                        else:
                            if angle2_min is None:
                                new_interactions_info[interaction_type].append('{} - {} \ndist: {} angle: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, np.round(dist, 4), round(angle1, 4))
                                + '{}:{}:{}\t{} atom {}'.format(res_name.split(':')[1], res_name.split(':')[0], debug_dict_rna[(rna_atom[0][0], rna_atom[0][1], rna_atom[0][2])], debug_dict_ligand[ligand_name][ligand_atom[1]][0], ligand_name))
                            else:
                                new_interactions_info[interaction_type].append('{} - {} \ndist: {} angle1: {} angle2: {}\n'.format(filename_RNA.split(sys_sep)[-1], ligand_name, np.round(dist, 4), round(angle1, 4), round(angle2, 4))
                                + '{}:{}:{}\t{} atom {}'.format(res_name.split(':')[1], res_name.split(':')[0], debug_dict_rna[(rna_atom[0][0], rna_atom[0][1], rna_atom[0][2])], debug_dict_ligand[ligand_name][ligand_atom[1]][0], ligand_name))

                    if interaction_found and detail:
                        global detail_list
                        ligand_name_detail = ligand_name.split('^')[0]
                        ligand_pose_detail = ligand_name.split('^')[1]
                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], debug_dict_ligand[ligand_name][ligand_atom[1]][1], interaction_type,
                        debug_dict_ligand[ligand_name][ligand_atom[1]][0],
                        round(ligand_atom[1][0], 4), round(ligand_atom[1][1], 4),  round(ligand_atom[1][2], 4),
                        res_name.split(':')[2], res_name.split(':')[0], res_name.split(':')[1], debug_dict_rna[(rna_atom[0][0], rna_atom[0][1], rna_atom[0][2])],
                        rna_atom[0][0], rna_atom[0][1], rna_atom[0][1],
                        np.round(dist, 4)])

                    if interaction_found and not detail:
                        flag = False
                        break
        else:
            break

    return result

###################################################################################### MAIN ######################################################################################

if __name__ == "__main__":

    welcome_mssg = '# Welcome to fingeRNAt! #'
    columns = shutil.get_terminal_size().columns
    print('')
    print(('#'*len(welcome_mssg)).center(columns))
    print(welcome_mssg.center(columns))
    print(('#'*len(welcome_mssg)).center(columns))
    print('')

    #######################
    #  ARGUMENTS PARSING  #
    #######################

    parser = argparse.ArgumentParser(description = '''Software calculating Structural Interaction Fingerprint (SIFt) in RNA/DNA-ligand complexes.''',
                                     epilog = 'If parameter -o was not passed, fingeRNAt will create outputs/ directory in the working directory and save there SIFt in tsv format.',
                                     add_help = False,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter,
                                     conflict_handler = 'resolve')

    required_arguments = parser.add_argument_group('Required arguments')
    required_arguments.add_argument('-r', help='pass RNA/DNA structure in pdb format', required=True, metavar='RNA/DNA', default=argparse.SUPPRESS)

    optional_arguments = parser.add_argument_group('Optional arguments')
    optional_arguments.add_argument('-l', help='pass ligands file in sdf format', metavar='LIGANDS')
    optional_arguments.add_argument('-f', help='pass fingerprint type, available types are SIMPLE, PBS & FULL', default='FULL', metavar='TYPE')
    optional_arguments.add_argument('-addH', help="pass module name to add hydrogens to ligands' structures, available modules are OpenBabel, RDKit, None", default='OpenBabel', metavar='MODULE')
    optional_arguments.add_argument('-wrapper', help='pass results wrapper types (multiple types possible at once, but have to be comma separated)', metavar='WRAPPER')
    optional_arguments.add_argument('-o', help='pass output path', metavar='NAME')
    optional_arguments.add_argument('-h2o', help='consider water-mediated nucleic acid - ligand interactions', action='store_true')
    optional_arguments.add_argument('-dha', help='consider Donor-Hydrogen-Acceptor angle in hydrogen bonds calculation', action='store_true')
    optional_arguments.add_argument('-custom', help='add user-defined interactions in SMARTS format', metavar='YAML with SMARTS')
    optional_arguments.add_argument('-fingerDISt', help='run directly fingerDISt on the obtained SIFts output', metavar='METRICS')
    optional_arguments.add_argument('-print', help='print found interactions on screen', action='store_true')
    optional_arguments.add_argument('-detail', help='generate an additional file with detailed data on detected interactions (used for PyMOL visualization)', action='store_true')
    optional_arguments.add_argument('-verbose', help='provides additional details about calculations performed at the given moment', action='store_true')
    optional_arguments.add_argument('-debug', help='enter debug mode', action='store_true')
    optional_arguments.add_argument('-h', action = 'help', help = 'show this help message and exit')
    optional_arguments.add_argument('--help', '-h', action = 'help', help = 'show this help message and exit')

    args = vars(parser.parse_args())
    sys_sep = os.sep

    filename_RNA = args['r']
    extension_structure = ".".join(filename_RNA.split(sys_sep)[-1].split('.')[1:])
    if extension_structure != 'pdb':
        raise Exception('RNA/DNA structure has to be in pdb format!')

    filename_ligand=args['l']
    if filename_ligand:
        extension_ligand = ".".join(filename_ligand.split(sys_sep)[-1].split('.')[1:])
        if extension_ligand != 'sdf':
            raise Exception("Ligands' structures have to be in sdf format!")

    fingerprint = args['f']
    if not filename_ligand:
        if fingerprint not in ['SIMPLE', 'PBS']:
            raise Exception("If not specyfing ligands path, only fingerprint SIMPLE or PBS can be calculated!")

    output = args['o']
    consider_dha = args['dha']
    consider_H2O = args['h2o']
    how_addH = args['addH']
    new_interactions = args['custom']
    run_fingerDISt = args['fingerDISt']

    if new_interactions and fingerprint != 'FULL':
        raise Exception("User-defined interactions can only be calculated for fingerprint FULL!")

    if how_addH not in ['OpenBabel', 'RDKit', 'None']: raise Exception('Unknown module to add hydrogens!')

    try:
        wrapper = args['wrapper'].split(',')
    except AttributeError:
        wrapper = None

    print_flag = args['print']
    detail = args['detail']
    verbose = args['verbose']
    debug = args['debug']

    #########################
    #  FINGERPRINT CALLING  #
    #########################

    FUNCTIONS = {'SIMPLE' : 1, 'PBS' : 3, 'FULL' : 12}
    WRAPPERS = {'ACUG' : 4, 'PuPy' : 2, 'Counter' : 1 }
    ANALYSIS_NAME = []
    IONS_DIST_CUTOFFS = {'MG' : config.MAX_MAGNESIUM_DIST, 'K' : config.MAX_POTASSIUM_DIST,
                         'NA' : config.MAX_SODIUM_DIST}

    if fingerprint not in FUNCTIONS.keys(): raise Exception('Unknown fingerprint')
    if wrapper:
        for w in wrapper:
            if w not in WRAPPERS.keys(): raise Exception('Unknown wrapper(s)')

    ANALYSIS_NAME.append(fingerprint)
    if wrapper:
        ANALYSIS_NAME.extend(wrapper)

    # Create dictionary of all residues with their calculated interactions
    RESULTS = {}

    # Parse nucleic acid using OpenBabel
    structure = next(pybel.readfile(extension_structure,filename_RNA))

    RNA_residues_objects = []
    Ions_objects  = []
    HOH_objects  = []

    for residue in openbabel.OBResidueIter(structure.OBMol):
        if residue.GetNumAtoms() > 3:
            RNA_residues_objects.append(residue)
        elif residue.GetName() == 'HOH' and consider_H2O:
            HOH_objects.append(residue)
        elif residue.GetNumAtoms() == 1 and residue.GetName() != 'HOH':
            Ions_objects.append(residue)
        else:
            pass

    # Get number of nucleic acid's residues
    RNA_LENGTH = len(RNA_residues_objects )
    # Create empty list of all nucleic acid's residues
    RNA_residues = []
    # Create empty list of all nucleic acid's residues
    RNA_nucleotides = []

    if not filename_ligand:

        if verbose:
            mssg = '# Calculating nucleic acid - ions interactions type {} #'.format(fingerprint)
            print('#'*len(mssg))
            print(mssg)
            print('#'*len(mssg))

        Inorganic_ions_dict = {}

        for inorganic_ion in Ions_objects:
            ion_id = inorganic_ion.GetName() + ':' + str(inorganic_ion.GetNum()) + ':' + inorganic_ion.GetChain()

            # Fill the RESULTS dictionary of keys - ions ids and values - lists of 0
            RESULTS[ion_id] = [0] * RNA_LENGTH  * FUNCTIONS[fingerprint]

            for atom in openbabel.OBResidueAtomIter(inorganic_ion):
                Inorganic_ions_dict[ion_id] = np.array([[atom.GetX(), atom.GetY(), atom.GetZ()]])

        ####################################################################

        if debug or detail:
            # Create dicts with atoms' coords as keys and their indices as values
            debug_dict_rna = rna_coords_atom_index_dict(structure)
            debug_dict_ligand = {}
            for inorganic_ion in Ions_objects:
                ion_id = inorganic_ion.GetName() + ':' + str(inorganic_ion.GetNum()) + ':' + inorganic_ion.GetChain()
                debug_dict_ligand[ion_id] = {}
                for atom in openbabel.OBResidueAtomIter(inorganic_ion):
                    debug_dict_ligand[ion_id][(atom.GetX(), atom.GetY(), atom.GetZ())] = (str(inorganic_ion.GetNum()) + ':' + inorganic_ion.GetChain(), 0)
            detail_list = []

            if debug:
                debug_mssg = '# Entering DEBUG MODE of  nucleic acid - ions interactions type {} #'.format(fingerprint)
                print(('#'*len(debug_mssg)).center(columns))
                print((debug_mssg).center(columns))
                print(('#'*len(debug_mssg)).center(columns))
                print()

        ####################################################################

        for residue in RNA_residues_objects: # Loop over all nucleic acid residues
            RNA_residues.append(str(residue.GetNum())+ ':' + str(residue.GetChain()))
            RNA_nucleotides.append(str(residue.GetName()))

            for inorganic_ion in Inorganic_ions_dict.keys():
                inorganic_ion_name = inorganic_ion.split(':')[0]
                if inorganic_ion_name in IONS_DIST_CUTOFFS.keys():
                    ion_cutoff = IONS_DIST_CUTOFFS[inorganic_ion_name]
                else:
                    ion_cutoff = config.MAX_OTHER_ION_DIST

                if fingerprint == 'SIMPLE':
                    result = calculate_SIMPLE(residue, inorganic_ion, Inorganic_ions_dict[inorganic_ion], Inorganic_ions_dict[inorganic_ion], ion_cutoff)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0)
                else:
                    result = calculate_PBS(residue, inorganic_ion, Inorganic_ions_dict[inorganic_ion], Inorganic_ions_dict[inorganic_ion], ion_cutoff)
                    if sum(result[-3:]) != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0, 'PBS') # Assign each of 3 (P/B/S) residue-ligand interaction
    else:

        if verbose:
            mssg = '# Calculating fingerprint type {} #'.format(fingerprint)
            print('#'*len(mssg))
            print(mssg)
            print('#'*len(mssg))

        if fingerprint == 'SIMPLE' or fingerprint == 'PBS':

            # Read all ligands
            ligands_mols = list(pybel.readfile(extension_ligand, filename_ligand))

            ############################ For debug mode #########################

            if debug or detail:
                # Create dicts with atoms' coords as keys and their indices as values
                debug_dict_rna = rna_coords_atom_index_dict(structure)
                debug_dict_ligand = ligands_coords_atom_index_dict(ligands_mols)
                detail_list = []
                if debug:
                    debug_mssg = '# Entering DEBUG MODE of {} #'.format(fingerprint)
                    print(('#'*len(debug_mssg)).center(columns))
                    print((debug_mssg).center(columns))
                    print(('#'*len(debug_mssg)).center(columns))
                    print()

            ####################################################################

            # Create ligands all atoms (except hydrogens) dictionary
            ligands_all_atoms = find_ligands_all_atoms(ligands_mols)

            # Fill the RESULTS dictionary of keys - ligand ids and values - lists of 0
            for ligand_name in ligands_all_atoms.keys():
                RESULTS[ligand_name] = [0] * RNA_LENGTH * FUNCTIONS[fingerprint]

            for residue in RNA_residues_objects: # Loop over all nucleic acid residues
                RNA_residues.append(str(residue.GetNum())+ ':' + str(residue.GetChain()))
                RNA_nucleotides.append(str(residue.GetName()))

                for ligand_name, ligand_atoms in ligands_all_atoms.items():

                    centroid_ligand = centroid(ligand_atoms)

                    if fingerprint == 'SIMPLE':
                        result = calculate_SIMPLE(residue, ligand_name, ligand_atoms, centroid_ligand, config.CUT_OFF_SIMPLE)
                        if result[-1] != 0:
                            assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0)
                    else:
                        result = calculate_PBS(residue, ligand_name, ligand_atoms, centroid_ligand, config.CUT_OFF_SIMPLE)
                        if sum(result[-3:]) != 0:
                            assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0, 'PBS') # Assign each of 3 (P/B/S) residue-ligand interaction

        else: # fingerprint == 'FULL'

            # Add missing hydrogens according to -addH parameter
            if how_addH == 'OpenBabel':
                # Read all ligands
                ligands_mols = list(pybel.readfile(extension_ligand, filename_ligand))

                save_ligands_addedH = ''
                for i in range(len(ligands_mols)):
                    ligands_mols[i].OBMol.AddHydrogens()
                    save_ligands_addedH += ligands_mols[i].write('sdf')

                f=open(filename_ligand[:-4] + '_OB_addedH.sdf', 'w')
                f.write(save_ligands_addedH)
                f.close()
            elif how_addH == 'RDKit':
                addHwithRDKit(filename_ligand, filename_ligand[:-4] + '_RDKit_addedH.sdf')
                # Read all modified by RDKit ligands
                ligands_mols = list(pybel.readfile(extension_ligand, filename_ligand[:-4] + '_RDKit_addedH.sdf'))
            else:
                # Read all ligands
                ligands_mols = list(pybel.readfile(extension_ligand, filename_ligand))

            # Create dictionary of positively charged ions with their coords
            Inorganic_ions_dict = {}
            for inorganic_ion in Ions_objects:
                if inorganic_ion.GetName() in config.POS_CHARGED_IONS:
                    for atom in openbabel.OBResidueAtomIter(inorganic_ion):
                        ion_id = inorganic_ion.GetName() + ':' + str(inorganic_ion.GetNum()) + ':' +  str(inorganic_ion.GetChain())
                        Inorganic_ions_dict[ion_id] = np.array([[atom.GetX(), atom.GetY(), atom.GetZ()]])

            # Create dictionary of water molecules with their coords
            if consider_H2O:
                Water_dict = {}
                for water in HOH_objects:
                    for atom in openbabel.OBResidueAtomIter(water):
                        if atom.GetAtomicNum() == config.OXYGEN_NUM:
                            water_id = 'HOH' + ':' + str(water.GetNum()) + ':' +  str(water.GetChain())
                            Water_dict[water_id] = np.array([[atom.GetX(), atom.GetY(), atom.GetZ()]])

            ############################ For debug mode #########################

            if debug or detail:

                # Create dicts with atoms' coords as keys and their indices as values
                debug_dict_rna = rna_coords_atom_index_dict(structure)
                debug_dict_ligand = ligands_coords_atom_index_dict(ligands_mols)
                detail_list = []

                if debug:
                    debug_mssg = '# Entering DEBUG MODE of {} #'.format(fingerprint)
                    print(('#'*len(debug_mssg)).center(columns))
                    print((debug_mssg).center(columns))
                    print(('#'*len(debug_mssg)).center(columns))
                    print()

                    arom_ring_ligands_info = {}
                    RNA_HB_acc_don_info = {}
                    RNA_anion_info = {}
                    arom_RNA_ligands_info = {}
                    HB_RNA_acc_info = ''
                    HB_RNA_donor_info = ''
                    user_def_receptor_atoms = {}
                    user_def_ligands_atoms = {}
                    HAL_info = ''
                    Cation_Anion_info = ''
                    Anion_Pi_info = ''
                    Pi_Cation_info = ''
                    Pi_Anion_info = ''
                    Sandwich_Displaced_info = ''
                    T_shaped_info = ''
                    ion_mediated_info = ''
                    water_mediated_info = ''
                    lipophilic_info = ''
                    new_interactions_info = {}

            ####################################################################

            # Create ligands' all atoms dictionary
            ligands_all_atoms = find_ligands_all_atoms(ligands_mols)
            # Create ligands acceptors' & donors' dictionary
            ligands_hba_hbd = find_ligands_HBA_HBD(ligands_mols, verbose)
            # Create ligands halogens' donors dictionary
            ligands_HAL = find_ligands_HAL_don(ligands_mols, verbose)
            # Create ligands cations' & anions' dictionary
            ligands_CA = find_ligands_CA(ligands_mols, verbose)
            # Create dictionary of ligands in contacts with ions
            ligands_ions = find_ligands_ions(ligands_mols, Inorganic_ions_dict, verbose)
            # Create dictionary of ligands in contacts with water molecules
            if consider_H2O: ligands_water = find_ligands_water(ligands_hba_hbd, Water_dict, verbose)
            # Create dictionary of hydrophobic regions of ligand
            ligands_lipophilic = find_ligands_lipophilic(ligands_mols, verbose)
            # Find all ligands' rings
            ligands_aromatic_rings = findAromaticRingsWithRDKit(filename_ligand)
            # Find all RNA rings
            rings_RNA = find_RNA_rings(structure)
            # Get dictionary with atoms from user-defined interactions
            if new_interactions:
                additionalInteractions = parseYaml(new_interactions)
                for add_int in additionalInteractions: # Check if there are two SMARTS for both receptor & ligand when Angle2 is defined
                    if 'Angle2' in additionalInteractions[add_int].keys():
                        if len(additionalInteractions[add_int]['Receptor_SMARTS']) != 2 or len(additionalInteractions[add_int]['Ligand_SMARTS']) != 2:
                            raise Exception('Improperly defined {} in the yaml file!'.format(add_int))

                FUNCTIONS[fingerprint] += len(additionalInteractions.keys())
                dict_rna_coords_to_res_ids = rna_coords_residue_index_dict(structure)
                user_defined_all_atoms = find_atoms_from_SMARTS(structure, ligands_mols, additionalInteractions, dict_rna_coords_to_res_ids, verbose)
                print(user_defined_all_atoms)
                if debug:
                    for k in user_defined_all_atoms.keys():
                        user_def_receptor_atoms[k] = {}
                        user_def_ligands_atoms[k] = {}

                        for residue in user_defined_all_atoms[k]['Receptor'].keys():
                            chain = residue.split(':')[1]
                            if chain not in user_def_receptor_atoms[k].keys():
                                user_def_receptor_atoms[k][chain] = {}
                            res_no = residue.split(':')[0]
                            user_def_receptor_atoms[k][chain][res_no] = []
                            for el in user_defined_all_atoms[k]['Receptor'][residue]:
                                if el[1] is not None:
                                    user_def_receptor_atoms[k][chain][res_no].append((debug_dict_rna[el[0]], debug_dict_rna[el[1]]))
                                else:
                                    user_def_receptor_atoms[k][chain][res_no].append(debug_dict_rna[el[0]])

                        for lig_name in user_defined_all_atoms[k]['Ligands'].keys():
                            if lig_name not in user_def_ligands_atoms[k].keys():
                                user_def_ligands_atoms[k][lig_name] = []
                            for l_a in user_defined_all_atoms[k]['Ligands'][lig_name]:
                                if l_a[0] is not None:
                                    user_def_ligands_atoms[k][lig_name].append((debug_dict_ligand[lig_name][l_a[1]], debug_dict_ligand[lig_name][l_a[0]]))
                                else:
                                    user_def_ligands_atoms[k][lig_name].append(debug_dict_ligand[lig_name][l_a[1]])

            # Fill the RESULTS dictionary of keys - ligand ids and values - lists of 0
            for ligand_name in ligands_all_atoms.keys():
                RESULTS[ligand_name] = [0] * RNA_LENGTH * FUNCTIONS[fingerprint]

            # Loop over all nucleic acid residue to calculate all except Pi-interactions
            for residue in RNA_residues_objects:

                res_name = str(residue.GetNum())+ ':' + str(residue.GetChain())
                RNA_residues.append(res_name)
                RNA_nucleotides.append(str(residue.GetName()))
                acceptors_RNA, donors_RNA = find_RNA_HB_HAL_acc_don(residue)
                anions_RNA = find_RNA_anions(residue)

                if debug:
                    if residue.GetChain() not in RNA_HB_acc_don_info.keys():
                        RNA_HB_acc_don_info[residue.GetChain()] = {}
                    if residue.GetNum() not in RNA_HB_acc_don_info[residue.GetChain()].keys():
                        RNA_HB_acc_don_info[residue.GetChain()][residue.GetNum()] =  [[],[]]
                    for a in acceptors_RNA:
                        RNA_HB_acc_don_info[residue.GetChain()][residue.GetNum()][0].append(debug_dict_rna[(a[0].GetX(), a[0].GetY(), a[0].GetZ())])
                    for d in donors_RNA:
                        RNA_HB_acc_don_info[residue.GetChain()][residue.GetNum()][1].append(debug_dict_rna[(d[0].GetX(), d[0].GetY(), d[0].GetZ())])
                    if residue.GetChain() not in RNA_anion_info.keys():
                        RNA_anion_info[residue.GetChain()] = {}
                    if residue.GetNum() not in RNA_anion_info[residue.GetChain()].keys():
                        RNA_anion_info[residue.GetChain()][residue.GetNum()] =  []
                    for an in anions_RNA:
                        RNA_anion_info[residue.GetChain()][residue.GetNum()].append(debug_dict_rna[(an.GetX(), an.GetY(), an.GetZ())])

                residue_atoms_coords = []
                residue_carbon_atoms_coords = []
                residue_oxygen_nitrogen_atoms_coords = []
                for atom in openbabel.OBResidueAtomIter(residue):
                    if atom.GetAtomicNum() != config.HYDROGEN_NUM:
                        residue_atoms_coords.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))
                    if atom.GetAtomicNum() == config.CARBON_NUM:
                        residue_carbon_atoms_coords.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))
                    if atom.GetAtomicNum() == config.OXYGEN_NUM or atom.GetAtomicNum() == config.NITROGEN_NUM:
                        residue_oxygen_nitrogen_atoms_coords.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))

                for ligand_name_HB, ligand_values_HB in ligands_hba_hbd.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_HB])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_HB(residue, acceptors_RNA, donors_RNA, ligand_name_HB, ligand_values_HB, consider_dha)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0)

                for ligand_name_HAL, ligand_values_HAL in ligands_HAL.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_HAL])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_HAL(residue, acceptors_RNA, ligand_name_HAL, ligand_values_HAL)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 1)

                for ligand_name_CA, ligand_values_CA in ligands_CA.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_CA])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_CATION_ANION(residue, anions_RNA, ligand_name_CA, ligand_values_CA[0])
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 2)

                for ligand_name_ring, ligand_values_ring in ligands_aromatic_rings.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_ring])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_ANION_PI(residue, anions_RNA, ligand_name_ring, ligand_values_ring)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 4)

                for ligand_name, ions in ligands_ions.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_ION_MEDIATED(residue, residue_oxygen_nitrogen_atoms_coords, ligand_name, ions, Inorganic_ions_dict)
                    if sum(result[-4:]) != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 6, 'FULL')

                if consider_H2O:
                    for ligand_name, water_molecules in ligands_water.items():
                        if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                            continue

                        result = calculate_WATER_MEDIATED(residue, acceptors_RNA, donors_RNA, ligand_name, water_molecules, Water_dict)
                        if result[-1] != 0:
                            assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 10)

                for ligand_name_lipophilic, ligand_values_lipophilic in ligands_lipophilic.items():
                    if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_lipophilic])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                        continue

                    result = calculate_lipophilic_interactions(residue, residue_carbon_atoms_coords, ligand_name_lipophilic, ligand_values_lipophilic)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 11)

                if new_interactions:
                    c = 0
                    for new_interaction_type in user_defined_all_atoms.keys():
                        c += 1
                        if res_name in user_defined_all_atoms[new_interaction_type]['Receptor'].keys() and user_defined_all_atoms[new_interaction_type]['Ligands'].keys():
                            distance_min = additionalInteractions[new_interaction_type]['Distance']['min']
                            distance_max = additionalInteractions[new_interaction_type]['Distance']['max']
                            if 'Angle1' not in additionalInteractions[new_interaction_type].keys():
                                for ligand_name in user_defined_all_atoms[new_interaction_type]['Ligands'].keys():
                                    result = detect_user_def_interaction(res_name + ':{}'.format(residue.GetName()), user_defined_all_atoms[new_interaction_type]['Receptor'][res_name], ligand_name, user_defined_all_atoms[new_interaction_type]['Ligands'][ligand_name], new_interaction_type, distance_min, distance_max)
                                    if result[-1] != 0:
                                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 11+c)
                            else:
                                angle1_min = additionalInteractions[new_interaction_type]['Angle1']['min']
                                angle1_max = additionalInteractions[new_interaction_type]['Angle1']['max']

                                if 'Angle2' in additionalInteractions[new_interaction_type].keys():
                                    angle2_min = additionalInteractions[new_interaction_type]['Angle2']['min']
                                    angle2_max = additionalInteractions[new_interaction_type]['Angle2']['max']
                                else:
                                    angle2_min = None
                                    angle2_max = None

                                for ligand_name in user_defined_all_atoms[new_interaction_type]['Ligands'].keys():
                                    result = detect_user_def_interaction(res_name + ':{}'.format(residue.GetName()), user_defined_all_atoms[new_interaction_type]['Receptor'][res_name], ligand_name, user_defined_all_atoms[new_interaction_type]['Ligands'][ligand_name], new_interaction_type, distance_min, distance_max, angle1_min, angle1_max, angle2_min, angle2_max)
                                    if result[-1] != 0:
                                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 11+c)

            PI_INTERACTIONS = calculate_PI_INTERACTIONS(rings_RNA, structure.atoms, ligands_CA, ligands_aromatic_rings) # Calculate Pi-cation, Pi-anion & Pi-stacking interactions

            for i in range(3):

                result = PI_INTERACTIONS[i] # i=0 Pi-cation interaction; i=1 Pi-anion interaction; i=2 Pi-stacking interaction
                for res in result:
                    if res[-1] != 0: # We assign only 1
                        assign_interactions_results(res, RESULTS, RNA_LENGTH, RNA_residues.index(res[1]), FUNCTIONS[fingerprint], i+3)
            if debug:
                if not consider_H2O: ligands_water = None
                print_debug_info(ligands_hba_hbd, ligands_HAL, ligands_CA, ligands_ions, ligands_water, ligands_lipophilic, arom_ring_ligands_info,
                user_def_ligands_atoms, debug_dict_ligand, RNA_HB_acc_don_info, RNA_anion_info, arom_RNA_ligands_info, user_def_receptor_atoms,
                HB_RNA_acc_info, HB_RNA_donor_info,HAL_info, Cation_Anion_info, Pi_Cation_info, Pi_Anion_info, Anion_Pi_info,
                Sandwich_Displaced_info, T_shaped_info, ion_mediated_info, water_mediated_info, lipophilic_info, new_interactions_info, columns)

    # Wrap results if wrapper was passed
    if wrapper:

        WRAP_RESULTS = {}
        for w in wrapper:
            if verbose: print('Wrapping fingerprint type {} results to {} wrapper'.format(fingerprint, w))
            WRAP_RESULTS[w] = wrap_results(w, RESULTS, RNA_nucleotides, FUNCTIONS[fingerprint], WRAPPERS[w])

    # Create dataframe
    columns = {'SIMPLE': ['SIMPLE'],
               'PBS': ['P', 'B', 'S'],
               'FULL': ['HB', 'HAL', 'CA', 'Pi_Cation', 'Pi_Anion', 'Pi_Stacking',
               'Mg_mediated', 'K_mediated', 'Na_mediated', 'Other_mediated',
               'Water_mediated', 'Lipophilic']
              }

    if new_interactions:
        columns['FULL'] = columns['FULL'] + [key for key in additionalInteractions.keys()]

    is_structure_RNA = check_if_RNA(RNA_nucleotides)

    if is_structure_RNA:
        nucleotides_letters = ['A','C','U','G']
    else:
        nucleotides_letters = ['A','C','T','G']

    if verbose: print('Calculations completed, saving the results...')

    if detail:
        detail_already_saved = False #save DETAIL file only once as it is identical for fingerprint and all its wrappers

    for analysis in ANALYSIS_NAME:

        if analysis == 'ACUG':
            DF_COLUMNS = [filename_RNA.split(sys_sep)[-1] + '#' + res + '#' + fing_type for res in nucleotides_letters for fing_type in columns[fingerprint]]
        elif analysis == 'PuPy':
            DF_COLUMNS = [filename_RNA.split(sys_sep)[-1] + '#' + res + '#' + fing_type for res in ['Purines','Pyrimidynes'] for fing_type in columns[fingerprint]]
        elif analysis == 'Counter':
             DF_COLUMNS = [filename_RNA.split(sys_sep)[-1] + '#' + fing_type for fing_type in columns[fingerprint]]
        else:
             DF_COLUMNS = [filename_RNA.split(sys_sep)[-1] + '#' + res + '#' + fing_type for res in RNA_residues for fing_type in columns[fingerprint]]

        DF_INDEXES = list(RESULTS.keys())
        ALL_FINGERPRINTS_DF = pd.DataFrame(index = DF_INDEXES, columns = DF_COLUMNS)
        ALL_FINGERPRINTS_DF.index.name = 'Ligand_name'

        if analysis in FUNCTIONS.keys():
            for index in DF_INDEXES:
                ALL_FINGERPRINTS_DF.loc[index] = RESULTS[index]
        else:
            for index in DF_INDEXES:
                ALL_FINGERPRINTS_DF.loc[index] = WRAP_RESULTS[analysis][index]

        if detail:
            if filename_ligand:
                detail_list.sort(key=lambda x: x[2])
            detail_columns = ['Ligand_name', 'Ligand_pose', 'Ligand_occurrence_in_sdf', 'Interaction',
                              'Ligand_Atom', 'Ligand_X', 'Ligand_Y', 'Ligand_Z',
                              'Receptor_Residue_Name', 'Receptor_Number', 'Receptor_Chain', 'Receptor_Atom',
                              'Receptor_X', 'Receptor_Y', 'Receptor_Z', 'Distance']
            detail_df = pd.DataFrame(detail_list, columns=detail_columns)

        if not consider_H2O:
            all_columns = ALL_FINGERPRINTS_DF.columns
            for col_name in all_columns:
                if 'Water' in col_name:
                    ALL_FINGERPRINTS_DF[col_name] = None

        # Save output as tsv
        if output:
            if not filename_ligand: filename_ligand='IONS'
            output_proper = output
            save_name = filename_RNA.split(sys_sep)[-1] + '_' + filename_ligand.split(sys_sep)[-1] + '_' + fingerprint
            if analysis in FUNCTIONS.keys():
                if output[-1] == sys_sep:
                    output_proper += save_name
                else:
                    output_proper += sys_sep + save_name
            else:
                if output[-1] == sys_sep:
                    output_proper += save_name + '_' + analysis
                else:
                    output_proper += sys_sep + save_name + '_' + analysis

            ALL_FINGERPRINTS_DF.to_csv('%s.tsv' %output_proper, sep='\t')

            if detail and not detail_already_saved:
                is_sep = True
                if sys_sep in output:
                    detail_save = output_proper.split(sys_sep)
                else:
                    detail_save = output
                    is_sep = False

                if is_sep:
                    detail_save[-1] = 'DETAIL_' + detail_save[-1]
                    detail_save = sys_sep.join(detail_save)
                else:
                    detail_save = os.path.join(detail_save, 'DETAIL_%s_%s_%s' %(filename_RNA.split(sys_sep)[-1], filename_ligand.split(sys_sep)[-1], fingerprint))
                detail_df.to_csv('%s.tsv' %detail_save, sep='\t')
                detail_already_saved = True

        else:
            if not filename_ligand: filename_ligand = 'IONS'
            if not os.path.exists('outputs'): os.makedirs('outputs')
            save_name = filename_RNA.split(sys_sep)[-1] + '_' + filename_ligand.split(sys_sep)[-1] + '_' + fingerprint
            if analysis in FUNCTIONS.keys():
                ALL_FINGERPRINTS_DF.to_csv(os.path.join('outputs', '%s.tsv' %(save_name)), sep='\t')
                if detail and not detail_already_saved:
                      detail_save = os.path.join('outputs', 'DETAIL_%s' %save_name)
                      detail_df.to_csv('%s.tsv' %detail_save, sep='\t' )
                      detail_already_saved = True
            else:
                ALL_FINGERPRINTS_DF.to_csv(os.path.join('outputs', '%s_%s.tsv' %(save_name, analysis)), sep='\t')

    # Print found interactions on screen
        if print_flag:
            interact_names = {'P':'Phosphate contact', 'B': 'Base contact', 'S' : 'Sugar contact', 'SIMPLE' : 'contact',
                              'HB': 'Hydrogen Bonds', 'HAL': 'Halogen Bonds', 'CA' : 'Cation-Anion', 'Pi_Cation' : 'Pi-Cation',
                              'Pi_Anion' : 'Pi-Anion', 'Pi_Stacking' : 'Pi-Stacking',
                              'Mg_mediated' : 'Magnesium ion-mediated', 'K_mediated' : 'Potassium ion-mediated', 'Na_mediated' : 'Sodium ion-mediated',
                              'Other_mediated' : 'Other ion-mediated', 'Water_mediated' : 'Water-mediated', 'Lipophilic' : 'Lipophilic'}
            if new_interactions:
                for key in additionalInteractions.keys():
                    interact_names[key] = key

            for index, row in ALL_FINGERPRINTS_DF.iterrows():
                print('# {} - {} #'.format(filename_RNA.split(sys_sep)[-1], index))
                print()
                for el in range(len(row)):
                    if row[el] is not None and row[el] > 0:
                        s = DF_COLUMNS[el].split('#')[1:]
                        print('{}\t{}\t{}'.format(s[0],interact_names[s[1]], row[el]))
                print()
            print_flag = False

        if run_fingerDISt:
            if output:
                saved_output_fingerprint = output
                if output[-1] == sys_sep:
                    saved_output_fingerprint += save_name + '.tsv'
                else:
                    saved_output_fingerprint += sys_sep + save_name + '.tsv'
            else:
                saved_output_fingerprint = os.path.join('outputs', '%s.tsv' %(save_name))

            if output:
                command = 'python %s/fingerDISt.py -i %s -m %s -o %s' %(os.path.dirname(os.path.realpath(__file__)), saved_output_fingerprint, run_fingerDISt, output)
            else:
                command = 'python %s/fingerDISt.py -i %s -m %s' %(os.path.dirname(os.path.realpath(__file__)), saved_output_fingerprint, run_fingerDISt)

        print('{} results saved successfully!'.format(analysis))

        if run_fingerDISt:
            import subprocess
            if subprocess.call(command, shell = True):
                print('fingerDISt result saved successfully!')
