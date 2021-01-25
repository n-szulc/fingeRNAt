#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
fingeRNAt is a software to calculate Structural Interaction Fingerprints (SIFs) in
nucleic acids - ligands complexes.

Authors:
Natalia A. Szulc, nszulc@iimcb.gov.pl
Filip Stefaniak, fstefaniak@genesilico.pl

If you use this software, please cite:
Natalia A. Szulc, Zuzanna Mackiewicz, Janusz M. Bujnicki, Filip Stefaniak
[in preparation]

Requires Python 3.5 - 3.8
'''

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os
import shutil
from openbabel import openbabel
from openbabel import pybel
from tqdm import tqdm

# Own modules
import config
from preprocessing import measure_distance, check_distance, vector, calculate_angle, modify_HB_result_list
from preprocessing import calculate_planar, centroid, get_ligand_name_pose, projection, find_ligands_all_atoms
from preprocessing import assign_interactions_results, wrap_results, find_ligands_HBA_HBD, find_ligands_HAL_don
from preprocessing import find_ligands_CA, find_RNA_rings, find_RNA_HB_HAL_acc_don, find_RNA_anions, check_if_RNA
from preprocessing import findAromaticRingsWithRDKit, rna_coords_atom_index_dict, ligands_coords_atom_index_dict, print_debug_info


##################################################
#  FUNCTIONS CALCULATING MOLECULAR INTERACTIONS  #
##################################################

def calculate_SIMPLE(residue, ligand_name, ligand_atoms, centroid_ligand):
    """Calculates SIMPLE interaction between residue - ligand pair:
            1. Check nucleic acid residue - ligand distance
            2. Compare the distance to CUTOFF:
                - write down 1 if the distance <= CUTOFF
                - write down 0 if the distance > CUTOFF

        :param residue: residue as OpenBabel object
        :param ligand_name: ligand_name^pose_number
        :param ligand_atoms: coordinates of all ligand's atoms (except hydrogens)
        :type residue: openbabel.OBResidue
        :type ligand_name: str
        :type ligand_atoms: list
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
                if check_distance(ligand_atom, rna_atom, config.CUT_OFF_SIMPLE):
                        result[-1]=1 # Condition met; write down 1

                        if debug:
                            print('### {} - {} first below cutoff dist: {} ###'.format(filename_RNA.split('/')[-1], ligand_name, np.round(np.linalg.norm(ligand_atom - rna_atom), 4)))
                            print('    between\t{}:{}:{}\t {} atom {}'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])], debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])], ligand_name))

                        if detail:
                            global detail_list
                            detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'SIMPLE',
                            debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])],
                            ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                            residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom[0], rna_atom[1], rna_atom[2])],
                            rna_atom[0], rna_atom[1], rna_atom[2],
                            np.round(np.linalg.norm(ligand_atom - rna_atom), 4)])

                        if not detail:
                            flag=False
                            break
        else:
            break

    return result


def calculate_PBS(residue, ligand_name, ligand_atoms, centroid_ligand):
    """Calculates PBS interaction between residue - ligand pair:
            1. Divide nucleic acid's residue into 3 groups Phosphate/Base/Sugar (P/B/S)
            2. Check each group - ligand distance
            3. Compare the distance to CUTOFF:
                - write down 1 if the distance <= CUTOFF
                - write down 0 if the distance > CUTOFF

        :param residue: residue as OpenBabel object
        :param ligand_name: ligand_name^pose_number
        :param ligand_atoms: coordinates of all ligand's atoms (except hydrogens)
        :type residue: openbabel.OBResidue
        :type ligand_name: str
        :type ligand_atoms: list
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
                        if check_distance(ligand_atom, rna_atom_coords, config.CUT_OFF_SIMPLE):
                                result[g-3]= 1 # Condition met; write down 1

                                if debug:
                                    print('### {} - {} first below cutoff {} GROUP dist: {} ###'.format(filename_RNA.split('/')[-1], ligand_name, config.WHICH_GROUP[g], np.round(np.linalg.norm(ligand_atom - rna_atom_coords), 4)))
                                    print('    between\t{}:{}:{}\t {} atom {}'.format(residue.GetChain(), residue.GetNum(), debug_dict_rna[(rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2])], debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])], ligand_name))

                                if detail:
                                    global detail_list
                                    detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'PBS',
                                    debug_dict_ligand[ligand_name][(ligand_atom[0],ligand_atom[1],ligand_atom[2])],
                                    ligand_atom[0], ligand_atom[1],  ligand_atom[2],
                                    residue.GetNum(), residue.GetChain(), debug_dict_rna[(rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2])],
                                    rna_atom_coords[0], rna_atom_coords[1], rna_atom_coords[2],
                                    np.round(np.linalg.norm(ligand_atom - rna_atom_coords), 4)])

                                if not detail:
                                    flags[g]=False
                                    break

        if not atom_group_type:
            raise Exception('Unknown atom type %s' %rna_atom_name)

    return result


def calculate_HB(residue, acceptors_RNA, donors_RNA, ligand_name, ligand_donors_acceptors, precision):
    """Calculates hydrogen bond between residue - ligand pair.
        Simplified graphical representation:\n
        A \***** H --- D\n
            where:
            A  - hydrogen bond acceptor\n
            H  - hydrogen\n
            D  - hydrogen bond acceptor\n
            \* - hydrogen bond\n
        Geometric Rules are:
            1. D-A distance < 3.9 A
            2. 100 < D-H-A angle < 260\n
            If fingerprint's type is FULL:
                - simply checks if there is any hydrogen bonding in residue - ligand pair
            If fingerprint's type is XP:
                - calculates total number of hydrogen bondings in residue - ligand pair
                - assigns type for each hydrogen bonding (strong/moderate/weak), depending on the D-A distance, and calculates number of each type's occurence
                    - strong hydrogen bond:   2.2 A < D-A < 2.5 A
                    - moderate hydrogen bond: 2.5 A < D-A < 3.5 A
                    - weak hydrogen bond:     3.5 A < D-A < 3.9 A

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors
        :param donors_RNA: all tuples (D, H) of the residue
        :param ligand_name: ligand_name^pose_number
        :param ligand_donors_acceptors: [[sublist of acceptors coords], [sublist of tuples (D, H) coords]]
        :param precision: fingerprint type
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type donors_RNA: list
        :type ligand_name: str
        :type ligand_donors_acceptors: list
        :type precision: str
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    if precision == 'XP':
        result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0, 0, 0, 0]
    else:
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

        RNA_acc = RNA_acceptor_set[0] # We do not need coords of acceptor's neighbours

        if searching_flag:
            RNA_acc_coords = np.array([RNA_acc.GetX(), RNA_acc.GetY(), RNA_acc.GetZ()])

            for donor in ligand_donors_coords:
                dist = measure_distance(donor[0], RNA_acc_coords) # Measure D-A distance

                if config.MIN_DIST < dist < config.MAX_HB_DIST:
                    dh = vector(donor[0], donor[1])
                    ha = vector(RNA_acc_coords, donor[1])
                    angle = calculate_angle(dh, ha)

                    if config.MIN_HB_ANGLE < angle < config.MAX_HB_ANGLE:
                        modify_HB_result_list(precision, result, dist)

                        if debug:
                            global HB_RNA_acc_info
                            HB_RNA_acc_info += '***\n'
                            HB_RNA_acc_info +=('{} acceptor - {} donor\ndist: {}; angle: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4), round(angle, 4)))
                            HB_RNA_acc_info +=('{}:{}:{}\t{} atom of {}\n'.format(RNA_acc.GetResidue().GetChain(), RNA_acc.GetResidue().GetNum(), RNA_acc.GetResidue().GetAtomID(RNA_acc).strip(), debug_dict_ligand[ligand_name][donor[0]], ligand_name))

                        if detail:
                            detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'HB',
                            debug_dict_ligand[ligand_name][donor[0]],
                            donor[0][0], donor[0][1],  donor[0][2],
                            RNA_acc.GetResidue().GetNum(), RNA_acc.GetResidue().GetChain(), RNA_acc.GetResidue().GetAtomID(RNA_acc).strip(),
                            RNA_acc_coords[0], RNA_acc_coords[1], RNA_acc_coords[2],
                            dist])

                        if precision == 'FULL' and not detail:
                            searching_flag = False
                            break
        else: break

    if searching_flag:

        for RNA_don in donors_RNA:
            if searching_flag:

                RNA_don_coords = np.array([RNA_don[0].GetX(), RNA_don[0].GetY(), RNA_don[0].GetZ()])
                RNA_donH_coords = np.array([RNA_don[1].GetX(), RNA_don[1].GetY(), RNA_don[1].GetZ()])

                for acceptor in ligand_acceptors_coords:

                    dist = measure_distance(RNA_don_coords,acceptor) # Measure D-A distance

                    if config.MIN_DIST < dist < config.MAX_HB_DIST:
                        dh = vector(RNA_don_coords, RNA_donH_coords)
                        ha = vector(acceptor, RNA_donH_coords)
                        angle = calculate_angle(dh, ha)


                        if config.MIN_HB_ANGLE < angle < config.MAX_HB_ANGLE:
                            modify_HB_result_list(precision, result, dist)

                            if debug:
                                global HB_RNA_donor_info
                                HB_RNA_donor_info += '***\n'
                                HB_RNA_donor_info += ('{} donor - {} acceptor\ndist: {}; angle: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4), round(angle, 4)))
                                HB_RNA_donor_info += ('{}:{}:{}\t{} atom of {}\n'.format(RNA_don[0].GetResidue().GetNum(), RNA_don[0].GetResidue().GetChain(), RNA_don[0].GetResidue().GetAtomID(RNA_don[0]), debug_dict_ligand[ligand_name][acceptor], ligand_name))

                            if detail:
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'HB',
                                debug_dict_ligand[ligand_name][acceptor],
                                acceptor[0], acceptor[1],  acceptor[2],
                                RNA_don[0].GetResidue().GetNum(), RNA_don[0].GetResidue().GetChain(), RNA_don[0].GetResidue().GetAtomID(RNA_don[0]),
                                RNA_don_coords[0], RNA_don_coords[1], RNA_don_coords[2],
                                dist])

                            if precision == 'FULL' and not detail:
                                searching_flag = False
                                break

            else: break

    return result


def calculate_HB_no_dha(residue, acceptors_RNA, donors_RNA, ligand_name, ligand_donors_acceptors, precision):
    """Calculates hydrogen bond between residue - ligand pair.
        Simplified graphical representation:\n
        A \***** H --- D\n
            where:
            A  - hydrogen bond acceptor\n
            H  - hydrogen\n
            D  - hydrogen bond acceptor\n
            \* - hydrogen bond\n
        Geometric Rule is:
            1. D-A distance < 3.9 A\n
            If fingerprint's type is FULL:
                - simply checks if there is any hydrogen bonding in residue - ligand pair
            If fingerprint's type is XP:
                - calculates total number of hydrogen bondings in residue - ligand pair
                - assigns type for each hydrogen bonding (strong/moderate/weak), depending on the D-A distance, and calculates number of each type's occurence
                    - strong hydrogen bond:   2.2 A < D-A < 2.5 A
                    - moderate hydrogen bond: 2.5 A < D-A < 3.5 A
                    - weak hydrogen bond:     3.5 A < D-A < 3.9 A

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors
        :param donors_RNA: all tuples (D, H) of the residue
        :param ligand_name: ligand_name^pose_number
        :param ligand_donors_acceptors: [[sublist of acceptors coords], [sublist of tuples (D, H) coords]]
        :param precision: fingerprint type
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type donors_RNA: list
        :type ligand_name: str
        :type ligand_donors_acceptors: list
        :type precision: str
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    if precision == 'XP':
        result = [ligand_name, str(residue.GetNum())+ ':' + str(residue.GetChain()), 0, 0, 0, 0]
    else:
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

                # internal debug mode
                # if RNA_acc.GetResidue().GetNum() == 9: print('res {} {} dist to lig atom {} is {}'.format(RNA_acc.GetResidue().GetNum(), RNA_acc.GetResidue().GetAtomID(RNA_acc), debug_dict[ligand_name][donor[0]], dist))

                if config.MIN_DIST < dist < config.MAX_HB_DIST:
                    modify_HB_result_list(precision, result, dist)

                    if debug:
                        global HB_RNA_acc_info
                        HB_RNA_acc_info += '***\n'
                        HB_RNA_acc_info += ('{} acceptor - {} donor\ndist: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4)))
                        HB_RNA_acc_info += ('{}:{}:{}\t{} atom of {}\n'.format(RNA_acc.GetResidue().GetChain(), RNA_acc.GetResidue().GetNum(), RNA_acc.GetResidue().GetAtomID(RNA_acc).strip(), debug_dict_ligand[ligand_name][donor[0]], ligand_name))

                    if detail:
                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'HB',
                        debug_dict_ligand[ligand_name][donor[0]],
                        donor[0][0], donor[0][1],  donor[0][2],
                        RNA_acc.GetResidue().GetNum(), RNA_acc.GetResidue().GetChain(), RNA_acc.GetResidue().GetAtomID(RNA_acc).strip(),
                        RNA_acc_coords[0], RNA_acc_coords[1], RNA_acc_coords[2],
                        dist])

                    if precision == 'FULL' and not detail:
                        searching_flag = False
                        break
        else: break

    if searching_flag:

        for RNA_don in donors_RNA:

            if searching_flag:

                RNA_don_coords = np.array([RNA_don[0].GetX(), RNA_don[0].GetY(), RNA_don[0].GetZ()])

                for acceptor in ligand_acceptors_coords:
                    dist = measure_distance(RNA_don_coords,acceptor) # Measure D-A distance

                    if config.MIN_DIST < dist < config.MAX_HB_DIST:
                        modify_HB_result_list(precision, result, dist)

                        if debug:
                            global HB_RNA_donor_info
                            HB_RNA_donor_info += '***\n'
                            HB_RNA_donor_info +=('{} donor - {} acceptor\ndist: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4)))
                            HB_RNA_donor_info +=('{}:{}:{}\t{} atom of {}\n'.format(RNA_don[0].GetResidue().GetNum(), RNA_don[0].GetResidue().GetChain(), RNA_don[0].GetResidue().GetAtomID(RNA_don[0]), debug_dict_ligand[ligand_name][acceptor], ligand_name))

                        if detail:
                            detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'HB',
                            debug_dict_ligand[ligand_name][acceptor],
                            acceptor[0], acceptor[1],  acceptor[2],
                            RNA_don[0].GetResidue().GetNum(), RNA_don[0].GetResidue().GetChain(), RNA_don[0].GetResidue().GetAtomID(RNA_don[0]),
                            RNA_don_coords[0], RNA_don_coords[1], RNA_don_coords[2],
                            dist])

                        if precision == 'FULL' and not detail:
                            searching_flag = False
                            break

            else: break

    return result


def calculate_HAL(residue, acceptors_RNA, ligand_name, ligand_donors_coords, precision):
    """Calculates halogen bond between residue - ligand pair.
        Simplified graphical representation:\n
        Y --- O \***** X --- C\n
            where:\n
            Y  - acceptor'; atom covalently bond to the acceptor\n
            O  - halogen bond acceptor\n
            X  - halogen [F, Br, I, Cl]\n
            C  - halogen donor: Carbon\n
            \* - halogen bond\n
        .. note::
            There may be two Ys, if O is part of the ring. If so, we need to apply above rules to both Ys.\n
        Geometric Rules are:
            - X-O distance < 4.0 A
            - C-X-O angle ~ 165 +/- 30
            - X-O-Y angle ~ 120 +/- 30\n
        - If fingerprint's type is FULL:
            - simply checks if there is any halogen bonding in residue - ligand pair
        - If fingerprint's type is XP:
            - calculates total number of halogen bondings in residue - ligand pair

        :param residue: residue as OpenBabel object
        :param acceptors_RNA: residue's hydrogen bond acceptors
        :param ligand_name: ligand_name^pose_number
        :param ligand_donors_coords: [list of tuples (C, halogen)]
        :param precision: fingerprint type
        :type residue: openbabel.OBResidue
        :type acceptors_RNA: list
        :type ligand_name: str
        :type ligand_donors_coords: list
        :type precision: str
        :return: calculated interaction for particular ligand - residue
        :rtype: list
    """

    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0]

    # Important for 'FULL' fingerprint as we are searching for only the first halogen bond
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

                                if precision == 'XP':
                                    result[-1] += 1
                                else:
                                    if not detail:
                                        searching_flag = False # Just found first halogen bond, no need to search further
                                    result[-1] = 1

                                if debug:
                                    global HAL_info
                                    HAL_info += '***\n'
                                    HAL_info += ('{} acceptor - {} donor\ndist: {}; C-X-O angle: {}; X-O-Y angle: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4), round(angle_acc, 4), round(angle_don, 4)))
                                    HAL_info += ('{}:{}:{}\t{} atom of {}\n'.format(RNA_acceptor_set[0].GetResidue().GetChain(), RNA_acceptor_set[0].GetResidue().GetNum(), RNA_acceptor_set[0].GetResidue().GetAtomID(RNA_acceptor_set[0]).strip(), debug_dict_ligand[ligand_name][(donor[0][0], donor[0][1], donor[0][2])], ligand_name))

                                if detail:
                                    global detail_list
                                    detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'HAL',
                                    debug_dict_ligand[ligand_name][(donor[0][0], donor[0][1], donor[0][2])],
                                    donor[0][0], donor[0][1],  donor[0][2],
                                    RNA_acceptor_set[0].GetResidue().GetNum(), RNA_acceptor_set[0].GetResidue().GetChain(), RNA_acceptor_set[0].GetResidue().GetAtomID(RNA_acceptor_set[0]).strip(),
                                    RNA_acc_coords[0], RNA_acc_coords[1], RNA_acc_coords[2],
                                    dist])

                                # If we found halogen bond for one O-Y pair, there is no need to check angles for another Y of the same O (if O has 2 neighbours)
                                acc_bond_found = True
                                break
        else:
           break

    return result


def calculate_CATION_ANION(residue, RNA_anions, ligand_name, ligand_cation_coords, precision):
    """ Calculates cation-anion interaction between residue - ligand pair.
        Simplified graphical representation:\n
        C ***** A\n
            where:\n
            C  - cation\n
            A  - anion\n
        Geometric Rule is:
            - 0.5 A < cation-anion distance < 5.5 A\n
        - If fingerprint's type is FULL:
            - simply checks if there is any cation-anion interaction in residue - ligand pair
        - If fingerprint's type is XP:
            - calculates total number of cation-anion interactions in residue - ligand pair

        :param residue: residue as OpenBabel object
        :param RNA_anions: residue's anions coordinates [OP1, OP2]
        :param ligand_name: ligand_name^pose_number
        :param ligand_cation_coords: list of ligand's cations coords
        :param precision: fingerprint type
        :type residue: openbabel.OBResidue
        :type RNA_anions: list
        :type ligand_name: str
        :type ligand_cation_coords: list
        :type precision: str
        :return: calculated interaction for particular ligand - residue
        :rtype: list

    """

    result = [ligand_name, str(residue.GetNum()) + ':' + str(residue.GetChain()), 0]

    # Important for 'FULL' fingerprint as we are searching for only the first occurance of cation-anion interaction
    searching_flag = True

    for anion in anions_RNA:

        if searching_flag:
            RNA_anion_coords = np.array([anion.GetX(), anion.GetY(), anion.GetZ()])

            for cation in ligand_cation_coords:

                dist = measure_distance(cation, RNA_anion_coords) # Measure cation-anion distance

                if config.MIN_DIST < dist < config.MAX_CA_DIST:

                    if debug:
                        global Cation_Anion_info
                        Cation_Anion_info += '***\n'
                        Cation_Anion_info += ('{} - {}\ndist: {}\n'.format(filename_RNA.split('/')[-1], ligand_name, round(dist, 4)))
                        Cation_Anion_info += ('{}:{}:{}\t{} atom of {}\n'.format(anion.GetResidue().GetChain(), anion.GetResidue().GetNum(), anion.GetResidue().GetAtomID(anion).strip(), debug_dict_ligand[ligand_name][cation], ligand_name))

                    if detail:
                        global detail_list
                        detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'CA',
                        debug_dict_ligand[ligand_name][cation],
                        cation[0], cation[1],  cation[2],
                        anion.GetResidue().GetNum(), anion.GetResidue().GetChain(), anion.GetResidue().GetAtomID(anion).strip(),
                        RNA_anion_coords[0], RNA_anion_coords[1], RNA_anion_coords[2],
                        dist])

                    if precision == 'XP':
                        result[-1]+=1

                    else:
                        result[-1]=1
                        if not detail:
                            searching_flag = False # Just found first cation-anion interaction, no need to search further
                            break
        else:
            break

    return result


def calculate_PI_INTERACTIONS(RNA_rings, RNA_all_atoms, all_ligands_CA_dict, filename_ligand, extension_ligand, precision):
    """Calculates Pi-cation, Pi-anion & Pi-stacking interactions between all residues' aromatic rings - all ligands' pairs.\n
        .. note::
            Important note: each ring of purines is considered separately.\n
        Pi-cation & Pi-anion Geometric Rules:
            - Cation/anion - aromatic ring's center distance < 6.0 A
            - Angle between aromatic ring's planar and cation/anion ~ 90 +/- 30\n
        Three types of Pi-stacking interactions:
            - Sandwich
            - Parallel Displaced
            - T-shaped\n
        Pi-stacking Geometric Rules:
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
        :param filename_ligand: path to ligands input file
        :param extension_ligand: extension of ligands input file
        :param precision: fingerprint type
        :type RNA_rings: list
        :type RNA_all_atoms: list
        :type all_ligands_CA_dict: dict
        :type filename_ligand: str
        :type extension_ligand: str
        :type precision: str
        :return: calculated 3 Pi-interaction for RNA - all ligands
        :rtype: list
       """

    #################################################################
    #  Create dictionary of all ligands' aromatic rings with RDKit  #
    ################################################################

    all_ligands_rings_dict  = findAromaticRingsWithRDKit(filename_ligand)

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

                    if config.MIN_DIST < measure_distance(ion, ring_center_RNA) < config.PI_ION_DISTANCE: # Measure ring center-cation/anion distance

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
                                    Pi_Cation_info +=("{} - {}\ndist: {}; ring's planar to Cation angle: {}\n".format(filename_RNA.split('/')[-1], ligand_name, round(measure_distance(ion, ring_center_RNA), 4), round(pi_ion_angle, 4)))
                                    Pi_Cation_info +=('{}:{}:{}\t{} atom of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_dict_ligand[ligand_name][ion], ligand_name))
                                else:
                                    global Pi_Anion_info
                                    Pi_Anion_info += '***\n'
                                    Pi_Anion_info +=("{} - {}\ndist: {}; ring's planar to Anion angle: {}\n".format(filename_RNA.split('/')[-1], ligand_name, round(measure_distance(ion, ring_center_RNA), 4), round(pi_ion_angle, 4)))
                                    Pi_Anion_info +=('{}:{}:{}\t{} atom {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), debug_dict_ligand[ligand_name][ion], ligand_name))

                            if detail:
                                if j == 0:
                                    detail_inter = 'Pi_Cation'
                                else:
                                    detail_inter = 'Pi_Anion'
                                detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], detail_inter,
                                debug_dict_ligand[ligand_name][ion],
                                ion[0], ion[1],  ion[2],
                                residue.GetNum(), residue.GetChain(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]),
                                ring_center_RNA[0], ring_center_RNA[1], ring_center_RNA[2],
                                measure_distance(ion, ring_center_RNA)])

        ########################################
        #  Calculate Pi-Stacking interactions  #
        ########################################

        for ligand_name in all_ligands_rings_dict.keys():

            results[2].append([ligand_name,  str(residue.GetNum()) + ':' + str(residue.GetChain()), 0])

            ligand_rings = all_ligands_rings_dict[ligand_name] # Take list of ligand's aromatic rings atoms coords

            if debug:
                global arom_ring_ligands_info
                arom_ring_ligands_info[ligand_name] = str(len(all_ligands_rings_dict[ligand_name]))

            for aring in ligand_rings:

                ring_atoms_ligand = aring
                atoms_creating_planar_space_ligand = np.array([ring_atoms_ligand[0],ring_atoms_ligand[1],ring_atoms_ligand[2]], dtype=np.longdouble) # Add 3 atoms (we do not need more) from ring to calculate planar
                planar_ligand = calculate_planar(atoms_creating_planar_space_ligand)
                ring_center_ligand = centroid([ra for ra in ring_atoms_ligand])
                centroid_distance = measure_distance(ring_center_RNA,ring_center_ligand) # Measure ring - ring distance
                planar_angle = calculate_angle(planar_RNA,planar_ligand)

                # Calculate aromatic rings' center offset (project each ring center into the other ring)
                proj1 = projection(planar_ligand, ring_center_ligand, ring_center_RNA)
                proj2 = projection(planar_RNA, ring_center_RNA, ring_center_ligand)
                offset = min(measure_distance(proj1, ring_center_ligand), measure_distance(proj2, ring_center_RNA))

                if  config.MIN_DIST < centroid_distance < config.RING_RING_MAX and offset < config.PISTACK_OFFSET_MAX:
                    planar_angle = calculate_angle(planar_RNA, planar_ligand) # Calculate planar - planar angle, which is equal to angle between normal vectors of both planes

                    if planar_angle > 90: planar_angle = 180 - planar_angle # We are 'on the other side'

                    if planar_angle < config.PI_ANGLE_DISPLACED: # Sandwich & Displaced Pi-stacking interactions
                        results[2][-1][-1] = 1

                        if debug:
                            global Sandwich_Displaced_info
                            Sandwich_Displaced_info += '***\n'
                            Sandwich_Displaced_info += "{} - {}\nrings center dist: {}; rings offset: {}; rings planars angle: {}\n".format(filename_RNA.split('/')[-1], ligand_name, round(centroid_distance, 4), round(offset, 4), round(planar_angle, 4))
                            Sandwich_Displaced_info += '{}:{}:{}\t{} atoms of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), ring_atoms_ligand, ligand_name)

                        if detail:
                            debug_ligand = ''
                            for i in range(len(ring_atoms_ligand)):
                                debug_ligand += str(debug_dict_ligand[ligand_name][ring_atoms_ligand[i]])
                                if i != len(ring_atoms_ligand)-1:
                                    debug_ligand += ','

                            detail_list.append([ligand_name.split('^')[0], ligand_name.split('^')[1], 'Pi_Stacking',
                            debug_ligand,
                            ring_center_ligand[0], ring_center_ligand[1],  ring_center_ligand[2],
                            residue.GetNum(), residue.GetChain(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]),
                            ring_center_RNA[0], ring_center_RNA[1], ring_center_RNA[2],
                            centroid_distance])

                    elif abs(config.PLANAR_ANGLE_TSHAPED - planar_angle) < config.PLANAR_ANGLE_TSHAPED_DEV: # T-shaped Pi-stacking interaction
                        results[2][-1][-1] = 1

                        if debug:
                            global T_shaped_info
                            T_shaped_info += '***\n'
                            T_shaped_info += "{} - {}\nrings center dist: {}; rings offset: {}; rings planars angle: {}\n".format(filename_RNA.split('/')[-1], ligand_name, round(centroid_distance, 4), round(offset, 4), round(planar_angle, 4))
                            T_shaped_info += '{}:{}:{}\t{} atoms of {}\n'.format(residue.GetChain(), residue.GetNum(), ','.join([debug_dict_rna[ring_atoms_RNA[i].coords] for i in range(len(ring_atoms_RNA))]), ring_atoms_ligand, ligand_name)

                    else:
                        continue

        ###########################################################
        #  Merge all the Pi-interactions results and return them  #
        ###########################################################

        for j in range(3): # Append results from calculated 3 Pi-interactions: Pi-cation, Pi-anion, Pi-stacking
            RESULTS[j].extend(results[j])

    return RESULTS


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

    parser = argparse.ArgumentParser(description = '''Script calculating Structural Interaction Fingerprints (SIFs) in RNA/DNA - ligand complexes.''',
                                     epilog = 'If no optional -o parameter was passed, script will create outputs/ directory in the current working directory and save there SIFs in tsv format.',
                                     add_help = False,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter,
                                     conflict_handler = 'resolve')

    required_arguments = parser.add_argument_group('Required arguments')
    required_arguments.add_argument('-r', help='pass RNA/DNA structure in pdb/mol2 format', required=True, metavar='RNA/DNA', default=argparse.SUPPRESS)
    required_arguments.add_argument('-l', help='pass ligands file in sdf format', required=True, metavar='LIGANDS', default=argparse.SUPPRESS)

    optional_arguments = parser.add_argument_group('Optional arguments')
    optional_arguments.add_argument('-f', help='pass fingerprint type, available types are SIMPLE, PBS, FULL, XP', default='FULL', metavar='TYPE')
    optional_arguments.add_argument('-o', help='pass output name', metavar='NAME')
    optional_arguments.add_argument('-dha', help='consider Donor-Hydrogen-Acceptor angle in hydrogen bonds calculation', action='store_true')
    optional_arguments.add_argument('-wrapper', help='pass results wrapper types (multiple types possible at once, but have to be comma separated)', metavar='WRAPPER')
    optional_arguments.add_argument('-print', help='print found interactions on screen', action='store_true')
    optional_arguments.add_argument('-detail', help='generate an additional file with detailed data on detected interactions (used for detail visualization)', action='store_true')
    optional_arguments.add_argument('-vis', help='make heatmap visualization', action='store_true')
    optional_arguments.add_argument('-verbose', help='provides additional details about calculations performed at the given moment', action='store_true')
    optional_arguments.add_argument('-debug', help='enter debug mode', action='store_true')
    optional_arguments.add_argument('-h', action = 'help', help = 'show this help message and exit')
    optional_arguments.add_argument('--help', '-h', action = 'help', help = 'show this help message and exit')

    args = vars(parser.parse_args())

    filename_RNA = args['r']
    filename_ligand=args['l']
    extension_ligand = ".".join(filename_ligand.split('/')[-1].split('.')[1:])
    extension_structure = ".".join(filename_RNA.split('/')[-1].split('.')[1:])
    fingerprint = args['f']
    output = args['o']
    consider_dha = args['dha']

    try:
        wrapper = args['wrapper'].split(',')
    except AttributeError:
        wrapper = None

    print_flag = args['print']
    detail = args['detail']
    visualization = args['vis']
    verbose = args['verbose']
    debug = args['debug']

    if extension_structure not in ['pdb', 'mol2']:
        raise Exception('RNA/DNA structure has to be in pdb or mol2 format!')

    if extension_ligand != 'sdf':
        raise Exception("Ligands' structures have to be in sdf format!")

    #########################
    #  FINGERPRINT CALLING  #
    #########################

    FUNCTIONS = {'SIMPLE': 1, 'PBS': 3, 'FULL': 6, 'XP': 9}
    WRAPPERS = {'ACUG': 4, 'PuPy':2, 'Counter': 1 }
    ANALYSIS_NAME = []

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
    # Create empty list of all nucleic acid residues
    RNA_residues = []
    # Create empty list of all nucleic acid residues
    RNA_nucleotides = []
    RNA_LENGTH = structure.OBMol.NumResidues()

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

        for residue in openbabel.OBResidueIter(structure.OBMol): # Loop over all nucleic acid residues
            RNA_residues.append(str(residue.GetNum())+ ':' + str(residue.GetChain()))
            RNA_nucleotides.append(str(residue.GetName()))

            for ligand_name, ligand_atoms in ligands_all_atoms.items():

                centroid_ligand = centroid(ligand_atoms)

                if fingerprint == 'SIMPLE':
                    result = calculate_SIMPLE(residue, ligand_name, ligand_atoms, centroid_ligand)
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0)
                else:
                    result = calculate_PBS(residue, ligand_name, ligand_atoms, centroid_ligand)
                    if sum(result[-3:]) != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0, True) # Assign each of 3 (P/B/S) residue-ligand interaction

    else: # fingerprint == 'FULL' or fingerprint == 'XP'

        # Read all ligands
        ligands_mols = list(pybel.readfile(extension_ligand, filename_ligand))

        # Add missing hydrogens to all ligands
        for i in range(len(ligands_mols)):
            ligands_mols[i].OBMol.AddHydrogens()

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
                HAL_info = ''
                Cation_Anion_info = ''
                Pi_Cation_info = ''
                Pi_Anion_info = ''
                Sandwich_Displaced_info = ''
                T_shaped_info = ''

        ####################################################################

        # Create ligands' all atoms dictionary
        ligands_all_atoms = find_ligands_all_atoms(ligands_mols)
        # Create ligands acceptors' & donors' dictionary
        ligands_hba_hbd = find_ligands_HBA_HBD(ligands_mols, verbose)
        # Create ligands halogens' donors dictionary
        ligands_HAL = find_ligands_HAL_don(ligands_mols, verbose)
        # Create ligands cations' & anions' dictionary
        ligands_CA = find_ligands_CA(ligands_mols, verbose)
        # Find all RNA rings
        rings_RNA = find_RNA_rings(structure, extension_structure)

        # Fill the RESULTS dictionary of keys - ligand ids and values - lists of 0
        for ligand_name in ligands_hba_hbd.keys():
            RESULTS[ligand_name] = [0] * RNA_LENGTH * FUNCTIONS[fingerprint]

        for residue in openbabel.OBResidueIter(structure.OBMol): # Loop over all nucleic acid residue to calculate hydrogen bondings, halogen bondings & cation-anion interactions
            RNA_residues.append(str(residue.GetNum())+ ':' + str(residue.GetChain()))
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
            for atom in openbabel.OBResidueAtomIter(residue):
                if atom.GetAtomicNum() != config.HYDROGEN_NUM:
                    residue_atoms_coords.append(np.array([atom.GetX(), atom.GetY(), atom.GetZ()]))

            for ligand_name_HB, ligand_values_HB in ligands_hba_hbd.items():
                if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_HB])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                    continue

                if consider_dha:
                    result = calculate_HB(residue, acceptors_RNA, donors_RNA, ligand_name_HB, ligand_values_HB, fingerprint)
                else:
                    result = calculate_HB_no_dha(residue, acceptors_RNA, donors_RNA, ligand_name_HB, ligand_values_HB, fingerprint)

                if fingerprint == 'FULL':
                    if result[-1] != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0)
                else:
                    if sum(result[-4:]) != 0:
                        assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 0, True) # XP fingerprint holds information about 4 hydrogen bonds types (total number/strong number/moderate number/weak number)

            for ligand_name_HAL, ligand_values_HAL in ligands_HAL.items():
                if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_HAL])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                    continue

                result = calculate_HAL(residue, acceptors_RNA, ligand_name_HAL, ligand_values_HAL, fingerprint)

                if result[-1] != 0:
                    assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 1)

            for ligand_name_CA, ligand_values_CA in ligands_CA.items():
                if measure_distance(centroid(residue_atoms_coords), centroid(ligands_all_atoms[ligand_name_CA])) > config.RES_LIGAND_MAX_DIST: # nucleic acid residue centroid and ligand centroid are futher than declared threshold, no chance for any contact
                    continue

                result = calculate_CATION_ANION(residue, anions_RNA, ligand_name_CA, ligand_values_CA[0], fingerprint)
                if result[-1] != 0:
                    assign_interactions_results(result, RESULTS, RNA_LENGTH, len(RNA_residues)-1, FUNCTIONS[fingerprint], 2)


        PI_INTERACTIONS = calculate_PI_INTERACTIONS(rings_RNA, structure.atoms, ligands_CA, filename_ligand, extension_ligand, fingerprint) # Calculate Pi-cation, Pi-anion & Pi-stacking interactions

        for i in range(3):

            result = PI_INTERACTIONS[i] # i=0 Pi-cation interaction; i=1 Pi-anion interaction; i=2 Pi-stacking interaction
            for res in result:
                if res[-1] != 0: # We assign only 1
                    assign_interactions_results(res, RESULTS, RNA_LENGTH, RNA_residues.index(res[1]), FUNCTIONS[fingerprint], i+3)

        if debug:
            print_debug_info(ligands_hba_hbd, ligands_HAL, ligands_CA, arom_ring_ligands_info, debug_dict_ligand,
            RNA_HB_acc_don_info, RNA_anion_info, arom_RNA_ligands_info, HB_RNA_acc_info, HB_RNA_donor_info,
            HAL_info, Cation_Anion_info, Pi_Cation_info, Pi_Anion_info, Sandwich_Displaced_info, T_shaped_info, columns)

    # Wrap results if wrapper was passed
    if wrapper:

        WRAP_RESULTS = {}
        for w in wrapper:
            if verbose: print('Wrapping fingerprint type {} results to {} wrapper'.format(fingerprint, w))
            WRAP_RESULTS[w] = wrap_results(w, RESULTS, RNA_nucleotides, FUNCTIONS[fingerprint], WRAPPERS[w])

    # Create dataframe
    columns = {'SIMPLE': ['SIMPLE'],
               'PBS': ['P', 'B', 'S'],
               'FULL': ['HB', 'HAL', 'CA', 'Pi_Cation', 'Pi_Anion', 'Pi_Stacking'],
               'XP': ['How_many_HB','HB_Strong','HB_Moderate','HB_Weak','How_many_HAL','How_many_CA','Pi_Cation','Pi_Anion', 'Pi_Stacking']}

    is_structure_RNA = check_if_RNA(RNA_nucleotides)

    if is_structure_RNA:
        nucleotides_letters = ['A','C','U','G']
    else:
        nucleotides_letters = ['A','C','T','G']

    if verbose: print('Calculations completed, saving the results...')

    for analysis in ANALYSIS_NAME:

        if analysis == 'ACUG':
            DF_COLUMNS = [filename_RNA.split('/')[-1] + '#' + res + '#' + fing_type for res in nucleotides_letters for fing_type in columns[fingerprint]]
        elif analysis == 'PuPy':
            DF_COLUMNS = [filename_RNA.split('/')[-1] + '#' + res + '#' + fing_type for res in ['Purines','Pyrimidynes'] for fing_type in columns[fingerprint]]
        elif analysis == 'Counter':
             DF_COLUMNS = [filename_RNA.split('/')[-1] + '#' + fing_type for fing_type in columns[fingerprint]]
        else:
            DF_COLUMNS = [filename_RNA.split('/')[-1] + '#' + res + '#' + fing_type for res in RNA_residues for fing_type in columns[fingerprint]]

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
            detail_columns = ['Ligand_name', 'Ligand_pose', 'Interaction', 'Ligand_Atom',
                             'Ligand_X', 'Ligand_Y', 'Ligand_Z',
                             'Receptor_Number', 'Receptor_Chain', 'Receptor_Atom',
                             'Receptor_X', 'Receptor_Y', 'Receptor_Z', 'Distance']
            detail_df = pd.DataFrame(detail_list, columns=detail_columns)

        # Save output as tsv
        if output:
            output_proper = output
            if analysis in FUNCTIONS.keys():
                if output[-1] == '/' or output[-1] == '\\': # default output name, location specified
                    output_proper += filename_RNA.split('/')[-1] + '_' + filename_ligand.split('/')[-1] + '_' + fingerprint
            else:
                if output[-1] == '/' or output[-1] == '\\': # default output name, location specified
                    output_proper += filename_RNA.split('/')[-1] + '_' + filename_ligand.split('/')[-1] + '_' + fingerprint + '_' + analysis

            ALL_FINGERPRINTS_DF.to_csv('%s.tsv' %output_proper, sep='\t')
            if detail:
                  detail_save = output + 'detail_' + filename_RNA.split('/')[-1] + '_' + filename_ligand.split('/')[-1] + '_' + fingerprint
                  detail_df.to_csv('%s.tsv' %detail_save, sep='\t' )

        else:
            if not os.path.exists('outputs'): os.makedirs('outputs')
            if analysis in FUNCTIONS.keys():
                ALL_FINGERPRINTS_DF.to_csv('outputs/%s_%s_%s.tsv' %(filename_RNA.split('/')[-1],filename_ligand.split('/')[-1], fingerprint), sep='\t')
            else:
                ALL_FINGERPRINTS_DF.to_csv('outputs/%s_%s_%s_%s.tsv' %(filename_RNA.split('/')[-1],filename_ligand.split('/')[-1], fingerprint, analysis), sep='\t')

            if detail:
                  detail_save = 'outputs/detail_%s_%s_%s.tsv' %(filename_RNA.split('/')[-1],filename_ligand.split('/')[-1], fingerprint)
                  detail_df.to_csv('%s.tsv' %detail_save, sep='\t' )

    # Print found interactions on screen
        if print_flag:
            interact_names = {'HB': 'Hydrogen Bonds', 'HAL': 'Halogen Bonds', 'CA': 'Cation-Anion', 'Pi_Cation': 'Pi-Cation',
                              'Pi-Anion': 'Pi-Anion', 'Pi_Stacking':'Pi-Stacking', 'How_many_HB': 'Total Hydrogen Bonds',
                              'HB_Strong': 'No of strong Hydrogen Bonds', 'HB_Moderate': 'No of moderate Hydrogen Bonds',
                              'HB_Weak': 'No of weak Hydrogen Bonds', 'How_many_HAL': 'Total Halogen Bonds', 'How_many_CA': 'Total Cation-Anion',
                              'P':'Phosphate contact', 'B': 'Base contact', 'S':'Sugar contact', 'SIMPLE':'contact'}
            for index, row in ALL_FINGERPRINTS_DF.iterrows():
                print('# {} - {} #'.format(filename_RNA.split('/')[-1], index))
                print()
                for el in range(len(row)):
                    if row[el] > 0:
                        s = DF_COLUMNS[el].split('#')[1:]
                        print('{}\t{}\t{}'.format(s[0],interact_names[s[1]], row[el]))
                print()
            print_flag = False

    # Visualize output as heatmap
        if visualization:

            fig = plt.figure()
            height = len(DF_INDEXES)*3

            width_multipliers = {'SIMPLE' : 1, 'PBS' : 1.5, 'FULL' : 3, 'XP' : 3.5}

            if analysis in FUNCTIONS.keys():
                width = len(RNA_residues) * width_multipliers[fingerprint]
            else:
                if analysis == 'ACUG':
                    width = 18
                elif analysis == 'PuPy':
                    width = 14
                else:
                    width = 12

            plt.figure(figsize=(width, height))
            plt.axes(aspect='equal')

            if fingerprint == 'XP' or analysis == 'Counter' :
                cmap = colors.ListedColormap(['cornflowerblue','turquoise','bisque', 'orange', 'lightcoral', 'darkred'])
            else:
                cmap = colors.ListedColormap(['cornflowerblue','turquoise'])

            max_value = ALL_FINGERPRINTS_DF.to_numpy().max()
            bounds = np.arange(-0.5, max_value + 1, 1)
            norm = colors.BoundaryNorm(bounds,cmap.N)

            heatmap = plt.pcolormesh(ALL_FINGERPRINTS_DF, cmap=cmap, norm=norm, edgecolors='silver')

            plt.yticks(np.arange(0.5, len(ALL_FINGERPRINTS_DF.index), 1), ALL_FINGERPRINTS_DF.index, fontsize = 14)

            if analysis in FUNCTIONS.keys(): # no wrapper
                x = list(res + '#' + fing_type for res in RNA_residues for fing_type in columns[fingerprint])
                plt.xticks(np.arange(0.5, len(x), 1), x, fontsize = 14, rotation = 90)

            else:
                if analysis == 'PuPy':
                    x = list(res + '#' + fing_type for res in ['Purines','Pyrimidynes'] for fing_type in columns[fingerprint])
                    plt.xticks(np.arange(0.5, len(x), 1), x, fontsize = 14, rotation = 90)
                elif analysis == 'Counter':
                    x = list(res + '#' + fing_type for res in ['Total'] for fing_type in columns[fingerprint])
                    plt.xticks(np.arange(0.5, len(x), 1), x, fontsize = 14, rotation = 90)
                else:
                    x = list(res + '#' + fing_type for res in nucleotides_letters for fing_type in columns[fingerprint])
                    plt.xticks(np.arange(0.5, len(x), 1), x, fontsize = 14, rotation = 90)

            plt.colorbar(heatmap,ticks = range(max_value+1), shrink = 0.2)

            if output:
                if analysis in FUNCTIONS.keys():
                    plt.tight_layout()
                    plt.savefig('%s_%s.png' %(output, fingerprint), dpi = 300)
                else:
                    plt.tight_layout()
                    plt.savefig('%s_%s_%s.png' %(output, fingerprint, analysis), dpi = 300)
            else:
                if analysis in FUNCTIONS.keys():
                    plt.tight_layout()
                    plt.savefig('outputs/%s_%s_%s.png' %(filename_RNA.split('/')[-1],filename_ligand.split('/')[-1],fingerprint), dpi = 300)
                else:
                    plt.tight_layout()
                    plt.savefig('outputs/%s_%s_%s_%s.png' %(filename_RNA.split('/')[-1],filename_ligand.split('/')[-1], fingerprint, analysis), dpi = 300)

        print('{} results saved successfully!'.format(analysis))
