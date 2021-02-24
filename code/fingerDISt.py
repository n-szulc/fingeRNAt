#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
fingerDIST is a software to calculate different Distance Metrics based on Structural Interaction Fingerprint (SIFt).

Authors:
Natalia A. Szulc, nszulc@iimcb.gov.pl
Filip Stefaniak, fstefaniak@genesilico.pl

If you use this software, please cite:
Natalia A. Szulc, Zuzanna Mackiewicz, Janusz M. Bujnicki, Filip Stefaniak
[in preparation]

Requires Python 3.5 - 3.8
'''

import argparse
import os
import pandas as pd
import numpy as np
import shutil
import DistanceMetrics as DM
from tqdm import tqdm

np.set_printoptions(suppress=True)

metric_name_to_function = {'manhattan' : (DM.Similarity, 'manhattan_distance'), 'square_euclidean' : (DM.Similarity, 'square_euclidean_distance'), \
                           'euclidean' : (DM.Similarity, 'euclidean_distance'), 'half_square_euclidean' : (DM.Similarity, 'half_square_euclidean_distance'), \
                           'cosine_similarity' : (DM.Similarity, 'cosine_similarity'), 'tanimoto' : (DM.Similarity, 'tanimoto_coefficient'), \
                           'soergel_distance' : (DM.Similarity, 'soergel_distance')}

if __name__ == "__main__":

    #######################
    #  ARGUMENTS PARSING  #
    #######################


    welcome_mssg = '# Welcome to fingerDISt! #'
    columns = shutil.get_terminal_size().columns
    print('')
    print(('#'*len(welcome_mssg)).center(columns))
    print(welcome_mssg.center(columns))
    print(('#'*len(welcome_mssg)).center(columns))

    parser = argparse.ArgumentParser(description = '''Script calculating different Distance Metrics between Structural Interaction Fingerprint (SIFt).''',
                                     epilog = 'If no optional -o parameter was passed, script will create outputs/ directory in the current working directory and save there calculated Distance Metrics in tsv format.',
                                     add_help = False,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter,
                                     conflict_handler = 'resolve')

    required_arguments = parser.add_argument_group('Required arguments')
    required_arguments.add_argument('-i', help='pass SIFt in tsv/csv format', required=True, metavar='SIFt', default=argparse.SUPPRESS)
    required_arguments.add_argument('-m', help='pass types of desired Distance Metrics, available types are: tanimoto, cosine_similarity, manhattan, euclidean, square_euclidean, half_square_euclidean, soergel_distance', required=True, metavar='METRICS', default=argparse.SUPPRESS)

    optional_arguments = parser.add_argument_group('Optional arguments')
    optional_arguments.add_argument('-o', help='pass output path or name', metavar='NAME')
    optional_arguments.add_argument('-verbose', help='prints calculated Distance Metrics on the screen', action='store_true')
    optional_arguments.add_argument('-h', action = 'help', help = 'show this help message and exit')
    optional_arguments.add_argument('--help', '-h', action = 'help', help = 'show this help message and exit')

    args = vars(parser.parse_args())
    filename_SIFt = args['i']
    metrics=args['m'].split(',')
    extension_SIFt= filename_SIFt.split('/')[-1].split('.')[-1]
    output = args['o']
    verbose = args['verbose']

    if extension_SIFt == 'csv':
        sep = ','
    elif extension_SIFt == 'tsv':
        sep = '\t'
    else:
        raise Exception('Unknown SIFt extension')

    for m in metrics:
        if m not in metric_name_to_function.keys():
            raise Exception('Unknown metrics')

    # Define helper function

    def check_binary_values(matrix):
        """Checks if given matrix has only binary values

        :param matrix: matrix of SIFt
        :type matrix: numpy.ndarray
        :return: True if matrix has only binary values, False otherwise
        :rtype: bool
        """

        result = np.array_equal(matrix, matrix.astype(bool))
        return result

    # Read file with SIFt and prepare it
    df = pd.read_csv(filename_SIFt, sep=sep)
    df.fillna(value=0, inplace=True)
    ligands = list(df['Ligand_name'])
    df = df.drop(['Ligand_name'], axis=1)
    fingerprint_matrix = df.to_numpy()

    # Distance Metrics calculations

    for m in metrics:
        if m == 'tanimoto':
            check_tanimoto = check_binary_values(fingerprint_matrix)
            if not check_tanimoto:
                raise Exception('Incorrect values for Tanimoto coefficient!')

        results = np.zeros((len(fingerprint_matrix),len(fingerprint_matrix)))
        for i in tqdm(range(len(fingerprint_matrix))):
            for j in range(len(fingerprint_matrix)):
                cls, method = metric_name_to_function[m]
                results[i][j] = getattr(cls(0.0001), method)(fingerprint_matrix[i], fingerprint_matrix[j])

        if verbose:
            print(('#'*(len(m)+4)))
            print('# {} #'.format(m))
            print(('#'*(len(m)+4)))
            print()
            print('\n'.join('\t'.join(str(cell) for cell in row) for row in np.round(results, 4)))

        # Save results to DataFrame

        df_results = pd.DataFrame(data=np.round(results, 4), index=ligands, columns=ligands)
        df_results.index.name = 'Ligand_name'

        if output:
            output_proper = output
            save_name = filename_SIFt.split('/')[-1] + '_' + m
            if output[-1] == '/' or output[-1] == '\\': # default output name, location specified
                output_proper += save_name
            else:
                output_proper += '/' + save_name
            df_results.to_csv('%s.%s' %(output_proper, extension_SIFt), sep=sep)
        else:
            if not os.path.exists('outputs'): os.makedirs('outputs')
            output_proper = filename_SIFt.split('/')[-1] + '_' + m
            df_results.to_csv('outputs/%s.%s' %(output_proper, extension_SIFt), sep=sep)

        print()
        print('{} scores saved successfully!'.format(m))
        print('')
