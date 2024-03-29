'''
The following code is based on
github.com/varun-parthasarathy/crux-fr-sprint/blob/master/DistanceMetrics.py
(MIT license)
with few additional metrics.
'''


import math
import numpy as np
import hashlib


class Similarity:
    """
    This class contains instances of similarity / distance metrics. These are used in centroid based clustering
    algorithms to identify similar patterns and put them into the same homogeneous sub sets
    :param minimum: the minimum distance between two patterns (so you don't divide by 0)
    """

    def __init__(self, minimum):
        self.e = minimum
        self.vector_operators = VectorOperations()

    def manhattan_distance(self, p_vec, q_vec):
        """
        This method implements the manhattan distance metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the manhattan distance between vector one and two
        """
        return max(np.sum(np.fabs(p_vec - q_vec)), self.e)

    def square_euclidean_distance(self, p_vec, q_vec):
        """
        This method implements the squared euclidean distance metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the squared euclidean distance between vector one and two
        """
        diff = p_vec - q_vec
        return max(np.sum(diff**2), self.e)

    def euclidean_distance(self, p_vec, q_vec):
        """
        This method implements the euclidean distance metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the euclidean distance between vector one and two
        """
        return max(math.sqrt(self.square_euclidean_distance(p_vec, q_vec)), self.e)

    def half_square_euclidean_distance(self, p_vec, q_vec):
        """
        This method implements the half squared euclidean distance metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the half squared euclidean distance between vector one and two
        """
        return max(0.5 * self.square_euclidean_distance(p_vec, q_vec), self.e)

    def cosine_similarity(self, p_vec, q_vec):
        """
        This method implements the cosine similarity metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the cosine similarity between vector one and two
        """
        pq = self.vector_operators.product(p_vec, q_vec)
        p_norm = self.vector_operators.norm(p_vec)
        q_norm = self.vector_operators.norm(q_vec)
        return max(pq / (p_norm * q_norm), self.e)

    def tanimoto_coefficient(self, p_vec, q_vec):
        """
        This method implements the cosine tanimoto coefficient metric
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the tanimoto coefficient between vector one and two
        """
        pq = self.vector_operators.product(p_vec, q_vec)
        p_sum = self.vector_operators.vec_sum(p_vec)
        q_sum = self.vector_operators.vec_sum(q_vec)
        return max(pq / (p_sum + q_sum - pq), self.e)

    def tversky(self, p_vec, q_vec):
        """
        Implemented based on https://docs.eyesopen.com/toolkits/python/graphsimtk/measure.html
        Alpha, beta values taken from Leung et al.
        "SuCOS is Better than RMSD for Evaluating Fragment Elaboration and Docking Poses'"
        """

        pq = self.vector_operators.product(p_vec, q_vec)
        # p_sum = self.vector_operators.vec_sum(p_vec)
        q_sum = self.vector_operators.vec_sum(q_vec)

        # use this formula if you want to manupulate alpha and beta
        # alpha=1
        # beta=0
        # return max(pq / (alpha*(p_sum-pq) + beta*(q_sum - pq) + pq), self.e)

        # default:
        # faster formula for:         alpha=1                beta=0
        # return max(pq / p_sum, self.e) # reference in rows
        return max(pq / q_sum, self.e) # reference is in columns, compared molecule in rows

    def soergel(self, p_vec, q_vec):
        """
        This method implements the Soergel distance metric
        See: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0069-3
        and: https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics_OLCC_(2019)/6%3A_Molecular_Similarity/6.2%3A_Similarity_Coefficients
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the Soergel distance between vector one and two
        """
        pq = self.vector_operators.product(p_vec, q_vec)
        p_sum = self.vector_operators.vec_sum(p_vec)
        q_sum = self.vector_operators.vec_sum(q_vec)
        return max( (p_sum + q_sum - 2*pq) / (p_sum + q_sum - pq), self.e)


    @staticmethod
    def get_key(p_vec, q_vec):
        """
        This method returns a unique hash value for two vectors. The hash value is equal to the concatenated string of
        the hash value for vector one and vector two. E.g. is hash(p_vec) = 1234 and hash(q_vec) = 5678 then get_key(
        p_vec, q_vec) = 12345678. Memoization improved the speed of this algorithm 400%.
        :param p_vec: vector one
        :param q_vec: vector two
        :return: a unique hash
        """
        # return str(hash(tuple(p_vec))) + str(hash(tuple(q_vec)))
        return str(hashlib.sha1(p_vec)) + str(hashlib.sha1(q_vec))


class VectorOperations():
    """
    This class contains useful implementations of methods which can be performed on vectors
    """

    @staticmethod
    def product(p_vec, q_vec):
        """
        This method returns the product of two lists / vectors
        :param p_vec: vector one
        :param q_vec: vector two
        :return: the product of p_vec and q_vec
        """
        return np.dot(p_vec, q_vec)

    @staticmethod
    def vec_sum(p_vec):
        """
        This method returns the square of a vector
        :param p_vec: the vector to be squared
        :return: the squared value of the vector
        """
        return sum(p_vec)

    @staticmethod
    def norm(p_vec):
        """
        This method returns the norm value of a vector
        :param p_vec: the vector to be normed
        :return: the norm value of the vector
        """
        return np.linalg.norm(p_vec)
