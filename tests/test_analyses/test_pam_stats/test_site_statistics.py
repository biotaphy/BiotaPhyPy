"""Module for testing site statistic computations.

Notes:
    * Uses pytest style testing.
"""
import os

import dendropy
import numpy as np
import pytest

from lmpy import TreeWrapper

import biotaphy.analyses.helpers.data_readers as data_readers
from random import shuffle


# .............................................................................
class Test_calculate_tree_site_statistics:
    """Test calculate tree site statistics."""
    # ................................
    def test_calculate_tree_site_statistics_simple(self):
        """Perform a simple test of calculate_tree_site_statistics.

        This is a simple test of the calculate_tree_site_statistics function.
        It creates a simple tree and PAM and tests that the metrics generated
        are the same as the expected values.
        """
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        squid_dict = {
            'A': 'squidA',
            'B': 'squidB',
            'C': 'squidC',
            'D': 'squidD',
            'E': 'squidE'
        }
        tree.annotate_tree_tips('squid', squid_dict)
        pam = Matrix(
            np.array([
                [1, 0, 1, 0, 1],
                [1, 1, 0, 1, 0],
                [0, 1, 0, 0, 1],
                [1, 1, 1, 1, 0],
                [0, 0, 1, 1, 1],
                ]),
            headers={
                '0': ['site_0', 'site_1', 'site_2', 'site_3', 'site_4'],
                '1': ['squidA', 'squidB', 'squidC', 'squidD', 'squidE']})
        test_stats = Matrix(
            np.array([
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                ]))
        assert np.all(
            np.isclose(test_stats, calculate_tree_site_statistics(pam, tree)))

# .............................................................................
class Test_calculate_tree_statistics:
    """Test calculate_tree_statistics."""
    # ................................
    def test_calculate_tree_statistics_simple(self):
        """Perform a simple test of calculate_tree_statistics.

        This is a simple test of the calculate_tree_statistics function.  It
        creates a simple tree with known metric values and ensures the results
        are what are expected.
        """
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert calculate_tree_statistics(tree) == (
            35.0, 35.0, 20.75, 27.5, 42.5, 49.25,
            23.75, 20.0, 10.0, 10.0, 32.5, 48.25)


# .............................................................................
class Test_get_subtree:
    """Tests the get_subtree function"""
    # ................................
    def test_get_subtree_simple(self):
        """Perform a simple test of get_subtree.

        This is a simple test that will create a tree with squids at some of
        the tips and some tips without.  The test will then subset tree and
        check that the resulting tree only includes the tips with the specified
        squids.
        """
        keep_percentage = 0.5
        # Tip pool
        tip_pool = [
            ('A', 'squid_A'),
            ('B', 'squid_B'),
            ('C', None),
            ('D', 'squid_D'),
            ('E', 'squid_E'),
            ('F', None),
            ('G', 'squid_G'),
            ('H', 'squid_H'),
            ('I', None),
            ('J', 'squid_J')
            ]
        # Build a tree
        tree_str = '({});'.format(','.join([tip[0] for tip in tip_pool]))
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        tree.resolve_polytomies()
        # Add some squids
        squid_dict = {}
        for tip_name, squid in tip_pool:
            if squid is not None:
                squid_dict[tip_name] = squid
        tree.annotate_tree_tips('squid', squid_dict)

        # Determine which tips to keep
        keep_tips = []
        while len(keep_tips) < 3:
            keep_tips = []
            for tip_idx in np.where(
                    np.random.random((len(tip_pool),)) < keep_percentage)[0]:
                if tip_pool[tip_idx][1] is not None:
                    keep_tips.append(tip_pool[tip_idx])
        # Get subtree
        subtree = get_subtree(tree, [tip[1] for tip in keep_tips], 'squid')
        # Go through subtree tips and make sure they match keep tips
        kept_labels = [taxon.label for taxon in tree.taxon_namespace]
        assert len(kept_labels) == len(keep_tips)
        for tip in keep_tips:
            assert tip[0] in kept_labels


# .............................................................................
def calculate_tree_site_statistics(pam, tree):
    """Calculate tree statistics for each site in a PAM.

    Args:
        pam (Matrix): A PAM to compute statistics for.
        tree (TreeWrapper): A corresponding tree to use for statistics.

    Returns:
        Matrix - A site by statistic matrix of tree pam statistics.
    """
    tree_stats_mtx = lmpy.Matrix(
        np.zeros((pam.shape[0], 6), dtype=np.float),
        headers={
            '0': pam.get_row_headers(),
            '1': ['Median Node Height', 'Node Height 75th percentile',
                  'Node Height 90th percentile', 'Median Tip Length',
                  'Tip Length 75th percentile', 'Tip Length 90th percentile']})

    squids = pam.get_column_headers()
    for site_index in range(pam.shape[0]):
        present_species = np.where(pam[site_index] == 1)[0]
        present_squids = []
        for squid_idx in present_species:
            present_squids.append(squids[squid_idx])
        subtree = get_subtree(tree, present_squids, 'squid')
        tree_stats_mtx[site_index] = calculate_tree_statistics(subtree)
    return tree_stats_mtx


# .............................................................................
def calculate_tree_statistics(tree):
    """Calculate statistic for a (sub)tree.

    Args:
        tree (TreeWrapper): A tree or subtree to calculate statistics for.

    Returns:
        tuple - Tree statistics of interest
    """
    node_heights = tree.internal_node_ages()
    tip_lengths = []
    for node in tree.nodes():
        if node.is_leaf():
            tip_lengths.append(node.edge_length)
    return (
        np.median(node_heights), np.percentile(node_heights, 75),
        np.percentile(node_heights, 90), np.median(tip_lengths),
        np.percentile(tip_lengths, 75), np.percentile(tip_lengths, 90))


# .............................................................................
class Test_calculate_continuous_ancestral_states(object):
    """Tests ancestral state reconstruction.
    """
    # .....................................
    def test_package_invalid(self, invalid_ancestral_state_package):
        """Test calculate_continusous_ancestral_states with invalid data.

        Args:
            invalid_ancestral_state_package (pytest.fixture): A parameterized
                pytest fixture defined in conftest.py that provides an invalid
                test package.

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label.

        Raises:
            IOError: When the alignment or tree cannot be loaded for the
                specified file extension.
        """
        # Get the data files
        (tree_filename, alignment_filename, results_filename
         ) = invalid_ancestral_state_package
        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise IOError(
                'Cannot handle tree with extension: {}'.format(tree_ext))
        tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)

        # Process the alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_csv_alignment_flo(
                    align_file)
        elif align_ext == '.json':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_json_alignment_flo(
                    align_file)
        elif align_ext == '.phylip':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_phylip_alignment_flo(
                    align_file)
        elif align_ext == '.tbl':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_table_alignment_flo(
                    align_file)
        else:
            raise IOError(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
                    sequences)
        # Run analysis
        with pytest.raises(Exception):
            anc_dp.calculate_continuous_ancestral_states(tree, char_mtx)

    # .....................................
    def test_package_valid(self, valid_ancestral_state_package):
        """Tests the calculate_continusous_ancestral_states method.

        Args:
            valid_ancestral_state_package (pytest.fixture): A parameterized
                pytest fixture defined in conftest.py that provides a valid
                test package.

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label.

        Raises:
            IOError: When the tree or alignment cannot be loaded for the
                specified file extension.
            Exception: When a specified successful result value cannot be
                found.
        """
        # Get the data files
        (tree_filename, alignment_filename, results_filename
         ) = valid_ancestral_state_package

        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise IOError(
                'Cannot handle tree with extension: {}'.format(tree_ext))
        # tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)
        tree = TreeWrapper.get(path=tree_filename, schema=tree_schema)

        # Process the alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_csv_alignment_flo(
                    align_file)
        elif align_ext == '.json':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_json_alignment_flo(
                    align_file)
        elif align_ext == '.phylip':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_phylip_alignment_flo(
                    align_file)
        elif align_ext == '.tbl':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_table_alignment_flo(
                    align_file)
        else:
            raise IOError(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
            sequences)
        # Run analysis
        _, anc_mtx = anc_dp.calculate_continuous_ancestral_states(
            tree, char_mtx, calc_std_err=True, sum_to_one=False)

        # New testing method
        # (For now) assume that results file is csv with row headers for
        #    node labels and column headers for variables
        results = []
        h = None
        with open(results_filename) as results_file:
            for line in results_file:
                if h is None:
                    # Get headers
                    h = line.strip().split(',')[1:]
                else:
                    # Add result (without label) to list
                    node_result = [
                        float(i) for i in line.strip().split(',')[1:]]
                    results.append(np.array(node_result, dtype=float))

        # Look for all results (only maximum likelihood)
        for row in anc_mtx[:, :, 0]:
            found = False
            for i in range(len(results)):
                # Allow for some wiggle room with decimal precision
                if np.all(np.isclose(row, results[i])):
                    found = True
                    results.pop(i)
                    break
            if not found:
                raise Exception(
                    'Could not find expected result: {} in results'.format(
                        row))


# .............................................................................
class Test_ancestral_distribution(object):
    """Tests ancestral distribution reconstruction.
    """
    # .....................................
    def test_package_invalid(self, invalid_ancestral_distribution_package):
        """Test calculate_ancestral_distributions method with invalid data.

        Args:
            invalid_ancestral_distribution_package (pytest.fixture): A pytest
                fixture that is parametrized to provide invalid ancestral
                distributions, one at a time, so that there are multiple test
                functions defined for each invalid package.

        Note:
            * This test will need to evolve as the output format changes.  It
                will probably be better to return a data structure with various
                values for each node rather than assigning the value to the
                node label.

        Raises:
            IOError: When the tree or alignment cannot be loaded for the
                specified file extension.
        """
        # Get the data files
        (tree_filename, alignment_filename, results_filename
         ) = invalid_ancestral_distribution_package
        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise IOError(
                'Cannot handle tree with extension: {}'.format(tree_ext))
        tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)

        # Process the alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_csv_alignment_flo(
                    align_file)
        elif align_ext == '.json':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_json_alignment_flo(
                    align_file)
        elif align_ext == '.phylip':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_phylip_alignment_flo(
                    align_file)
        elif align_ext == '.tbl':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_table_alignment_flo(
                    align_file)
        else:
            raise IOError(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
                    sequences)
        # Run analysis
        with pytest.raises(Exception):
            anc_dp.calculate_ancestral_distributions(tree, char_mtx)

    # .....................................
    def test_package_valid(self, valid_ancestral_distribution_package):
        """Tests the calculate_ancestral_distributions method.

        Args:
            invalid_ancestral_distribution_package (pytest.fixture): A pytest
                fixture that is parametrized to provide invalid ancestral
                distributions, one at a time, so that there are multiple test
                functions defined for each invalid package.

        Raises:
            IOError: When the tree or alignment cannot be loaded for the
                specified file extension.
            Exception: When a specified successful result value cannot be
                found.
        """
        # Get the data files
        (tree_filename, alignment_filename, results_filename
         ) = valid_ancestral_distribution_package
        # Process the tree file
        _, tree_ext = os.path.splitext(tree_filename)
        if tree_ext == '.nex':
            tree_schema = 'nexus'
        elif tree_ext == '.xml':
            tree_schema = 'nexml'
        elif tree_ext == '.tre':
            tree_schema = 'newick'
        else:
            raise IOError(
                'Cannot handle tree with extension: {}'.format(tree_ext))
        # tree = dendropy.Tree.get(path=tree_filename, schema=tree_schema)
        tree = TreeWrapper.get(path=tree_filename, schema=tree_schema)

        # Process the alignment file
        _, align_ext = os.path.splitext(alignment_filename)
        if align_ext == '.csv':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_csv_alignment_flo(
                    align_file)
        elif align_ext == '.json':
            with open(alignment_filename) as align_file:
                sequences, headers = data_readers.read_json_alignment_flo(
                    align_file)
        elif align_ext == '.phylip':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_phylip_alignment_flo(
                    align_file)
        elif align_ext == '.tbl':
            with open(alignment_filename) as align_file:
                sequences = data_readers.read_table_alignment_flo(
                    align_file)
        else:
            raise IOError(
                'Cannot handle alignments with extension: {}'.format(
                    align_ext))

        char_mtx = data_readers.get_character_matrix_from_sequences_list(
            sequences)
        # Run analysis
        _, anc_mtx = anc_dp.calculate_ancestral_distributions(
            tree, char_mtx)

        # Testing method
        # Assume that the results file is a csv with row headers for node
        #    labels and output layer (maximum_likeliehood / standard_error)
        #    and column headers for variables
        ml_results = []
        std_err_results = []
        h = None
        with open(results_filename) as results_file:
            for line in results_file:
                if h is None:
                    # Get headers
                    h = line.strip().split(',')[1:]
                else:
                    # Add result (without label) to appropriate list
                    parts = line.strip().split(',')
                    layer = parts[1].lower()
                    values = np.array(
                        [float(i) for i in parts[2:]], dtype=float)
                    if layer == 'maximum_likelihood':
                        ml_results.append(values)
                    else:
                        std_err_results.append(values)
        assert(len(ml_results) == len(std_err_results))
        print('ml results')
        print(ml_results)
        print('std err results')
        print(std_err_results)

        # Look for all results (ml and std err results should match rows)
        for row_idx in range(anc_mtx.shape[0]):
            found = False
            # Get rows from data
            ml_row = anc_mtx[row_idx, :, 0]
            std_err_row = anc_mtx[row_idx, :, 1]

            for i in range(len(ml_results)):
                print(ml_results[i])
                print(std_err_results[i])
                if np.all(np.isclose(ml_row, ml_results[i])) and \
                        np.all(np.isclose(
                            std_err_row, std_err_results[i])):
                    found = True
                    ml_results.pop(i)
                    std_err_results.pop(i)
                    break
            if not found:
                raise Exception(
                    'Could not find {}, {} in results'.format(
                        ml_row, std_err_row))
