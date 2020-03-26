"""Module for testing site statistic computations.

Notes:
    * Uses pytest style testing.
"""
import os

import numpy as np
import pytest

from lmpy import Matrix, TreeWrapper

from biotaphy.analyses.pam_stats.site_statistics import (
    calculate_tree_site_statistics, calculate_tree_statistics, get_subtree)


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
                [40.0, 40.0, 30.5, 35.0, 45.0, 49.5,
                    36.66666667, 30.0, 30.0, 30.0, 40.0, 49.0],
                [45.0, 45.0, 40.25, 42.5, 47.5, 49.75,
                    43.33333333, 40.0, 40.0, 40.0, 45.0, 49.5],
                [30.0, 30.0, 20.5, 25.0, 35.0, 39.5,
                    27.5, 25.0, 20.0, 20.0, 32.5, 39.25],
                [35.0, 35.0, 20.75, 27.5, 42.5, 49.25,
                    32.0, 30.0, 20.0, 20.0, 40.0, 49.0],
                [35.0, 35.0, 20.75, 27.5, 42.5, 49.25,
                    32.0, 30.0, 20.0, 20.0, 40.0, 49.0]
                ]))
        tmp = calculate_tree_site_statistics(pam, tree)
        print(tmp)
        print(test_stats)
        assert np.all(np.isclose(test_stats, tmp))


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
            32.0, 30.0, 20.0, 20.0, 40.0, 49.0)


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
        print(kept_labels)
        print(keep_tips)
        assert len(kept_labels) == len(keep_tips)
        for tip in keep_tips:
            assert tip[0] in kept_labels
