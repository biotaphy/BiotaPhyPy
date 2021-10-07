"""Module for testing site statistic computations."""
import numpy as np

from lmpy import Matrix, TreeWrapper

from biotaphy.analyses.pam_stats.site_statistics import (
    get_tip_lengths, mean_node_height, median_node_height,
    node_height_percentile_2_5, node_height_percentile_25,
    node_height_percentile_75, node_height_percentile_97_5,
    mean_tip_length, median_tip_length, tip_length_percentile_2_5,
    tip_length_percentile_25, tip_length_percentile_75,
    tip_length_percentile_97_5, calculate_tree_site_statistics
)


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
                [3.0, 0.6, 9.0, 0.6, 0.0, 40.0, 40.0, 30.5, 30.5, 45.0, 49.5,
                    36.6666667, 30.0, 30.0, 30.0, 40.0, 49.0, 73.3333333,
                    86.6666667, 260.0, 1.0],
                [3.0, 0.6, 9.0, 0.6, 0.0, 45.0, 45.0, 40.25, 40.25, 47.5,
                    49.75, 43.3333333, 40.0, 40.0, 40.0, 45.0, 49.5,
                    86.6666667, 93.3333333, 280.0, 1.0],
                [2.0, 0.4, 6.0, 0.6, 90.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
                    40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 80.0, 80.0, 80.0, 0.0],
                [4.0, 0.8, 12.0, 0.6, 0.0, 40.0, 40.0, 30.5, 30.5, 45.0, 49.5,
                    37.5, 35.0, 30.0, 30.0, 42.5, 49.25, 75.0, 86.6666667,
                    520.0, 0.2],
                [3.0, 0.6, 9.0, 0.6, 100.0, 25.0, 25.0, 20.25, 20.25, 27.5,
                    29.75, 23.3333333, 20.0, 20.0, 20.0, 25.0, 29.5,
                    46.6666667, 53.3333333, 160.000000, 1.0]]))
        tmp = calculate_tree_site_statistics(pam, tree)
        print(tmp)
        print(test_stats)
        assert np.all(np.isclose(test_stats, tmp))


# .............................................................................
class Test_individual_stats:
    """Test individual stats."""
    # ................................
    def test_get_tip_lengths(self):
        """Simple test of get_tip_lengths."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        tip_lengths = get_tip_lengths(tree)
        assert sorted(tip_lengths) == sorted([50, 40, 30, 20, 20])

    # ................................
    def test_mean_node_height(self):
        """Simple test of mean_node_height."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert mean_node_height(tree) == 35

    # ................................
    def test_median_node_height(self):
        """Simple test of median_node_height."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert median_node_height(tree) == 35

    # ................................
    def test_node_height_percentile_2_5(self):
        """Simple test of node_height_percentile_2_5."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert node_height_percentile_2_5(tree) == 20.75

    # ................................
    def test_node_height_percentile_25(self):
        """Simple test of node_height_percentile_25."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert node_height_percentile_25(tree) == 27.5

    # ................................
    def test_node_height_percentile_75(self):
        """Simple test of node_height_percentile_75."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert node_height_percentile_75(tree) == 42.5

    # ................................
    def test_node_height_percentile_97_5(self):
        """Simple test of node_height_percentile_97_5."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert node_height_percentile_97_5(tree) == 49.25

    # ................................
    def test_mean_tip_length(self):
        """Simple test of mean_tip_length."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert mean_tip_length(tree) == 32

    # ................................
    def test_median_tip_length(self):
        """Simple test of median_tip_length."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert median_tip_length(tree) == 30

    # ................................
    def test_tip_length_percentile_2_5(self):
        """Simple test of tip_length_percentile_2_5."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert tip_length_percentile_2_5(tree) == 20

    # ................................
    def test_tip_length_percentile_25(self):
        """Simple test of tip_length_percentile_25."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert tip_length_percentile_25(tree) == 20

    # ................................
    def test_tip_length_percentile_75(self):
        """Simple test of tip_length_percentile_75."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert tip_length_percentile_75(tree) == 40

    # ................................
    def test_tip_length_percentile_97_5(self):
        """Simple test of tip_length_percentile_97_5."""
        tree_str = '(A:50,(B:40,(C:30,(D:20,E:20):10):10):10);'
        tree = TreeWrapper.get(data=tree_str, schema='newick')
        assert tip_length_percentile_97_5(tree) == 49
