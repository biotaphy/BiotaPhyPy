"""Module containing functions for computing site statistics."""
import numpy as np

from lmpy.statistics.pam_stats import PamStats, TreeMetric


TOLERANCE = 0.005


# .............................................................................
def get_tip_lengths(tree):
    """Get the tip lengths in a tree.

    Args:
        tree (Tree): A tree object to get tip lengths from.

    Returns:
        list: A list of tip lengths in the tree.
    """
    tip_lengths = []
    for node in tree.nodes():
        if node.is_leaf():
            tip_lengths.append(node.edge_length)
    return tip_lengths


# .............................................................................
@TreeMetric
def mean_node_height(tree):
    """Get the mean node height for the tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The mean height of the nodes in the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.mean(node_heights))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def median_node_height(tree):
    """Get the median node height for the tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The median height of the nodes in the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.median(node_heights))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def node_height_percentile_2_5(tree):
    """Determine the 2.5 percentile node height for a tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The 2.5 percentile node height for the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.percentile(node_heights, 2.5))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def node_height_percentile_25(tree):
    """Determine the 25 percentile node height for a tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The 25 percentile node height for the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.percentile(node_heights, 25))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def node_height_percentile_75(tree):
    """Determine the 75 percentile node height for a tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The 75 percentile node height for the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.percentile(node_heights, 75))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def node_height_percentile_97_5(tree):
    """Determine the 97.5 percentile node height for a tree.

    Args:
        tree (Tree): A tree object used to measure node heights.

    Returns:
        float: The 97.5 percentile node height for the tree.
    """
    try:
        node_heights = tree.internal_node_ages(ultrametricity_precision=TOLERANCE)
        return float(np.percentile(node_heights, 97.5))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def mean_tip_length(tree):
    """Determine the mean tip length for the tips in a tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The mean tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.mean(tip_lengths))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def median_tip_length(tree):
    """Determine the median tip length for the tips in a tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The median tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.median(tip_lengths))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def tip_length_percentile_2_5(tree):
    """Determine the 2.5 percentile for tip length for the tips of the tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The 2.5 percentile tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.percentile(tip_lengths, 2.5))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def tip_length_percentile_25(tree):
    """Determine the 25 percentile for tip length for the tips of the tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The 25 percentile tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.percentile(tip_lengths, 25))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def tip_length_percentile_75(tree):
    """Determine the 75 percentile for tip length for the tips of the tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The 75 percentile tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.percentile(tip_lengths, 75))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
@TreeMetric
def tip_length_percentile_97_5(tree):
    """Determine the 97.5 percentile for tip length for the tips of the tree.

    Args:
        tree (Tree): A tree object used to retrieve tip lengths.

    Returns:
        float: The 97.5 percentile tip length for the tips in the tree.
    """
    try:
        tip_lengths = get_tip_lengths(tree)
        return float(np.percentile(tip_lengths, 97.5))
    except Exception:  # pragma: no cover
        return 0.0


# .............................................................................
def calculate_tree_site_statistics(pam, tree):
    """Calculate tree statistics for each site in a PAM.

    Args:
        pam (Matrix): A PAM to compute statistics for.
        tree (TreeWrapper): A corresponding tree to use for statistics.

    Returns:
        Matrix - A site by statistic matrix of tree pam statistics.
    """
    ps = PamStats(pam, tree=tree)
    reg_metrics = [
        ('Mean Node Height', mean_node_height),
        ('Median Node Height', median_node_height),
        ('Node Height 2.5th percentile', node_height_percentile_2_5),
        ('Node Height 25th percentile', node_height_percentile_2_5),
        ('Node Height 75th percentile', node_height_percentile_75),
        ('Node Height 97.5th percentile', node_height_percentile_97_5),
        ('Mean Tip Length', mean_tip_length),
        ('Median Tip Length', median_tip_length),
        ('Tip Length 2.5th percentile', tip_length_percentile_2_5),
        ('Tip Length 25th percentile', tip_length_percentile_25),
        ('Tip Length 75th percentile', tip_length_percentile_75),
        ('Tip Length 97.5th percentile', tip_length_percentile_97_5)
    ]
    for name, metric in reg_metrics:
        ps.register_metric(name, metric)
    site_stats_mtx = ps.calculate_site_statistics()
    return site_stats_mtx