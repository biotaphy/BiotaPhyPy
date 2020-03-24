"""Module containing functions for computing site statistics."""
from copy import copy

import numpy as np

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
def get_subtree(tree, present_squids, squid_attribute):
    """Get a subtree that only includes the provided taxa.

    Args:
        tree (TreeWrapper): The original tree.
        present_squids (:obj:`list` of :obj:`str`): A list of squids to include
            in the subtree.

    Returns:
        TreeWrapper: A subtree that only includes the tips specified by the
            squids.
    """
    subtree = copy(tree)
    prune_taxa = []
    for taxon in subtree.taxon_namespace:
        squid = taxon.annotations.get_value(squid_attribute)
        if val is None or val not in present_squids:
            prune_taxa.append(taxon)
    subtree.prune_taxa(prune_taxa)
    subtree.purge_taxon_namespace()
    return subtree
