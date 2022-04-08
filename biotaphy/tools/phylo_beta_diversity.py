#!/usr/bin/env python
"""Tool for computing beta diversity and phylo beta diversity.

Args:
1) file location of the PAM (presence absence matrix)
2) file location of the tree associsated with PAM
3) phylo beta diversity family to use
4) file location to write output
5) (optional) the number of permuatations to perform
6) (optional) alpha value to determine significance (defaul = 0.05)

Todo:
    * Constants.
    * Clean up help.
"""
import argparse
import os

from lmpy import TreeWrapper

from biotaphy.analyses.helpers import data_readers
from biotaphy.analyses.phylo_beta_diversity import phylo_beta_diversity

DESCRIPTION = """
Computes phylogenetic & ecological beta diversity components for Sorensen and
Jaccard Indices."""


# .....................................................................................
def cli():
    """Command line interface for the tool.

    Raises:
        ValueError: Raised for unknown format or missing family.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        'in_tree_filename', type=str, help='Path to the tree file')

    parser.add_argument(
        'pam_filename', type=str,
        help='Path to file with presence/ absence data (PAM)')
    parser.add_argument(
        'data_format', type=str, help='The format of the PAM',
        choices=['csv', 'json', 'phylip', 'table'])

    # Outputs
    # Annotated tree or trees
    # Plots
    # Matrix csv

    # Might use this for nodal beta diversity

    # parser.add_argument(
    #     'out_tree_filename', type=str,
    #     help='Path to write the resulting annotated tree')
    # parser.add_argument(
    #     'out_tree_schema', type=str,
    #     help='The format to use when writing the tree',
    #     choices=['newick', 'nexml', 'nexus'])
    # parser.add_argument(
    #     '-l', '--annotate_labels', type=str,
    #     help='If provided, annotate the tree labels with this data column')
    # parser.add_argument(
    #     '-p', '--plot_directory', type=str,
    #     help='If provided, write distribution plots to this directory')

    parser.add_argument(
        'family_name', type=str,
        help='Beta diversity family metric to calculate',
        choices=['sorensen', 'jaccard'])

    parser.add_argument(
        'out_foldername', type=str,
        help='Write the output of beta diversity calculations to this folder')

    args = parser.parse_args()

    # Read data
    if args.data_format == 'csv':
        with open(args.pam_filename) as in_file:
            sequences, headers = data_readers.read_csv_alignment_flo(
                in_file)
    elif args.data_format == 'json':  # pragma: no cover
        with open(args.pam_filename) as in_file:
            sequences, headers = data_readers.read_json_alignment_flo(
                in_file)
    elif args.data_format == 'phylip':  # pragma: no cover
        with open(args.pam_filename) as in_file:
            sequences = data_readers.read_phylip_alignment_flo(in_file)
        headers = None
    elif args.data_format == 'table':  # pragma: no cover
        with open(args.pam_filename) as in_file:
            sequences = data_readers.read_table_alignment_flo(in_file)
        headers = None
    else:  # pragma: no cover
        raise ValueError('Unknown data format: {}'.format(args.data_format))

    # Convert data to PAM format
    pam = data_readers.get_character_matrix_from_sequences_list(
        sequences, headers)

    # Read the tree
    tree = TreeWrapper.from_filename(args.in_tree_filename)

    print(args.family_name)
    # Run analysis
    if args.family_name == 'jaccard':
        results = phylo_beta_diversity.calculate_phylo_beta_diversity_jaccard(
            pam, tree)
        res_names = [
            'beta_jtu', 'phylo_beta_jtu', 'beta_jne', 'phylo_beta_jne',
            'beta_jac', 'phylo_beta_jac']
    elif args.family_name == 'sorensen':
        results = phylo_beta_diversity.calculate_phylo_beta_diversity_sorensen(
            pam, tree)
        res_names = [
            'beta_sim', 'phylo_beta_sim', 'beta_sne', 'phylo_beta_sne',
            'beta_sor', 'phylo_beta_sor']
    else:  # pragma: no cover
        raise ValueError('Could not find family name')

    # Write results to folder
    if not os.path.exists(args.out_foldername):  # pragma: no cover
        os.makedirs(args.out_foldername)

    for table in range(len(results)):
        # print table, res_names[table], "\n", results[table].data, "\n"
        with open(
            os.path.join(args.out_foldername, '{}.csv'.format(
                res_names[table])), 'w') as out_csv_f:
            results[table].write_csv(out_csv_f)


# .....................................................................................
if __name__ == '__main__':  # pragma: no cover
    cli()
