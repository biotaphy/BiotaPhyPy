#!/usr/bin/env python
"""Tool for performing ancestral distribution computations.

Todo:
    * Constants.
    * Clean up help.
"""
import argparse

from lmpy import TreeWrapper

from biotaphy.analyses import anc_dp
import biotaphy.common.annotators as annotators
import biotaphy.common.plots as tree_plots
from biotaphy.common import data_readers

DESCRIPTION = """\
Generates ancestral distribution estimations based on the environmental
 distributions at the tips of the tree"""


# .....................................................................................
def cli():
    """Command-line interface for the tool.

    Raises:
        ValueError: Raised if a column cannot be found for a label or bad format.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        'in_tree_filename', type=str, help='Path to the tree file')
    parser.add_argument(
        'in_tree_schema', type=str, help='The format of the tree',
        choices=['newick', 'nexml', 'nexus'])

    parser.add_argument(
        'data_filename', type=str,
        help='Path to file with character state data')
    parser.add_argument(
        'data_format', type=str, help='The format of the character data',
        choices=['csv', 'json', 'phylip', 'table'])

    # Outputs
    # Annotated tree or trees
    # Plots
    # Matrix csv
    parser.add_argument(
        'out_tree_filename', type=str,
        help='Path to write the resulting annotated tree')
    parser.add_argument(
        'out_tree_schema', type=str,
        help='The format to use when writing the tree',
        choices=['newick', 'nexml', 'nexus'])
    parser.add_argument(
        '-l', '--annotate_labels', type=str,
        help='If provided, annotate the tree labels with this data column')
    parser.add_argument(
        '-p', '--plot_directory', type=str,
        help='If provided, write distribution plots to this directory')
    parser.add_argument(
        '-c', '--out_csv_filename', type=str,
        help='If provided, write the output character matrix CSV '
             'to this file location')

    args = parser.parse_args()

    # Read the tree
    tree = TreeWrapper.get(
        path=args.in_tree_filename, schema=args.in_tree_schema)

    # Read data
    if args.data_format == 'csv':  # pragma: no cover
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_csv_alignment_flo(
                in_file)
    elif args.data_format == 'json':  # pragma: no cover
        with open(args.data_filename) as in_file:
            sequences, headers = data_readers.read_json_alignment_flo(in_file)
    elif args.data_format == 'phylip':  # pragma: no cover
        with open(args.data_filename) as in_file:
            sequences = data_readers.read_phylip_alignment_flo(in_file)
        headers = None
    elif args.data_format == 'table':
        with open(args.data_filename) as in_file:
            sequences = data_readers.read_table_alignment_flo(in_file)
        headers = None
    else:  # pragma: no cover
        raise ValueError('Unknown data format: {}'.format(args.data_format))

    # Get the label annotation column, or None
    label_column = None
    if args.annotate_labels is not None:  # pragma: no cover
        try:
            # Try looking for the string
            label_column = headers.index(args.annotate_labels)
        except Exception:
            try:
                # Treat it as an integer
                label_column = int(args.annotate_labels)
            except Exception:
                raise ValueError(
                    'Could not find column to use for labels.  '
                    'Check the name to make sure it matches or use column'
                    ' index.')

    # Get character matrix
    char_mtx = data_readers.get_character_matrix_from_sequences_list(
        sequences, var_headers=headers)
    # Run analysis
    tree, results = anc_dp.calculate_ancestral_distributions(tree, char_mtx)

    # Should we annotate the tree labels?
    if label_column is not None:  # pragma: no cover
        annotators.annotate_tree_with_label(
            tree, results, label_column=label_column)
    else:
        # Annotate tree
        annotators.add_all_annotations(tree, results, update=True)

    # Write the tree
    tree.write(path=args.out_tree_filename, schema=args.out_tree_schema)

    # CSV
    if args.out_csv_filename is not None:
        with open(args.out_csv_filename, 'w') as out_csv_f:
            results.write_csv(out_csv_f)
    # Plots
    if args.plot_directory is not None:  # pragma: no cover
        tree_plots.create_distribution_plots(
            tree, results, args.plot_directory)


# .....................................................................................
if __name__ == '__main__':  # pragma: no cover
    cli()
