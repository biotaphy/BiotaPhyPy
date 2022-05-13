"""This tool gets a tree from Open Tree from a species list."""
import argparse

from lmpy.species_list import SpeciesList
from lmpy.tree import TreeWrapper

from biotaphy.client.ot_service_wrapper.open_tree import (
    induced_subtree,
    resolve_names_otol,
)


DESCRIPTION = 'Tool to get a tree for a list of species.'


# .....................................................................................
def build_parser():
    """Get an argparse.ArgumentParser for the tool.

    Returns:
        argparse.ArgumentParser: An ArgumentParser for this tool.
    """
    parser = argparse.ArgumentParser(prog='get_tree', description=DESCRIPTION)
    parser.add_argument(
        'species_list_filename',
        type=str,
        help='A file containing a list of species.'
    )
    parser.add_argument(
        'tree_filename', type=str, help='File path to write the retrieved tree.'
    )
    parser.add_argument(
        'tree_schema', type=str, help='Schema format to write the tree.'
    )
    return parser


# .....................................................................................
def cli():
    """A command-line interface for the tool."""
    parser = build_parser()
    args = parser.parse_args()
    species_list = SpeciesList.from_file(args.species_list_filename)
    ott_ids = resolve_names_otol(species_list)
    tree_str = induced_subtree(ott_ids)
    raise Exception(tree_str)
    tree = TreeWrapper(data=tree_str, schema='newick')
    tree.write(path=args.tree_filename, schema=args.tree_schema)


# .....................................................................................
if __name__ == '__main__':
    cli()
