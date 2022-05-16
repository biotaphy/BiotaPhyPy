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
    tax_info = resolve_names_otol(species_list)
    ott_ids = []
    for tax in tax_info.values():
        if tax is not None and 'ott_id' in tax.keys() and tax['ott_id'] is not None:
            ott_ids.append(tax['ott_id'])
    tree_str = induced_subtree(ott_ids)
    tree = TreeWrapper.get(data=tree_str, schema='newick')
    tree.write(path=args.tree_filename, schema=args.tree_schema)


# .....................................................................................
if __name__ == '__main__':  # pragma: no cover
    cli()
