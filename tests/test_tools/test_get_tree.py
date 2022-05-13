"""This module tests the get_tree script."""
import os

from lmpy.species_list import SpeciesList
from lmpy.tree import TreeWrapper

from biotaphy.tools.get_tree import cli


# .....................................................................................
def test_get_tree(monkeypatch, valid_species_list_filename, tmpdir):
    """Test the get_tree tool.

    Args:
        monkeypatch (pytest.fixture): A fixture for monkeypatching the environment.
        valid_species_list_filename (pytest.Fixture): A fixture providing a species
            list filename.
        tmpdir (pytest.fixture): A built-in pytest fixture providing a
            temporary directory for this test function.
    """
    tree_filename = os.path.join(tmpdir.dirname, 'test_tree.tre')
    tree_schema = 'newick'

    params = [
        'get_tree.py',
        valid_species_list_filename,
        tree_filename,
        tree_schema,
    ]
    # Monkeypatch the command line arguments
    monkeypatch.setattr('sys.argv', params)
    # Call the command line interface function
    cli()

    # Load tree and check that there are some tips (<= len(species_list))
    species_list = SpeciesList.from_file(valid_species_list_filename)
    tree = TreeWrapper.get(path=tree_filename, schema=tree_schema)

    assert len(tree.taxon_namespace) <= len(species_list)
