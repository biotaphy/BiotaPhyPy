"""Test class instantiation via WranglerFactory."""
from lmpy.data_wrangling.factory import WranglerFactory

from biotaphy.data_wrangling.matrix.accepted_name_data_wrangler import (
    BiotaphyAcceptedNameMatrixWrangler
)
from biotaphy.data_wrangling.occurrence.accepted_name_data_wrangler import (
    BiotaphyAcceptedNameOccurrenceWrangler
)
from biotaphy.data_wrangling.species_list.accepted_name_data_wrangler import (
    BiotaphyAcceptedNameSpeciesListWrangler
)
from biotaphy.data_wrangling.tree.accepted_name_data_wrangler import (
    BiotaphyAcceptedNameTreeWrangler
)


# .....................................................................................
def test_factory_instantiation_matrix():
    """Test that the Biotaphy matrix accepted name wrangler can be created."""
    wrangler_configs = [
        {
            'module': 'biotaphy.data_wrangling.matrix.accepted_name_data_wrangler',
            'wrangler_type': 'BiotaphyAcceptedNameMatrixWrangler',
        }
    ]
    factory = WranglerFactory()
    wranglers = factory.get_wranglers(wrangler_configs)
    assert len(wranglers) == 1
    assert isinstance(wranglers[0], BiotaphyAcceptedNameMatrixWrangler)


# .....................................................................................
def test_factory_instantiation_occurrence():
    """Test that the Biotaphy occurrence accepted name wrangler can be created."""
    wrangler_configs = [
        {
            'module': 'biotaphy.data_wrangling.occurrence.accepted_name_data_wrangler',
            'wrangler_type': 'BiotaphyAcceptedNameOccurrenceWrangler',
        }
    ]
    factory = WranglerFactory()
    wranglers = factory.get_wranglers(wrangler_configs)
    assert len(wranglers) == 1
    assert isinstance(wranglers[0], BiotaphyAcceptedNameOccurrenceWrangler)


# .....................................................................................
def test_factory_instantiation_species_list():
    """Test that the Biotaphy species list accepted name wrangler can be created."""
    wrangler_configs = [
        {
            'module': (
                'biotaphy.data_wrangling.species_list.'
                'accepted_name_data_wrangler'
            ),
            'wrangler_type': 'BiotaphyAcceptedNameSpeciesListWrangler',
        }
    ]
    factory = WranglerFactory()
    wranglers = factory.get_wranglers(wrangler_configs)
    assert len(wranglers) == 1
    assert isinstance(wranglers[0], BiotaphyAcceptedNameSpeciesListWrangler)


# .....................................................................................
def test_factory_instantiation_tree():
    """Test that the Biotaphy tree accepted name wrangler can be created."""
    wrangler_configs = [
        {
            'module': 'biotaphy.data_wrangling.tree.accepted_name_data_wrangler',
            'wrangler_type': 'BiotaphyAcceptedNameTreeWrangler',
        }
    ]
    factory = WranglerFactory()
    wranglers = factory.get_wranglers(wrangler_configs)
    assert len(wranglers) == 1
    assert isinstance(wranglers[0], BiotaphyAcceptedNameTreeWrangler)
