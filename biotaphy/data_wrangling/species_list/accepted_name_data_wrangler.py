"""Accepted name wrangler for species list using OTOL names."""
from lmpy.data_wrangling.species_list.accepted_name_wrangler import (
    AcceptedNameSpeciesListWrangler
)

from biotaphy.client.ot_service_wrapper.open_tree import resolve_names_otol


# .....................................................................................
class BiotaphyAcceptedNameSpeciesListWrangler(AcceptedNameSpeciesListWrangler):
    """Modifies a species list to update taxon names to the "accepted" names."""
    name = 'BiotaphyAcceptedNameSpeciesListWrangler'
    version = '1.0'

    # .......................
    def __init__(
        self,
        name_map=None,
        name_resolver=None,
        purge_failures=True,
        out_map_filename=None,
        map_write_interval=100,
        out_map_format='json',
        **params
    ):
        """Constructor for BiotaphyAcceptedNameSpeciesListModifier class.

        Args:
            name_map (dict): A map of original name to accepted name.
            name_resolver (str or Method): If provided, use this method for getting new
                accepted names.  If set to 'gbif', use GBIF name resolution.
            purge_failures (bool): Should failures be purged from the tree.
            out_map_filename (str): A file location to write the updated name map.
            map_write_interval (int): Update the name map output file after each set of
                this many iterations.
            out_map_format (str): The format to write the names map (csv or json).
            **params (dict): Keyword parameters to pass to _TreeDataWrangler.
        """
        if isinstance(name_resolver, str) and name_resolver.lower() == 'otol':
            name_resolver = resolve_names_otol
        AcceptedNameSpeciesListWrangler.__init__(
            self,
            name_map=name_map,
            name_resolver=name_resolver,
            purge_failures=purge_failures,
            out_map_filename=out_map_filename,
            map_write_interval=map_write_interval,
            out_map_format=out_map_format,
            **params,
        )
