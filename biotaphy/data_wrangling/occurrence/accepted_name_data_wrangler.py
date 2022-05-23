"""Accepted name data wrangler for occurrences using OTOL name resolution."""
from lmpy.data_wrangling.occurrence.accepted_name_wrangler import (
    AcceptedNameOccurrenceWrangler,
)

from biotaphy.client.ot_service_wrapper.open_tree import resolve_names_otol


# .....................................................................................
class BiotaphyAcceptedNameOccurrenceWrangler(AcceptedNameOccurrenceWrangler):
    """Modifies the species_name to the "accepted" taxon name for the species."""
    name = 'BiotaphyAcceptedNameOccurrenceWrangler'
    version = '1.0'

    # .......................
    def __init__(
        self,
        name_map=None,
        name_resolver=None,
        store_original_attribute=None,
        out_map_filename=None,
        map_write_interval=100,
        out_map_format='json',
        **params
    ):
        """Constructor for BiotaphyAcceptedNameOccurrenceWrangler.

        Args:
            name_map (dict or str): A map of original name to accepted name.
            name_resolver (str or Method): If provided, use this method for getting new
                accepted names.  If set to 'gbif', use GBIF name resolution.
            store_original_attribute (str or None): A new attribute to store the
                original taxon name.
            out_map_filename (str): A file location to write the updated name map.
            map_write_interval (int): Update the name map output file after each set of
                this many iterations.
            out_map_format (str): The format to write the names map (csv or json).
            **params (dict): Keyword parameters to pass to _OccurrenceDataWrangler.
        """
        if isinstance(name_resolver, str) and name_resolver.lower() == 'otol':
            name_resolver = resolve_names_otol
        AcceptedNameOccurrenceWrangler.__init__(
            self,
            name_map=name_map,
            name_resolver=name_resolver,
            store_original_attribute=store_original_attribute,
            out_map_filename=out_map_filename,
            map_write_interval=map_write_interval,
            out_map_format=out_map_format,
            **params
        )
