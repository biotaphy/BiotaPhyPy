"""Module containing accepted name wrangler for a matrix using Open Tree of Life."""
from lmpy.data_wrangling.matrix.accepted_name_wrangler import (
    AcceptedNameMatrixWrangler
)

from biotaphy.client.ot_service_wrapper.open_tree import resolve_names_otol


# .....................................................................................
class BiotaphyAcceptedNameMatrixWrangler(AcceptedNameMatrixWrangler):
    """Modifies matrix columns to update taxon names to the "accepted" names."""
    name = 'BiotaphyAcceptedNameMatrixWrangler'
    version = '1.0'

    # .......................
    def __init__(
        self,
        name_map=None,
        name_resolver=None,
        taxon_axis=1,
        purge_failures=True,
        out_map_filename=None,
        map_write_interval=100,
        out_map_format='json',
        **params
    ):
        """Constructor for AcceptedNameMatrixModifier class.

        Args:
            name_map (dict): A map of original name to accepted name.
            name_resolver (str or Method): If provided, use this method for getting new
                accepted names.  If set to 'gbif', use GBIF name resolution.
            taxon_axis (int): The axis with taxon headers.
            purge_failures (bool): Should failures be purged from the matrix.
            out_map_filename (str): A file location to write the updated name map.
            map_write_interval (int): Update the name map output file after each set of
                this many iterations.
            out_map_format (str): The format to write the names map (csv or json).
            **params (dict): Keyword parameters to pass to _MatrixDataWrangler.
        """
        if isinstance(name_resolver, str) and name_resolver.lower() == 'otol':
            name_resolver = resolve_names_otol
        AcceptedNameMatrixWrangler.__init__(
            self,
            name_map=name_map,
            name_resolver=name_resolver,
            taxon_axis=taxon_axis,
            purge_failures=purge_failures,
            out_map_filename=out_map_filename,
            map_write_interval=map_write_interval,
            out_map_format=out_map_format,
            **params,
        )
