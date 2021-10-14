#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module for Open Tree of Life client."""
import json
from urllib.request import Request, urlopen


PRODUCTION_SERVER = 'https://api.opentreeoflife.org/v3'
INDUCED_SUBTREE_BASE_URL = '{}/tree_of_life/induced_subtree'.format(PRODUCTION_SERVER)
OTT_TAXON_INFO_URL = '{}/tnrs/match_names'.format(PRODUCTION_SERVER)

MAX_NAMES_PER_REQUEST = 1000


# .....................................................................................
class LABEL_FORMAT:
    """Represents the label format constants used when calling induced subtree."""
    NAME = 'name'
    ID = 'id'
    NAME_AND_ID = 'name_and_id'


# .....................................................................................
def get_info_for_names(names_list):
    """Get information from the OTL taxon match service for a list of names.

    Args:
        names_list (:obj:`list` of :obj:`str`): A list of taxon names to get
            information for.

    Returns:
        dict: A dictionary where keys are the searched taxon names and the values are
            dictionaries of values from Open Tree.
        list: A list of taxa that were not found.
    """
    taxa_info = {}
    not_found_taxa = []
    headers = {'Content-Type': 'application/json'}
    for i in range(0, len(names_list), MAX_NAMES_PER_REQUEST):
        request_body = {'names': names_list[i:i+MAX_NAMES_PER_REQUEST]}
        req = Request(
            OTT_TAXON_INFO_URL,
            data=json.dumps(request_body).encode('utf8'),
            headers=headers
            )
        resp = urlopen(req)
        resp_json = json.load(resp)
        # Add taxa that we didn't find to list
        not_found_taxa.extend(resp_json['unmatched_names'])
        for result in resp_json['results']:
            taxon = result['name']
            vals = {}
            cont = True
            for match in result['matches']:
                if cont:
                    if not match['is_synonym']:
                        cont = False
                    vals['ott_id'] = match['taxon']['ott_id']
                    for tax_source in match['taxon']['tax_sources']:
                        if tax_source.startswith('gbif'):
                            vals['gbif_id'] = int(tax_source[5:])
                    vals['accepted_name'] = match['taxon']['name']
                    vals['synonyms'] = match['taxon']['synonyms']
            taxa_info[taxon] = vals
    return taxa_info, not_found_taxa


# .....................................................................................
def induced_subtree(ott_ids, label_format=LABEL_FORMAT.NAME):
    """Retrieves a Newick tree containing the nodes represented by the ids.

    Calls the Open Tree 'induced_subtree' service to retrieve a tree, in Newick
    format, containing the nodes represented by the provided Open Tree IDs.

    Args:
        ott_ids (list) : A list of Open Tree IDs.  These will be converted to
            integers in the request.
        label_format (str) : The label string format to use when creating the
            tree on the server. (see: LABEL_FORMAT)

    Returns:
        dict: A dictionary of the subtree response after JSON processing.
    """
    # Ids need to be integers
    processed_ids = [int(ottid) for ottid in ott_ids]
    request_body = {
        'ott_ids': processed_ids,
        'label_format': label_format
    }

    headers = {
        'Content-Type': 'application/json'
    }
    req = Request(
        INDUCED_SUBTREE_BASE_URL,
        data=json.dumps(request_body).encode('utf-8'),
        headers=headers)

    resp_str = urlopen(req).read().decode('utf-8')
    return json.loads(resp_str)