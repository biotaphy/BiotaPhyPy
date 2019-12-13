#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module for Open Tree of Life client
"""
import json
try:  # pragma: no cover
    # Python 2
    from urllib2 import Request, urlopen
except:  # pragma: no cover
    # Python 3
    from urllib.request import Request, urlopen

DEV_SERVER = 'http://141.211.236.35:10999'
INDUCED_SUBTREE_BASE_URL = '{}/induced_subtree'.format(DEV_SERVER)
OTTIDS_FROM_GBIFIDS_URL = '{}/ottids_from_gbifids'.format(DEV_SERVER)


# .............................................................................
class LABEL_FORMAT(object):
    """Represents the label format constants used when calling induced subtree
    """
    NAME = 'name'
    ID = 'id'
    NAME_AND_ID = 'name_and_id'


# .............................................................................
def get_ottids_from_gbifids(gbif_ids):
    """Retrieves a GBIF ID : OTT_ID mapping dictionary

    Calls the Open Tree 'ottids_from_gbifids' service to retrieve a mapping
    dictionary from the Open Tree service where each key is one of the provided
    GBIF identifiers and the value is the corresponding OpenTree id.

    Args:
        gbif_ids (list) : A list of GBIF taxon ids.  They will be converted to
            integers in the request

    Note:
        * Any GBIF taxon id that was not found will have a value of None

    Returns:
        Mapping dictionary
    """
    # Ids need to be integers
    processed_ids = [int(gid) for gid in gbif_ids]

    request_body = {
        'gbif_ids': processed_ids
    }

    headers = {
        'Content-Type': 'application/json'
    }
    req = Request(
        OTTIDS_FROM_GBIFIDS_URL, data=json.dumps(request_body).encode('utf-8'),
        headers=headers)

    # Note: This is done for those versions of Python 3 where urlopen requires
    #    bytes and json can't handle bytes.  Could be changed if support for
    #    Python 3.4 and 3.5 is dropped
    resp_str = urlopen(req).read().decode('utf-8')
    resp = json.loads(resp_str)
    unmatchedIds = resp['unmatched_gbif_ids']

    id_map = resp['gbif_ott_id_map']

    for gid in unmatchedIds:
        id_map[gid] = None

    return id_map


# .............................................................................
def induced_subtree(ott_ids, label_format=LABEL_FORMAT.NAME):
    """Retrieves a Newick tree containing the nodes represented by the ids

    Calls the Open Tree 'induced_subtree' service to retrieve a tree, in Newick
    format, containing the nodes represented by the provided Open Tree IDs.

    Args:
        ott_ids (list) : A list of Open Tree IDs.  These will be converted to
            integers in the request.
        label_format (str) : The label string format to use when creating the
            tree on the server. (see: LABEL_FORMAT)
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
