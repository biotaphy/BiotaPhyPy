#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Test module for Open Tree of Life client

Note:
     * Assumes that dendropy is available for validating trees
     * Uses pytest style testing
"""
import dendropy
import pytest

from ot_service_wrapper import open_tree

# Testing constants
# ....................................
GOOD_OTT_IDS = [292466, 267845, 316878, 102710]
BAD_OTT_IDS = [999999999, 8888888888, 7777777777]

GBIF_ID_MAP = {
    '3032647': 3943087,  # Heuchera abramsii
    '3032648': 927401,   # Heuchera pilosissima
    '3032649': 3943094,  # Heuchera caespitosa
    '3032651': 3943092,  # Heuchera pulchella
    '3032652': 3943090,  # Heuchera duranii
    '3032653': 404246,   # Heuchera sanguinea
    '3032654': 98721,    # Heuchera glabra
    '3032655': 128615,   # Heuchera villosa
    '3032656': 3943088,  # Heuchera caroliniana
    '3032658': 590982,   # Heuchera grossulariifolia
    '3032660': 3943082,  # Heuchera alba
    '3032661': 3943081,  # Heuchera wootonii
    '3032662': 927418,   # Heuchera merriamii
    '3032664': 3943085,  # Heuchera parishii
    '3032665': 883164,   # Heuchera americana
    '3032666': 927422,   # Heuchera parvifolia
    '3032667': 3943083,  # Heuchera brevistaminea
    '3032668': 3943075,  # Heuchera pubescens
    '3032670': 98712,    # Heuchera cylindrica
    '3032671': 941383,   # Heuchera rubescens
    '3032672': 927426,   # Heuchera elegans
    '3032673': 890096,   # Heuchera hirsutissima
    '3032674': 3943080,  # Heuchera flabellifolia
    '3032675': 883158,   # Heuchera maxima
    '3032676': 3943078,  # Heuchera novamexicana
    '3032678': 3943084,  # Heuchera parviflora
    '3032679': 3943086,  # Heuchera bracteata
    '3032680': 883156,   # Heuchera chlorantha
    '3032681': 927404,   # Heuchera richardsonii
    '3032686': 3943095,  # Heuchera eastwoodiae
    '3032687': 3943076,  # Heuchera glomerulata
    '3032688': 3943091,  # Heuchera hallii
    '3032689': 3943093,  # Heuchera easthamii
    '3032690': 3943089,  # Heuchera longiflora
    '3752543': 6044332,  # Heuchera hemsleyana
    '3752610': 6044331,  # Heuchera halstedii
    '3753319': 6044329,  # Heuchera amoena
    '3753512': 5729581,  # Heuchera acutifolia
    '3754294': 6044339,  # Heuchera townsendii
    '3754395': 6044338,  # Heuchera sitgreavesii
    '3754671': 6044342,  # Heuchera ×rosea
    '3754743': 6044337,  # Heuchera reglensis
    '3755291': 6044334,  # Heuchera minutiflora
    '3755546': 5729586,  # Heuchera longipetala
    '4926214': 5729582,  # Heuchera woodsiaphila
    '7462054': 5729590,  # Heuchera soltisii
    '7516328': 5729585,  # Heuchera wellsiae
    '7551031': 5729584,  # Heuchera lakelae
    '7554971': 5729591,  # Heuchera mexicana
    '7588669': 5729589,  # Heuchera rosendahlii
    '8109411': 5729592,  # Heuchera inconstans
    '8280496': 6044341,  # Heuchera ×brizoides
    '8365087': 1035572   # Heuchera micrantha
}
BAD_GBIF_IDS = ['999999', '1231312222222', '8888888']


# .............................................................................
class Test_get_tree_from_gbif_ids(object):
    """Tests that a tree can be retrieved starting with a set of GBIF ids
    """
    # ...........................
    def test_good_gbif_ids(self):
        """Tests that tree can be retrieved when input is set of GBIF ids
        """
        gbif_ids = GBIF_ID_MAP.keys()
        id_map = open_tree.get_ottids_from_gbifids(gbif_ids)
        # Make sure we haven't lost any ids
        assert len(gbif_ids) == len(id_map.keys())

        resp = open_tree.induced_subtree(
            id_map.keys(), label_format=open_tree.LABEL_FORMAT.ID)

        newick = resp['newick']

        tree = dendropy.Tree.get(data=newick, schema='newick')

        # Number of labels in the tree plus the number of unmatched should
        #     be less than or equal the number of keys in id_map.  Less than if
        #     one of the ids matches somewhere other than a tiop

        assert len(tree.taxon_namespace) + len(resp['unmatched_ott_ids']
                                               ) <= len(id_map.keys())

        for taxon in tree.taxon_namespace:
            # Id map contains integers as of now, check to make sure each tip
            #     is in mapping
            # Taxon labels look like 'ott{ottid}'
            search_label = taxon.label.strip('ott')
            assert search_label in id_map.keys()

        # Make sure any unmatched ids are in mapping
        for unmatched in resp['unmatched_ott_ids']:
            assert str(unmatched) in id_map.keys()


# .............................................................................
class Test_get_ottids_from_gbif_ids(object):
    """Test the Open Tree of Life function for retrieving ottids from GBIF ids
    """
    # ...........................
    def test_all_bad_gbif_ids(self):
        """Test that the service responds correctly to "bad" GBIF ids.

        Note:
            * The wrapper will respond by setting the value of the dictionary
                entry is None
        """
        id_map = open_tree.get_ottids_from_gbifids(BAD_GBIF_IDS)
        for gid, ottid in id_map.items():
            # Value should be None
            assert not ottid

    # ...........................
    def test_all_good_gbif_ids(self):
        """Test that he service responds correctly and values match test values
        """
        test_gbif_ids = GBIF_ID_MAP.keys()

        id_map = open_tree.get_ottids_from_gbifids(test_gbif_ids)

        # Test that each of the returned ott ids are the same as in map
        for gid, ottid in id_map.items():
            assert ottid == GBIF_ID_MAP[gid]

    # ...........................
    def test_some_bad_gbif_ids(self):
        """Test that the service responds correctly to a mix of inputs.
        """
        test_gbif_ids = list(GBIF_ID_MAP.keys())
        test_gbif_ids.extend(BAD_GBIF_IDS)

        id_map = open_tree.get_ottids_from_gbifids(test_gbif_ids)

        # Test that each of the returned ott ids are the same as in map
        for gid, ottid in id_map.items():
            if str(gid) in BAD_GBIF_IDS:
                assert not ottid
            else:
                assert ottid == GBIF_ID_MAP[gid]


# .............................................................................
class Test_induced_subtree(object):
    """Test that the induced subtree service returns the tree we expect.
    """
    # ...........................
    def test_all_bad_ottids(self):
        """Test that the service responds appropriately when given bad data
        """
        resp = open_tree.induced_subtree(
            BAD_OTT_IDS, label_format=open_tree.LABEL_FORMAT.ID)
        newick = resp['newick']
        # Should fail if bad newick
        tree = dendropy.Tree.get(data=newick, schema='newick')
        # Make sure tree is not None
        assert tree is not None

    # ...........................
    def test_all_good_ottids(self):
        """Test that the service responds correctly when only good ids are used
        """
        resp = open_tree.induced_subtree(
            GOOD_OTT_IDS, label_format=open_tree.LABEL_FORMAT.ID)
        newick = resp['newick']
        # Should fail if bad newick
        tree = dendropy.Tree.get(data=newick, schema='newick')
        # Make sure tree is not None
        assert tree is not None

    # ...........................
    def test_some_bad_ottids(self):
        """Test that the service handles a mix of good and bad ott ids
        """
        test_ids = GOOD_OTT_IDS
        test_ids.extend(BAD_OTT_IDS)

        resp = open_tree.induced_subtree(
            test_ids, label_format=open_tree.LABEL_FORMAT.ID)
        newick = resp['newick']
        # Should fail if bad newick
        tree = dendropy.Tree.get(data=newick, schema='newick')
        # Make sure tree is not None
        assert tree is not None
