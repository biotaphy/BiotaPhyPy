#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Test module for Open Tree of Life client."""
import dendropy

from biotaphy.client.ot_service_wrapper import open_tree


# Testing constants
# ....................................
GOOD_OTT_IDS = [292466, 267845, 316878, 102710]
BAD_OTT_IDS = [999999999, 8888888888, 7777777777]

SANITIZE_NAMES_MAP = [
    (
        'Caspiomyzon hellenicus (Vladykov, Renaud, Kott & Economidis, 1982)',
        'Caspiomyzon hellenicus'
    ),
    ('Caspiomyzon wagneri (Kessler, 1870)', 'Caspiomyzon wagneri'),
    ('Entosphenus folletti Vladykov & Kott, 1976', 'Entosphenus folletti'),
    ('Eudontomyzon danfordi Regan, 1911', 'Eudontomyzon danfordi'),
    ('Eudontomyzon vladykovi Oliva & Zanandrea, 1959', 'Eudontomyzon vladykovi'),
    (
        'Lampetra alavariensis Mateus, Alves, Quintella & Almeida, 2013',
        'Lampetra alavariensis'
    ),
    ('Lampetra hubbsi (Vladykov & Kott, 1976)', 'Lampetra hubbsi'),
    ('Barbus cf. holotaenia', 'Barbus cf. holotaenia'),
    ('Heuchera ×brizoides', 'Heuchera ×brizoides'),
    ('Heuchera × brizoides', 'Heuchera × brizoides'),
    ('Heuchera abramsii', 'Heuchera abramsii'),
    ('Bacillus subtilis subsp. spizizenii', 'Bacillus subtilis subsp. spizizenii'),
    ('Somegenus somespecies var. somevar', 'Somegenus somespecies var. somevar'),
    ('Canis sp.', 'Canis sp.'),
    ('Canis spp.', 'Canis spp.'),
    ('Canis sp. A', 'Canis sp. A'),
    ('Canis spp. A, B', 'Canis spp. A, B'),
    ('Canis sp. A, B', 'Canis sp. A, B'),
    ('Canis sp. 1', 'Canis sp. 1'),
    ('Canis spp. 1, 2', 'Canis spp. 1, 2'),
    ('Canis sp. 1, 2', 'Canis sp. 1, 2'),
    ('Canis', 'Canis')
]

TAXON_NAMES = [
    'Heuchera abramsii',
    'Heuchera pilosissima',
    'Acer rubrum',
    'Boops boops'
]

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


# .....................................................................................
class Test_get_tree_from_taxa:
    """Tests that a tree can be retrieved starting with a set of taxon names."""
    # ...........................
    def test_good_taxa(self):
        """Tests that a tree can be retrieved starting with good taxon names."""
        taxa_info = open_tree.resolve_names_otol(TAXON_NAMES)
        unmatched_names = [tax for tax in taxa_info.keys() if taxa_info[tax] is None]
        assert len(unmatched_names) == 0
        ott_ids = [tax['ott_id'] for tax in taxa_info.values()]
        assert len(ott_ids) == len(TAXON_NAMES)
        resp = open_tree.induced_subtree(
            ott_ids, label_format=open_tree.LABEL_FORMAT.ID
        )
        newick = resp['newick']

        tree = dendropy.Tree.get(data=newick, schema='newick')

        assert len(tree.taxon_namespace) == len(ott_ids)

        for taxon in tree.taxon_namespace:
            # Id map contains integers as of now, check to make sure each tip
            #     is in mapping
            # Taxon labels look like 'ott{ottid}'
            search_label = taxon.label.replace('ott', '')
            assert int(search_label) in ott_ids


# .....................................................................................
class Test_sanitize_names:
    """Test that name strings are sanitized correctly."""
    # .........................
    def test_sanitation(self):
        """Test name sanitation."""
        for in_name, out_name in SANITIZE_NAMES_MAP:
            assert open_tree.sanitize_name(in_name) == out_name


# .....................................................................................
class Test_resolve_names_otol:
    """Test getting info for names."""
    # .........................
    def test_mixed(self):
        """Test good and bad taxon names."""
        bad_taxa = ['BBBBBAAAAADDDD', 'Bad taxon', 'Doesnot exist']
        test_names = bad_taxa + TAXON_NAMES
        tax_info = open_tree.resolve_names_otol(test_names)
        unmatched_names = [tax for tax in tax_info.keys() if tax_info[tax] is None]
        assert len(unmatched_names) == len(bad_taxa)
        for tax in bad_taxa:
            assert tax in unmatched_names


# .............................................................................
# class Test_induced_subtree(object):
#    """Test that the induced subtree service returns the tree we expect."""
#    # ...........................
#    def test_all_bad_ottids(self):
#        """Test that the service responds appropriately when given bad data."""
#        resp = open_tree.induced_subtree(
#            BAD_OTT_IDS, label_format=open_tree.LABEL_FORMAT.ID)
#        newick = resp['newick']
#        # Should fail if bad newick
#        tree = dendropy.Tree.get(data=newick, schema='newick')
#        # Make sure tree is not None
#        assert tree is not None
#
#    # ...........................
#    def test_all_good_ottids(self):
#        """Test that the service responds correctly when only good ids are used."""
#        resp = open_tree.induced_subtree(
#            GOOD_OTT_IDS, label_format=open_tree.LABEL_FORMAT.ID)
#        newick = resp['newick']
#        # Should fail if bad newick
#        tree = dendropy.Tree.get(data=newick, schema='newick')
#        # Make sure tree is not None
#        assert tree is not None
#
#    # ...........................
#    def test_some_bad_ottids(self):
#        """Test that the service handles a mix of good and bad ott ids."""
#        test_ids = GOOD_OTT_IDS
#        test_ids.extend(BAD_OTT_IDS)
#
#        resp = open_tree.induced_subtree(
#            test_ids, label_format=open_tree.LABEL_FORMAT.ID)
#        newick = resp['newick']
#        # Should fail if bad newick
#        tree = dendropy.Tree.get(data=newick, schema='newick')
#        # Make sure tree is not None
#        assert tree is not None
