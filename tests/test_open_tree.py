#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@summary: Test module for Open Tree of Life client
@author: CJ Grady
@version: 1.0.0
@status: beta

@license: gpl3
@copyright: Copyright (C) 2018, University of Kansas Center for Research

          Lifemapper Project, lifemapper [at] ku [dot] edu, 
          Biodiversity Institute,
          1345 Jayhawk Boulevard, Lawrence, Kansas, 66045, USA
   
          This program is free software; you can redistribute it and/or modify 
          it under the terms of the GNU General Public License as published by 
          the Free Software Foundation; either version 2 of the License, or (at 
          your option) any later version.
  
          This program is distributed in the hope that it will be useful, but 
          WITHOUT ANY WARRANTY; without even the implied warranty of 
          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
          General Public License for more details.
  
          You should have received a copy of the GNU General Public License 
          along with this program; if not, write to the Free Software 
          Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
          02110-1301, USA.
          
@note: Assumes that dendropy is available for validating trees
"""
import dendropy
import logging
import unittest

from open_tree import get_ottids_from_gbifids, induced_subtree, LABEL_FORMAT

# Testing constants
# ....................................
GOOD_OTT_IDS = [292466, 267845, 316878, 102710]
BAD_OTT_IDS = [999999999, 8888888888, 7777777777]

GBIF_ID_MAP = {
   "3032647" : 3943087, # Heuchera abramsii
   "3032648" : 927401,  # Heuchera pilosissima
   "3032649" : 3943094, # Heuchera caespitosa
   "3032651" : 3943092, # Heuchera pulchella
   "3032652" : 3943090, # Heuchera duranii
   "3032653" : 404246,  # Heuchera sanguinea
   "3032654" : 98721,   # Heuchera glabra
   "3032655" : 128615,  # Heuchera villosa
   "3032656" : 3943088, # Heuchera caroliniana
   "3032658" : 590982,  # Heuchera grossulariifolia
   "3032660" : 3943082, # Heuchera alba
   "3032661" : 3943081, # Heuchera wootonii
   "3032662" : 927418,  # Heuchera merriamii
   "3032664" : 3943085, # Heuchera parishii
   "3032665" : 883164,  # Heuchera americana
   "3032666" : 927422,  # Heuchera parvifolia
   "3032667" : 3943083, # Heuchera brevistaminea
   "3032668" : 3943075, # Heuchera pubescens
   "3032670" : 98712,   # Heuchera cylindrica
   "3032671" : 941383,  # Heuchera rubescens
   "3032672" : 927426,  # Heuchera elegans
   "3032673" : 890096,  # Heuchera hirsutissima
   "3032674" : 3943080, # Heuchera flabellifolia
   "3032675" : 883158,  # Heuchera maxima
   "3032676" : 3943078, # Heuchera novamexicana
   "3032678" : 3943084, # Heuchera parviflora
   "3032679" : 3943086, # Heuchera bracteata
   "3032680" : 883156,  # Heuchera chlorantha
   "3032681" : 927404,  # Heuchera richardsonii
   "3032686" : 3943095, # Heuchera eastwoodiae
   "3032687" : 3943076, # Heuchera glomerulata
   "3032688" : 3943091, # Heuchera hallii
   "3032689" : 3943093, # Heuchera easthamii
   "3032690" : 3943089, # Heuchera longiflora
   "3752543" : 6044332, # Heuchera hemsleyana
   "3752610" : 6044331, # Heuchera halstedii
   "3753319" : 6044329, # Heuchera amoena
   "3753512" : 5729581, # Heuchera acutifolia
   "3754294" : 6044339, # Heuchera townsendii
   "3754395" : 6044338, # Heuchera sitgreavesii
   "3754671" : 6044342, # Heuchera ×rosea
   "3754743" : 6044337, # Heuchera reglensis
   "3755291" : 6044334, # Heuchera minutiflora
   "3755546" : 5729586, # Heuchera longipetala
   "4926214" : 5729582, # Heuchera woodsiaphila
   "7462054" : 5729590, # Heuchera soltisii
   "7516328" : 5729585, # Heuchera wellsiae
   "7551031" : 5729584, # Heuchera lakelae
   "7554971" : 5729591, # Heuchera mexicana
   "7588669" : 5729589, # Heuchera rosendahlii
   "8109411" : 5729592, # Heuchera inconstans
   "8280496" : 6044341, # Heuchera ×brizoides
   "8365087" : 1035572  # Heuchera micrantha
}
BAD_GBIF_IDS = ["999999", "1231312222222", "8888888"]

# .............................................................................
class TestGetTreeFromGBIFIds(unittest.TestCase):
   """
   @summary: Tests that a tree can be retrieved starting with a set of GBIF ids
   """
   # ...........................
   def test_good_gbif_ids(self):
      """
      @summary: Combination test that makes sure a tree can be retrieved when
                   the input is a set of GBIF ids
      """
      gbif_ids = GBIF_ID_MAP.keys()
      id_map = get_ottids_from_gbifids(gbif_ids)
      # Make sure we haven't lost any ids
      self.assertEqual(len(gbif_ids), len(id_map.keys()))
      
      resp = induced_subtree(id_map.keys(), label_format=LABEL_FORMAT.ID)
      
      newick = resp['newick']
      
      tree = dendropy.Tree.get(data=newick, schema='newick')
      
      # Number of labels in the tree plus the number of unmatched should
      #    be less than or equal the number of keys in id_map.  Less than if
      #    one of the ids matches somewhere other than a tiop
      
      self.assertLessEqual(
         len(tree.taxon_namespace) + len(resp['unmatched_ott_ids']), 
         len(id_map.keys()))
      
      for taxon in tree.taxon_namespace:
         # Id map contains integers as of now, check to make sure each tip
         #    is in mapping
         # Taxon labels look like 'ott{ottid}'
         search_label = taxon.label.strip('ott')
         self.assertTrue(id_map.has_key(search_label))
         
      # Make sure any unmatched ids are in mapping
      for unmatched in resp['unmatched_ott_ids']:
         self.assertTrue(id_map.has_key(str(unmatched)))
      
# .............................................................................
class TestGetOTTIdsFromGBIFIds(unittest.TestCase):
   """
   @summary: Test the Open Tree of Life function for retrieving ottids from 
                GBIF ids
   """
   # ...........................
   def test_all_bad_gbif_ids(self):
      """
      @summary: Test that the service responds correctly when all of the GBIF
                   ids provided are "bad"
      @note: The wrapper will respond by setting the value of the dictionary
                entry is None
      """
      id_map = get_ottids_from_gbifids(BAD_GBIF_IDS)
      for gid, ottid in id_map.iteritems():
         # Value should be None
         self.assertFalse(ottid)
   
   # ...........................
   def test_all_good_gbif_ids(self):
      """
      @summary: Test that the service responds correctly and the returned ott 
                   ids match the test values
      """
      test_gbif_ids = GBIF_ID_MAP.keys()

      id_map = get_ottids_from_gbifids(test_gbif_ids)
      
      # Test that each of the returned ott ids are the same as the ones in map
      for gid, ottid in id_map.iteritems():
         self.assertEqual(ottid, GBIF_ID_MAP[gid])
   
   # ...........................
   def test_some_bad_gbif_ids(self):
      """
      @summary: Test that the service responds correctly when provided with 
                   good and bad gbif ids.  Good ids should return expected ott
                   ids and bad gbif ids return None
      """
      test_gbif_ids = GBIF_ID_MAP.keys()
      test_gbif_ids.extend(BAD_GBIF_IDS)

      id_map = get_ottids_from_gbifids(test_gbif_ids)
      
      # Test that each of the returned ott ids are the same as the ones in map
      for gid, ottid in id_map.iteritems():
         if str(gid) in BAD_GBIF_IDS:
            self.assertFalse(ottid)
         else:
            self.assertEqual(ottid, GBIF_ID_MAP[gid])

# .............................................................................
class TestInducedSubtree(unittest.TestCase):
   """
   @summary: Test that the Open Tree of Life Induced Subtree service returns 
                the tree we expect
   """
   # ...........................
   def test_all_bad_ottids(self):
      """
      @summary: Test that the service responds appropriately when given bad 
                   data
      """
      resp = induced_subtree(BAD_OTT_IDS, label_format=LABEL_FORMAT.ID)
      newick = resp['newick']
      # Should fail if bad newick
      tree = dendropy.Tree.get(data=newick, schema='newick')
      # Make sure tree is not None
      self.assertTrue(tree)
      
   # ...........................
   def test_all_good_ottids(self):
      """
      @summary: Test that the service responds appropriately when only good ids
                   are used
      """
      resp = induced_subtree(GOOD_OTT_IDS, label_format=LABEL_FORMAT.ID)
      newick = resp['newick']
      # Should fail if bad newick
      tree = dendropy.Tree.get(data=newick, schema='newick')
      # Make sure tree is not None
      self.assertTrue(tree)
      
   # ...........................
   def test_some_bad_ottids(self):
      """
      @summary: Test that the service responds appropriately when good and bad
                   ott ids are provided
      """
      test_ids = GOOD_OTT_IDS
      test_ids.extend(BAD_OTT_IDS)
      
      resp = induced_subtree(test_ids, label_format=LABEL_FORMAT.ID)
      newick = resp['newick']
      # Should fail if bad newick
      tree = dendropy.Tree.get(data=newick, schema='newick')
      # Make sure tree is not None
      self.assertTrue(tree)
      
# .............................................................................
def get_test_suites():
   """
   @summary: Gets the test suites for the module
   """
   loader = unittest.TestLoader()
   testSuites = []
   testSuites.append(loader.loadTestsFromTestCase(TestGetOTTIdsFromGBIFIds))
   testSuites.append(loader.loadTestsFromTestCase(TestInducedSubtree))
   testSuites.append(loader.loadTestsFromTestCase(TestGetTreeFromGBIFIds))
   return testSuites
   
# =============================================================================
# =                                   Main                                    =
# =============================================================================
if __name__ == '__main__':
   
   logging.basicConfig(level=logging.DEBUG)
   
   for suite in get_test_suites():
      unittest.TextTestResult(verbosity=2).run(suite)
   