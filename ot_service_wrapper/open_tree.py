"""
@summary: Module for Open Tree of Life client
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
"""
import json
import urllib2

DEV_SERVER = 'http://141.211.236.35:10999'
INDUCED_SUBTREE_BASE_URL = '{}/induced_subtree'.format(DEV_SERVER)
OTTIDS_FROM_GBIFIDS_URL = '{}/ottids_from_gbifids'.format(DEV_SERVER)

# .............................................................................
class LABEL_FORMAT(object):
   """
   @summary: This class represents label format constants that can be used 
                when calling the induced subtree function
   """
   NAME = 'name'
   ID = 'id'
   NAME_AND_ID = 'name_and_id'

# .............................................................................
def get_ottids_from_gbifids(gbif_ids):
   """
   @summary: Calls the Open Tree 'ottids_from_gbifids' service to retrieve a 
                mapping dictionary from the Open Tree service where each key is 
                one of the provided GBIF identifiers and the value is the 
                corresponding OpenTree id.
   @note: Any GBIF ID that was not found will have a value of None
   @param gbif_ids: A list of GBIF identifiers.  They will be converted to
                       integers in the request.
   """
   # Ids need to be integers
   processed_ids = [int(gid) for gid in gbif_ids]
      
   request_body = {
      "gbif_ids" : processed_ids
   }
   
   headers = {
      'Content-Type' : 'application/json'
   }
   req = urllib2.Request(OTTIDS_FROM_GBIFIDS_URL, 
                         data=json.dumps(request_body), headers=headers)
   
   resp = json.load(urllib2.urlopen(req))
   unmatchedIds = resp['unmatched_gbif_ids']
   
   id_map = resp["gbif_ott_id_map"]
   
   for gid in unmatchedIds:
      id_map[gid] = None
   
   return id_map

# .............................................................................
def induced_subtree(ott_ids, label_format=LABEL_FORMAT.NAME):
   """
   @summary: Calls the Open Tree 'induced_subtree' service to retrieve a tree,
                in Newick format, containing the nodes represented by the 
                provided Open Tree IDs
   @param ott_ids: A list of Open Tree IDs.  These will be converted to into
                      integers
   @param label_format: The label string format to use when creating the tree
                           on the server 
   """
   # Ids need to be integers
   processed_ids = [int(ottid) for ottid in ott_ids]
   request_body = {
      "ott_ids" : processed_ids,
      "label_format" : label_format
   }
   
   headers = {
      'Content-Type' : 'application/json'
   }
   req = urllib2.Request(INDUCED_SUBTREE_BASE_URL, 
                         data=json.dumps(request_body), headers=headers)
   
   resp_str = urllib2.urlopen(req).read()
   return json.loads(resp_str)

