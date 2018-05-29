"""
@summary: Module for Open Tree of Life client
@author: CJ Grady
@version: 1.0.0
@status: alpha

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

INDUCED_SUBTREE_BASE_URL = 'http://141.211.236.35:10999/induced_subtree'

# .............................................................................
class LABEL_FORMAT(object):
   NAME = 'name'
   ID = 'id'
   NAME_AND_ID = 'name_and_id'

# .............................................................................
def induced_subtree(ott_ids, label_format=LABEL_FORMAT.NAME, return_flo=False):
   """
   """
   request_body = {
      "ott_ids" : ott_ids,
      "label_format" : label_format
   }
   
   headers = {
      'Content-Type' : 'application/json'
   }
   req = urllib2.Request(INDUCED_SUBTREE_BASE_URL, 
                         data=json.dumps(request_body), headers=headers)
   
   if return_flo:
      return urllib2.urlopen(req)
   else:
      return urllib2.urlopen(req).read()

