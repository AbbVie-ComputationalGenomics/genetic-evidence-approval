from __future__ import print_function # for python2 compatibility

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 14:02:43 2018

@author: kingea
"""

# This script takes a zipped mesh xml file of the format available from the 
# ftp server and returns and creates a
# text file parsing the output to give all names associated with a UI and
# to create a file with mapped headings for supplementary concepts, if using
# supplementary concepts, and a file with tree number if using MeSH heading.

import xmltodict
import gzip 
import sys

# data_dir = '/Users/kingea/Downloads/'
# mesh_file = 'desc2018.gz'

data_dir = sys.argv[1]
mesh_file = sys.argv[2]

# Determine if mesh descriptors or mesh supplementary concept records

mesh_type = mesh_file[0:4]

if not (mesh_type=='desc' or mesh_type=='supp'):
    print('unrecognized mesh file name', file=sys.stderr)
    sys.exit()

# Year is used to label output automatically
    
mesh_year = mesh_file[4:8]

# check if extension is .gz

mesh_extension = mesh_file[8:11]

if not mesh_extension=='.gz':
    print('expect .gz extension', file=sys.stderr)
    sys.exit() 
    
with gzip.open(data_dir + mesh_file) as fd:
    doc = xmltodict.parse(fd.read())

# Fields accessed differ systematically between descriptor and supplementary concept

if mesh_type=='desc':
    type_prefix1 = 'Descriptor'
    type_prefix2 = 'Descriptor'
else:
    type_prefix1 = 'Supplemental'
    type_prefix2 = 'SupplementalRecord'

record_list = doc[type_prefix1 + 'RecordSet'][type_prefix1 + 'Record']

sr_name_out = open(data_dir + 'MeSH_' + mesh_type + '_' + mesh_year + '.tsv', 'w')
sr_name_out.write('UI' + '\t' + 'Name' + '\t' + 'Preferred' + '\n')

for record in record_list:
    UI = record[type_prefix2 + 'UI']
    Name = record[type_prefix2 + 'Name']['String']
#    sr_name_out.write(UI + '\t' + Name + '\t' + 'TRUE' + '\n')
    if not len(record['ConceptList'])==1:
        print("Unexpected length")
    if not isinstance(record['ConceptList']['Concept'], list):
        concept_list = [record['ConceptList']['Concept']]
    else:
        concept_list = record['ConceptList']['Concept']
    for concept in concept_list:
        concept_name = concept['ConceptName']['String']
        concept_preferred = concept['@PreferredConceptYN']
        if isinstance(concept['TermList']['Term'], list):
            term_list = concept['TermList']['Term']
        else:
            term_list = [concept['TermList']['Term']]
        for term in term_list:
            term_name = term['String']
            term_preferred = term['@RecordPreferredTermYN']
            if term_preferred == 'Y' and concept_preferred=='Y':
                preferred = 'Y'
            else:
                preferred = 'N'
            sr_name_out.write(UI + '\t' + term_name + '\t' + preferred + '\n')

if mesh_type=='supp':
    sr_mapped_out = open(data_dir + 'MeSH_' + mesh_type + 'mapped_' + mesh_year + '.tsv', 'w')
    sr_mapped_out.write('UI' + '\t' + 'MappedUI' + '\n')

    for record in record_list:
        UI = record['SupplementalRecordUI']
        if not len(record['HeadingMappedToList'])==1:
            print("Unexpected length")
        if not isinstance(record['HeadingMappedToList']['HeadingMappedTo'], list):
            mapped_list = [record['HeadingMappedToList']['HeadingMappedTo']]
        else:
            mapped_list = record['HeadingMappedToList']['HeadingMappedTo']
        for heading in mapped_list:
            mapped_UI = heading['DescriptorReferredTo']['DescriptorUI']
            sr_mapped_out.write(UI + '\t' + mapped_UI + '\n')

if mesh_type=='desc':
    sr_mapped_out = open(data_dir + 'MeSH_' + mesh_type + 'tree_' + mesh_year + '.tsv', 'w')
    sr_mapped_out.write('UI' + '\t' + 'TreeNumber' + '\n')
    for record in record_list:
        UI = record['DescriptorUI'] 
        # there are a few special terms that do not have a tree number
        if 'TreeNumberList' in record:
            if not len(record['TreeNumberList'])==1:
                print("Unexpected length")
            if not isinstance(record['TreeNumberList']['TreeNumber'], list):
                tree_list = [record['TreeNumberList']['TreeNumber']]
            else:
                tree_list = record['TreeNumberList']['TreeNumber']
            for tree_num in tree_list:
                sr_mapped_out.write(UI + '\t' + tree_num + '\n')