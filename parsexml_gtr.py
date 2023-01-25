"""
Script that parses an XML file and produces a TSV file with several columns, each one representing one field of the XML file.
"""


import os
import csv
import re
import time
import logging
import requests
from lxml import etree as ET


''' REAL XML STRUCTURE
GTRPublicData (root)
    GTRLabData (just a container, no info)
        GTRLab (info about the lab)
        GTRLabTest (info about the trial)
    .
    .
    .
'''

start_time = time.time()

destination_file = '/home/ricard/Documents/xml_parser/gtr_processed_data_v5.tsv' # lab
output_fh = open(os.path.join(destination_file), "w")
output_fh.write("LabName\tGTRAccession\tID\tCUI\tTraitSetCUIs\tMeasureSetGenes\tDiseaseName\tGeneSymbol\tGeneID\tdbSNPs\tVariants\tProtein\tEnzyme\tPurpose\tTestType\tPMIDs\tCategoryValue\tTestName\n")
print("Output file: "  + destination_file)

input_path = '/home/ricard/Documents/xml_parser' # lab
file_list = os.listdir(input_path)
input_file = os.path.join(input_path, 'gtr_ftp.xml') # cannot run the whole dataset from home cause it freezes
print("Input file: " + input_file)

try:
    tree = ET.parse(input_file)
    print(">>>> Input file loaded.")
    print(">>>> Parsing xml file ...")
except ET.XMLSyntaxError as e:
    print(">>>> There is some problem with the input file given !!!")
    print(e)

root = tree.getroot()
print("Number of elements in GTRLabData: " + str(len(root[0])))

keys = ['LabName','GTRAccession','ID','CUI','TraitSetCUIs','MeasureSetGenes', 'DiseaseName','GeneSymbol','GeneID','dbSNPs','Variants','Protein','Enzyme','Purpose','TestType','PMIDs','CatValue','TestName']

# Define writer to write output file
tsv_writer = csv.DictWriter(output_fh, keys, delimiter="\t", quoting=csv.QUOTE_MINIMAL)

lab_name = False
i = 0

# Parse the xml file storing only the desired fields
for labdata in root: # root = 'GTRPublicData'
    for elem in labdata: # labdata = 'GTRLabData'
        if elem.tag == 'GTRLab':
            lab_name = elem.find('Organization/Name').text

        elif elem.tag == 'GTRLabTest':
            id = elem.find('ClinVarSet').get('ID')

            gtr_accession = elem.get('GTRAccession')

            cui = elem.find('ClinVarSet/DescrSet').get('CUI')

            disease_name = elem.find('ClinVarSet/DescrSet').text if elem.find('ClinVarSet/DescrSet').get('Type') == 'preferred' else elem.find('ClinVarSet/DescrSet').text

            genesSymb_list = [subelem.find('Symbol').text for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure') if subelem.get('Type') == 'Gene' and subelem.find('Symbol') is not None]
            gene_symbol = ", ".join(genesSymb_list)
            
            geneID_list = [subelem.find('XRef').get('ID') for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure') if subelem.get('Type') == 'Gene' and subelem.find('XRef') is not None]
            geneID = ", ".join(geneID_list)

            # rsIDs
            for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure'):
                dbSNPs = ['rs' + name.get('ID') for name in subelem.iter('XRef') if name.get('Type') == 'rs' and name.get('DB') == 'dbSNP']
            dbSNPs = ", ".join(dbSNPs) if len(dbSNPs) > 1 else dbSNPs[0] if len(dbSNPs) == 1 else ''

            # variant name(s)
            variants_list = [subelem.find('Name').text for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure') if subelem.get('Type') == 'Variant' and subelem.find('Name') is not None]
            variants = ", ".join(variants_list)

            protein_list = [subelem.find('Name').text for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure') if subelem.get('Type') == 'Protein' and subelem.find('Name') is not None]
            protein = ", ".join(protein_list)

            enzyme_list = [subelem.find('Name').text for subelem in elem.findall('ClinVarSet/ClinVarAssertion/MeasureSet/Measure') if subelem.get('Type') == 'Enzyme' and subelem.find('Name') is not None]
            enzyme = ", ".join(enzyme_list)

            # purpose field
            if elem.find('Indications') is not None:
                purpose_list = [subelem.text for subelem in elem.find('Indications').iter('Purpose') if not subelem.text.isspace()]
                purpose = ", ".join(purpose_list) if len(purpose_list) > 1 else purpose_list[0] if len(purpose_list) == 1 else ''
            else:
                purpose = ''

            # testtype field
            if elem.find('Indications') is not None:
                testype_list = [subelem.text for subelem in elem.find('Indications').iter('TestType') if not subelem.text.isspace()]
                testype = ", ".join(testype_list) if len(testype_list) > 1 else testype_list[0] if len(testype_list) == 1 else ''
            else:
                testype = ''

            # if elem.find('MethodAdd/Protocol'):
            #     pmids = [subelem.text.strip("'") for subelem in elem.find('MethodAdd/Protocol').iter('PMID') if not subelem.text.isspace() and subelem.text is not None] # list
            #     pmids = ", ".join(pmids)
            #     # print('PMIDs:', pmids)
            # else:
            #     pmids = ''
            
            pmids_list = [pmid.text for pmid in elem.findall(".//PMID") if elem.findall(".//PMID")]
            pmids = ", ".join(pmids_list) if pmids_list is not None else ''

            # category field, value attribute
            if elem.find('Method/TopCategory/Category') is not None:
                catValue_list = [subelem.get('Value') for subelem in elem.findall('Method/TopCategory/Category')] # str
                cat_value = ", ".join(catValue_list) if len(catValue_list) > 1 else catValue_list[0] if len(catValue_list) == 1 else ''
            else:
                cat_value = ''

            if elem.find('TestName') is not None:
                test_name = elem.find('TestName').text
                if "\n" in test_name: # in some cases there was a bad string formatting in the XML input file therefore we have to check that
                    test_name = elem.find('TestName').text.replace("\n"," ").replace("        ","")
            else:
                test_name = ''
                
            # Here we want to obtain all the 'XRef' elements that are child of the same 'ClinVarAssertion', therefore we need another for loop because the 'ClinVarAssertion' element is not a direct child of 'GTRLabTest', in fact is its grandchild
            traitSet_dict = {}; measureSet_dict = {}; j = 0
            for cvAssertion in elem.find('ClinVarSet').iter('ClinVarAssertion'):
                j += 1
                trait_set_cuis = [subelem.get('ID') for subelem in cvAssertion.findall('TraitSet/Trait/XRef') if subelem.get('Type') == 'CUI']
                trait_set_cuis = ", ".join(trait_set_cuis) if trait_set_cuis is not None else ''
                traitSet_dict[j] = trait_set_cuis
                
                measure_set_geneIDs = [subelem.get('ID') for subelem in cvAssertion.findall('MeasureSet/Measure/XRef') if subelem.get('DB') == 'Gene']
                measure_set_geneIDs = ", ".join(measure_set_geneIDs) if measure_set_geneIDs is not None else ''
                measureSet_dict[j] = measure_set_geneIDs


            i += 1
            if i % 10000 == 0:
                print(">>>> Processed {} elements.".format(i))

            if lab_name:
                for trait_set_cuis, measure_set_geneIDs in zip(traitSet_dict.values(), measureSet_dict.values()): # we want to write one row for each 'ClinVarAssertion' element that was found in that 'GTRLabTest', therefore we need a for loop
                    row_dict = {
                        'LabName': lab_name,
                        'GTRAccession': gtr_accession,
                        'ID': id,
                        'CUI': cui,
                        'TraitSetCUIs': trait_set_cuis,
                        'MeasureSetGenes': measure_set_geneIDs,
                        'DiseaseName': disease_name,
                        'GeneSymbol': gene_symbol,
                        'GeneID': geneID,
                        'dbSNPs': dbSNPs,
                        'Variants': variants,
                        'Protein': protein,
                        'Enzyme': enzyme,
                        'Purpose': purpose,
                        'TestType': testype,
                        'PMIDs': pmids,
                        'CatValue': cat_value,
                        'TestName': test_name
                    }
                    tsv_writer.writerow(row_dict)
                gtr_accession, id, cui, trait_set_cuis, measure_set_geneIDs, disease_name, gene_symbol, geneID, dbSNPs, variants, protein, enzyme, purpose, testype, pmids, cat_value, test_name = '','','','','','','','','','','','','','','','',''

output_fh.close()

print(">>>> DONE !!!")
print(">>>> RUN TIME --- {} SECONDS ---".format(time.time()-start_time))

# # get all tags ordered in a list
# alist=[]
# for elem in tree.iter():
#     alist.append(elem.tag)
# alist = list(set(alist))
# print(sorted(alist))

