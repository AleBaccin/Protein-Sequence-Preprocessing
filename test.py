import pandas as pd
import re
from pprint import pprint
import numpy as np


item = "Cytoplasm {ECO:0000250|UniProtKB:Q00534}. Nucleus {ECO:0000250|UniProtKB:Q00534}. Cell projection, ruffle {ECO:0000250|UniProtKB:Q00534}. Cytoplasm, cytoskeleton, microtubule organizing center, centrosome {ECO:0000250|UniProtKB:Q00534}. Note=Localized to the ruffling edge of spreading fibroblasts. Kinase activity only in nucleus (By similarity). Present in the cytosol and in the nucleus in interphase cells and at the centrosome during mitosis from prophase to telophase (By similarity). Localized to the cytosol of neurons and showed prominent staining around either side of the nucleus. {ECO:0000250|UniProtKB:Q00534, ECO:0000269|PubMed:23918663}"#"Cytoplasmic vesicle, secretory vesicle lumen {ECO:0000250|UniProtKB:P26339}. Cytoplasmic vesicle, secretory vesicle membrane {ECO:0000250|UniProtKB:P26339}. Secreted {ECO:0000269|PubMed:11451958}. Note=Associated with the secretory granule membrane through direct interaction to SCG3 that in turn binds to cholesterol-enriched lipid rafts in intragranular conditions. {ECO:0000250|UniProtKB:P26339}.; Serpinin: Secreted {ECO:0000269|PubMed:10781584}. Cytoplasmic vesicle, secretory vesicle {ECO:0000250|UniProtKB:P26339}. Note=Pyroglutaminated serpinin localizes to secretory vesicle. {ECO:0000250|UniProtKB:P26339}."# "Cytoplasm, cytoskeleton. Nucleus {ECO:0000269|PubMed:23558171, ECO:0000269|PubMed:25759381}. Note=Localized in cytoplasmic mRNP granules containing untranslated mRNAs. {ECO:0000250|UniProtKB:P60709}."
dotdataset_encoding_positions = {'M': 0, 'I': 1, 'S': 2, 'F': 3, 'L': 4, 
'N': 5, 'C': 6, 'V': 7, 'D': 8, 'H': 9, 
'G': 10, 'Y': 11, 'K': 12, 'P': 13, 'Q': 14, 
'W': 15, 'A': 16, 'E': 17, 'R': 18, 'T': 19, 
'B': 20, 'J': 20, 'Z': 20, 'X': 21}

item = item.rstrip('.')
item = item.replace('. ', '; ').replace(' {', '{').replace('; ', ';').replace(': ', ';').replace('.', '')
# item = re.split('Note=', item)[0]
splitted = re.split(';', item)
print(splitted)


for string in splitted:
    if 'ECO:0000269' in string:
        string = re.sub('{[^}]+}', '', string)
        s = re.split(', ', string)
        print(s[0])
        break

line = "MPLPLLLAALCLAASPAPARACQLPSEWRPLSEGCRAELAETIVYAKVLALHPEVPGLYNYLPWQYQAGEGGLFYSAEVEMLCDQAWGSMLEVPAGSRLNLTGLGYFSCHSHTVVQDYSYFFFVRMDENYNLLPHGVNFQDAIFPDTQENRRMFSSLFQFANCSQGQQLTTFSSDWEVQEDNRLMCSSVQKALFEEEDHVKKLQQKVATLEKRNRQLRERVKKVKRSLRQARKNSRHLELVNQKLNEKLGASSAQQHINALGREPVRAPYLHG"

print(len(line))
encoding = ""
for c in line:
    for n in range(22):
        if dotdataset_encoding_positions[c] == n:
            encoding += '1 '
        else:
            encoding += '0 '

encoding = encoding.rstrip(' ')

print(encoding)


dict = {}

for c in line:
    if c in dict:
        dict[c] = dict[c] + 1
    else:
        dict[c] = 1

print(dict)
print(len(dict))
line = line.replace(' ', '')
print(len(line))

# f = open('test.txt', 'r')

# print(f.read())

# pprint(splitted)

