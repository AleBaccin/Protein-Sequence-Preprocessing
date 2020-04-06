import sys, getopt
import pandas as pd
import numpy as np
import random
import re
from pprint import pprint
from tqdm import tqdm
from pathlib import Path

title = r"""
      _____                                                                            .__                
    _/ ____\__.__. _____________   ____ _____________  ____   ____  ____   ______ _____|__| ____    ____  
    \   __<   |  | \____ \_  __ \_/ __ \\____ \_  __ \/  _ \_/ ___\/ __ \ /  ___//  ___/  |/    \  / ___\ 
     |  |  \___  | |  |_> >  | \/\  ___/|  |_> >  | \(  <_> )  \__\  ___/ \___ \ \___ \|  |   |  \/ /_/  >
     |__|  / ____| |   __/|__|    \___  >   __/|__|   \____/ \___  >___  >____  >____  >__|___|  /\___  / 
           \/      |__|               \/|__|                     \/    \/     \/     \/        \//_____/  
    """

instructions = """
How to:
prepro.py 
    -e <inputfile> extract all the localization labels with the respective counts and also create a .tab file contained the cleaning dataset
    -f <inputfile> given a cleaned dataset in .tab format, produce the .fasta file
    -a <inputfile> analyze a .fasta file counting the number of entries per class
    -l <inputfile> flatten a fasta, file converging every sequence to one liners, used to simplify encoding operations
    -d <inputfile> generate the .dataset encoding given a .fasta file
    --m <inputfolder> generate the .dataset encodings given a folder containing MSAs files.
"""

dotdataset_encoding_positions = {'M': 0, 'I': 1, 'S': 2, 'F': 3, 'L': 4, 
'N': 5, 'C': 6, 'U': 6, 'V': 7, 'D': 8, 'H': 9, 
'G': 10, 'Y': 11, 'K': 12, 'P': 13, 'Q': 14, 
'W': 15, 'A': 16, 'E': 17, 'R': 18, 'T': 19, 
'B': 20, 'J': 20, 'Z': 20, 'X': 21}
check = ['Cytoplasm', 'Nucleus', 'Secreted', 'Cell membrane', 'Plastid', 'Mitochondrion', 'Endoplasmic reticulum membrane']
fasta_encoding_lables = {'Cytoplasm': 1, 'Nucleus': 2, 'Secreted': 3, 'Cell membrane': 4, 'Plastid': 5, 'Mitochondrion': 6, 'Endoplasmic reticulum membrane': 7, 'Other': 8}
fasta_decoding_lables = {'1': 'Cytoplasm', '2': 'Nucleus', '3': 'Secreted', '4': 'Cell membrane', '5': 'Plastid', '6': 'Mitochondrion', '7':'Endoplasmic reticulum membrane', '8': 'Other'}
sub_loc = 'Subcellular location [CC]'
eco_269 = 'ECO:0000269'

def line_count(inputfile):
    return len(open(inputfile).readlines(  ))

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def clean(inputfile):
    df = pd.read_csv(inputfile, sep='\t')
    df[sub_loc].replace('', np.nan, inplace=True)
    df.dropna(subset=[sub_loc], inplace=True)
    df[sub_loc] = df[sub_loc].map(lambda x: re.split('Note=', x)[0])
    df = df[df[sub_loc].str.contains(eco_269) == True]
    df[sub_loc] = df[sub_loc].map(lambda x: x.replace('SUBCELLULAR LOCATION: ', ''))
    df = df.loc[(df['Length'] > 30) & (df['Length'] < 10000)]
    return df

def label_extraction(df):
    c = 0
    dict = {}

    for c, row in tqdm(df.iterrows(), total=df.shape[0]):
        item = row[sub_loc]
        item = item.rstrip('.')
        item = item.replace('. ', '; ').replace(' {', '{').replace('; ', ';').replace(': ', ';').replace('.', '')
        splitted = re.split(';', item)
        for string in splitted:
            if eco_269 in string:
                string = re.sub('{[^}]+}', '', string)
                s = re.split(', ', string)
            if s[0] in check:
                df.at[c, sub_loc] = s[0]
            else:
                df.at[c, sub_loc] = 'Other'
            if s[0] in dict:
                dict[s[0]] = dict[s[0]] + 1
            else:
                dict[s[0]] = 1
            break
    c += 1
    listoflabels = [(k, v) for k, v in sorted(dict.items(), key=lambda x: x[1], reverse=True)]
    return listoflabels, df

def save_outputs_of_extraction(listoflabels, df):
    df.to_csv('datasets/cooked.tab', columns = ['Entry', 'Sequence', 'Length', 'Subcellular location [CC]'], sep ='\t')
    with open('datasets/labels_extracted.txt', 'w') as f:
        for t in listoflabels:
            f.write(' '.join(str(s) for s in t) + '\n')

def fasta_encoding(inputfile):
    df = pd.read_csv(inputfile, sep='\t')
    f= open("datasets/entries.fasta","w+")

    c = 0
    for Entry, Sequence, Sub_loc in tqdm(zip(df.Entry, df.Sequence, df['Subcellular location [CC]']), total=df.shape[0]):
        f.write(">"+ str(c) + "|" + Entry + "|" + Sub_loc.replace(' ', '_') + "(" + str(fasta_encoding_lables[Sub_loc]) + ")\n")
        splitted = list(chunkstring(Sequence, 64))
        for x in splitted:
            f.write(x + "\n")
        c += 1

def analyze_fasta(inputfile):
    dict = {}

    l = line_count(inputfile)

    with open(inputfile) as f:
        protein_count = 0
        for line in tqdm(f, total=l):
            if line[0] == '>':
                label = re.search(r"\((\w+)\)", line).group(1)
                label = fasta_decoding_lables[label]
                if label in dict:
                    dict[label] = dict[label] + 1
                else:
                    dict[label] = 1
                protein_count += 1
        print(protein_count)
    listoflabels = [(k, v) for k, v in sorted(dict.items(), key=lambda x: x[1], reverse=True)]
    print(listoflabels)

def flatten_fasta_to_file(inputfile):
    ff= open("datasets/flattened_entries.fasta","w+")

    l = line_count(inputfile)

    with open(inputfile) as f:
        flattened_sequence = ''
        for line in tqdm(f, total=l):
            if line[0] == '>':
                if flattened_sequence != '':
                    flattened_sequence += '\n'

                ff.write(flattened_sequence + line)
                flattened_sequence = ''
            else:
                flattened_sequence += line.replace('\n', '')
        ff.write(flattened_sequence)

def add_msas(inputfolder):
    files = ['train.dataset', 'test.dataset', 'validation.dataset']

    for fs in files:
        print(f'Clipping MSA to {fs}')
        with open(f'datasets/NoMSAdataset/{fs}', 'r') as f_dataset_nomsa:
            f_dataset_msa = open(f'datasets/MSAdataset/{fs}','w+')
            l = line_count(f'datasets/NoMSAdataset/{fs}')
            count = 1
            sequence = ''
            f_dataset_msa.write(f'{int(l / 4)}\n22 8\n')
            for line in tqdm(f_dataset_nomsa, total=l):
                if count % 4 == 1:
                    sequence = line.rstrip('\n')

                if count % 4 == 3:
                    msa_sequence = Path(f'{inputfolder}\{sequence}.dataset')
                    if msa_sequence.is_file():
                        msa_sequence = open(f'{inputfolder}\{sequence}.dataset')
                        lines = msa_sequence.readlines()
                        f_dataset_msa.write(lines[4])
                        msa_sequence.close()
                    else:
                        f_dataset_msa.write('placeholder\n')
                else:
                    f_dataset_msa.write(line)
                count += 1
            f_dataset_msa.close()
        f_dataset_nomsa.close()

def encode_flattened_fasta(inputfile):
    train = open("datasets/NoMSAdataset/train.dataset","w+")
    test = open("datasets/NoMSAdataset/test.dataset","w+")
    val = open("datasets/NoMSAdataset/validation.dataset","w+")

    l = line_count(inputfile)

    with open(inputfile) as f:
        label = 0
        rdn = 0
        for line in tqdm(f, total=l):
            if rdn < 3 and rdn > -1:
                ff = train
            elif rdn < 4 and rdn > 2:
                ff = test
            else:
                ff = val

            if line[0] == '>':
                line = line.split('|')
                label = re.search(r"\((\w+)\)", line[2]).group(1)
                ff.write(line[1] + '\n')
            else:
                to_write = ''
                line = line.rstrip('\n')
                to_write += str(len(line)) + '\n'
                encoding = ""
                
                for c in line:
                    for n in range(22):
                        if dotdataset_encoding_positions[c] == n:
                            encoding += '1 '
                        else:
                            encoding += '0 '

                encoding = encoding.rstrip(' ')
                to_write += encoding + '\n'
                to_write += label + '\n'
                ff.write(to_write)
                rdn = random.randrange(5)

def main(argv):
    print(title)
    inputfile = ''
    try:
        opts, _ = getopt.getopt(argv,"he:f:c:a:l:d:m:",["ifile=", "ffile=", "afile=", "lfile=", "dfile=", "mfile=" ""])
    except getopt.GetoptError:
        print(instructions)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(instructions)
            sys.exit()
        elif opt in ("-e", "--ifile"):
            print("+----Extracting-Labels----+")
            inputfile = arg
            labels, df = label_extraction(clean(inputfile))
            save_outputs_of_extraction(labels, df)
        elif opt in ("-f", "--ffile="):
            print("+----Encoding-to-Fasta----+")
            inputfile = arg
            fasta_encoding(inputfile)
        elif opt in ("-a", "--afile="):
            print("+-----Analyzing-Fasta-----+")
            inputfile = arg
            analyze_fasta(inputfile)
        elif opt in ("-l", "--lfile="):
            print("+-----Flattening-Fasta----+")
            inputfile = arg
            flatten_fasta_to_file(inputfile)
        elif opt in ("-d", "--dfile="):
            print("+--dotdataset-from-Fasta--+")
            inputfile = arg
            encode_flattened_fasta(inputfile)
        elif opt in ("-m", "--mfil="):
            print("+--Clipping MSA to dotdataset--+")
            inputfolder = arg
            add_msas(inputfolder)


if __name__ == "__main__":
   main(sys.argv[1:])