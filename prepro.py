import sys, getopt
import pandas as pd
import numpy as np
import random
import re
import json
import os
from sklearn.model_selection import KFold
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
    -l | --lfile    <inputfile>     Extract, save and visualize all the present labels.
    -o | --model    <modelname>     Specify a model contained in settings.json
    -e | --efile    <inputfile>     Create a .tab file containing sequences of classes specified in the required model. 
    -f | --ffile    <inputfile>     Given a cleaned dataset in .tab format, produce the .fasta file.
    -a | --ffile    <inputfile>     Analyze a .fasta file counting the number of entries per class.
    -d | --afile    <inputfile>     Generate the three .dataset encoding files given a .fasta file
                                    Present MSAs are also taken into consideration, and sequences not preseting the MSA information are binned.
                                    Ouputs are stored into the generated "NoMSAdataset" folder.
    -m | --mfile                    Given the folder containing MSAs files specified in the settings.json, attach the MSA information to the .dataset files (new files are created).
                                    Ouputs are stored into the generated "MSAdataset" folder.
"""

dotdataset_encoding_positions = {'M': 0, 'I': 1, 'S': 2, 'F': 3, 'L': 4, 
'N': 5, 'C': 6, 'U': 6, 'V': 7, 'D': 8, 'H': 9, 
'G': 10, 'Y': 11, 'K': 12, 'P': 13, 'Q': 14,
'W': 15, 'A': 16, 'E': 17, 'R': 18, 'T': 19, 
'B': 20, 'J': 20, 'Z': 20, 'X': 21}
sub_loc = 'Subcellular location [CC]'
eco_269 = 'ECO:0000269'
n_letter_encoding_digits = 22
random.seed(41)

with open('settings.json', 'r') as json_file:
  data = json.load(json_file)

msas= data['msas_dir']
output_dir = data['model_one']['outputdir']
fasta_encoding_lables = data['model_one']['labels']
check = list(fasta_encoding_lables.keys())[:-1]

fasta_decoding_lables = {}
for k, v in fasta_encoding_lables.items():
    if str(v) in fasta_decoding_lables:
        fasta_decoding_lables[str(v)] += f' / {k}'
    else:
        fasta_decoding_lables[str(v)] = k

def line_count(inputfile):
    return len(open(f'{output_dir}/{inputfile}').readlines())

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def clean(inputfile):
    df = pd.read_csv(f'{output_dir}/{inputfile}', sep='\t')
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
            else:
                continue

            if s[0] in check:
                df.at[c, sub_loc] = s[0]
            else:
                if 'Other' in fasta_encoding_lables.keys():
                    df.at[c, sub_loc] = 'Other'
                else:
                    df.drop(c, inplace = True)
                
            if s[0] in dict:
                dict[s[0]] = dict[s[0]] + 1
            else:
                dict[s[0]] = 1
            break
    c += 1
    listoflabels = [(k, v) for k, v in sorted(dict.items(), key=lambda x: x[1], reverse=True)]

    return listoflabels, df

def fasta_encoding(inputfile):
    filename = inputfile.split('.')[0]
    df = pd.read_csv(f'{output_dir}/{inputfile}', sep='\t')
    f= open(f'{output_dir}/entries_{filename}.fasta','w+')

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

    with open(f'{output_dir}/{inputfile}') as f:
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
    print("+-----Flattening-Fasta----+")
    l = line_count(inputfile)

    ff = ""

    with open(f'{output_dir}/{inputfile}') as f:
        flattened_sequence = ''
        for line in tqdm(f, total=l):
            if line[0] == '>':
                if flattened_sequence != '':
                    flattened_sequence += '\n'

                ff += flattened_sequence + line
                flattened_sequence = ''
            else:
                flattened_sequence += line.replace('\n', '')
        ff += flattened_sequence

    return ff

def add_msas():
    files = ['train.dataset', 'test.dataset', 'validation.dataset']

    for fs in files:
        print(f'Clipping MSA to {fs}')
        with open(f'{output_dir}/NoMSAdataset/{fs}', 'r') as f_dataset_nomsa:
            f_dataset_msa = open(f'{output_dir}/MSAdataset/{fs}','w+')
            l = line_count(f'NoMSAdataset/{fs}')
            count = 1
            sequence = ''
            missed = 0

            f_dataset_msa.write(f'{int((l - 2) * 0.25)}\n{n_letter_encoding_digits} {len(fasta_decoding_lables)}\n')
            for line in tqdm(f_dataset_nomsa.readlines()[2:], total=l):
                if count % 4 == 1:
                    sequence = line.rstrip('\n')
                if count % 4 == 3:
                    msa_sequence = Path(f'{msas}\\{sequence}.dataset')
                    if msa_sequence.is_file():
                        msa_sequence = open(f'{msas}\\{sequence}.dataset')
                        lines = msa_sequence.readlines()
                        f_dataset_msa.write(lines[4])
                        msa_sequence.close()
                    else:
                        missed += 1
                else:
                    msa_sequence = Path(f'{msas}\\{sequence}.dataset')
                    if msa_sequence.is_file():
                        f_dataset_msa.write(line)
                count += 1

            f_dataset_msa.seek(0, 0)
            f_dataset_msa.write(f'{int((l - 2)  * 0.25 - missed)}\n')
            f_dataset_msa.close()
        f_dataset_nomsa.close()

def encode_flattened_fasta(flattened_fasta):
    train = open(f'{output_dir}/NoMSAdataset/train.dataset','w+')
    test = open(f'{output_dir}/NoMSAdataset/test.dataset','w+')
    val = open(f'{output_dir}/NoMSAdataset/validation.dataset','w+')
    sets = [1, 2, 3]
    weights = [0.80, 0.10, 0.10]

    label = 0
    rdn = random.choices(sets, weights)[0]

    train_count =  0
    val_count =  0
    test_count =  0

    rough_train_size = int(len(flattened_fasta.splitlines()) * 0.5 * 0.8)
    rough_test_size = int(len(flattened_fasta.splitlines()) * 0.5 * 0.1)
    rough_val_size = int(len(flattened_fasta.splitlines()) * 0.5 * 0.1)

    train.write(f'{rough_train_size}\n{n_letter_encoding_digits} {len(fasta_decoding_lables)}\n')
    test.write(f'{rough_test_size}\n{n_letter_encoding_digits} {len(fasta_decoding_lables)}\n')
    val.write(f'{rough_val_size}\n{n_letter_encoding_digits} {len(fasta_decoding_lables)}\n')
    msa_sequence = Path(msas)
    for line in tqdm(flattened_fasta.splitlines(), total=len(flattened_fasta.splitlines())):
        if rdn == 1:
            ff = train
        elif rdn == 2:
            ff = val
        elif rdn == 3:
            ff = test

        if line[0] == '>':
            line = line.split('|')
            label = re.search(r"\((\w+)\)", line[2]).group(1)
            msa_sequence = Path(msas)
            msa_sequence = Path(msa_sequence, f'{line[1]}.dataset')
            if msa_sequence.is_file():
                ff.write(line[1] + '\n')
        else:
            to_write = ''
            line = line.rstrip('\n')
            to_write += str(len(line)) + '\n'
            encoding = ""
            
            if msa_sequence.is_file():
                for c in line:
                    for n in range(n_letter_encoding_digits):
                        if dotdataset_encoding_positions[c] == n:
                            encoding += '1 '
                        else:
                            encoding += '0 '

                encoding = encoding.rstrip(' ')
                to_write += encoding + '\n'
                to_write += label + '\n'
                ff.write(to_write)

                if rdn == 1:
                    train_count += 1
                elif rdn == 2:
                    val_count += 1
                elif rdn == 3:
                    test_count += 1

            rdn = random.choices(sets, weights)[0]

    train.seek(0, 0)
    train.write(f'{int(train_count)}\n')
    train.close()

    test.seek(0, 0)
    test.write(f'{int(test_count)}\n')
    test.close()

    val.seek(0, 0)
    val.write(f'{int(val_count)}\n')
    val.close()

argv = sys.argv[1:]
print(title)
inputfile = ''
try:
    opts, _ = getopt.getopt(argv,"ho:l:e:f:a:d:m",["model=", "lfile=","efile=", "ffile=", "afile=", "dfile=", "mfile="])
except getopt.GetoptError:
    print(instructions)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(instructions)
        sys.exit()
    elif opt in ("-o", "--model"):
        model = arg
        output_dir = data[model]['outputdir']
        if not os.path.exists(os.path.join(output_dir, 'MSAdataset')):
            os.mkdir(os.path.join(output_dir, 'MSAdataset'))
        if not os.path.exists(os.path.join(output_dir, 'NoMSAdataset')):
            os.mkdir(os.path.join(output_dir, 'NoMSAdataset'))

        fasta_encoding_lables = data[model]['labels']

        if 'Other' in fasta_encoding_lables.keys():
            check = list(fasta_encoding_lables.keys())[:-1]
        else:
            check = list(fasta_encoding_lables.keys())

        fasta_decoding_lables = {}
        for k, v in fasta_encoding_lables.items():
            if str(v) in fasta_decoding_lables:
                fasta_decoding_lables[str(v)] += f' / {k}'
            else:
                fasta_decoding_lables[str(v)] = k

        print(f'+---Working on {output_dir}---+')
    elif opt in ("-l", "--lfile"):
        print("+----Extracting-Labels----+")
        inputfile = arg
        labels, _ = label_extraction(clean(inputfile))
        print(labels)
        with open(f'{output_dir}/labels_extracted.txt', 'w') as f:
            for t in labels:
                f.write(' '.join(str(s) for s in t) + '\n')
    elif opt in ("-e", "--efile"):
        print("+----Extracting-reduced-.tab----+")
        inputfile = arg
        _, df = label_extraction(clean(inputfile))

        filename = inputfile.split('.')[0]
        df.to_csv(f'{output_dir}/cooked_{filename}.tab', columns = ['Entry', 'Sequence', 'Length', 'Subcellular location [CC]'], sep ='\t')
    elif opt in ("-f", "--ffile="):
        print("+----Encoding-to-Fasta----+")
        inputfile = arg
        fasta_encoding(inputfile)
    elif opt in ("-a", "--afile="):
        print("+-----Analyzing-Fasta-----+")
        inputfile = arg
        analyze_fasta(inputfile)
    elif opt in ("-d", "--dfile="):
        print("+--dotdataset-from-Fasta--+")
        inputfile = arg
        encode_flattened_fasta(flatten_fasta_to_file(inputfile))
    elif opt in ("-m", "--mfile="):
        print("+--Clipping MSA to dotdataset--+")
        inputfolder = arg
        add_msas()