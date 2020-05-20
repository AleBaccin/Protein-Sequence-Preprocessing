import pandas as pd
import re
from pprint import pprint
import numpy as np
import os
from tqdm import tqdm

names = ['7-class-weights-1l', '7-class-weights-2l', '7-class-weights-2l_20_15', '7-class-weights-3l', '7-class_weights_3l_20', '7-class-weights-4l', '7-class-weights-3sl', '6-class-weights-2l', '5-class-weights-2l', '5-class-weights-3l',  '8-class-weights-2l']
results = open(os.path.join('weights', 'results.txt'),'w+')

for name in tqdm(names, total= len(names)):
    weights = open(os.path.join('weights', f'{name}.model'),'r')
    flat = open(os.path.join('weights', f'{name}_flat.txt'),'w+')

    for line in weights:
        line = line.replace(' ', '\n')
        flat.write(line)
    weights.close()
    flat.close()

    flat = open(os.path.join('weights', f'{name}_flat.txt'),'r')
    results.write(f'{name} # of weights: {len(flat.readlines())}\n')
    flat.close()