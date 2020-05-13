import pandas as pd
from tqdm import tqdm
from sklearn import metrics
import os

names = ['all_msa_1l', 'all_msa_2l', 'all_msa_3l', 'all_msa_3l_20', 'all_msa_4l', 'all_nomsa_1l', 'all_nomsa_2l', 'all_nomsa_3l', 'all+other_msa_2l', 'all_ems', 'all_hum_fun_MSA']
results = open(os.path.join('predictions', 'results.txt'),'w+')

for name in tqdm(names, total= len(names)):
    preds = open(os.path.join('predictions', f'{name}.predictions'),'r')

    df = pd.DataFrame(columns=['name', 'actual', 'predicted'])

    k = 0
    c = 1
    for line in preds.readlines():
        line=line.rstrip('\n')
        if c % 6 == 2:
            df.at[k, 'name'] = line
        if c % 6 == 4:
            df.at[k, 'actual'] = f"{int(line) + 1}"
        if c % 6 == 5:
            df.at[k, 'predicted'] = f"{int(line) + 1}"
            k += 1
        c += 1

    y_true = df['actual']
    y_pred = df['predicted']

    results.write(f'+-----------------{name}-----------------+\n')
    results.write("{} \n Weighted accuracy: {:2.3f} \n".format(str(metrics.classification_report(y_true, y_pred, digits=2)), metrics.balanced_accuracy_score(y_true, y_pred)))
    results.write('+-----------------END-----------------+\n')
    preds.close()
