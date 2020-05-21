Alessandro Baccin, 16724489, alessandro.baccin@ucdconnect.ie

# Protein-Sequence-Preprocessing
This repository contains scripts that I needed in order carry out the preprocessing of my datasets for my 4th year FYP project.

Note: metrics.ipynb and weights_count.py are used to: 
    - calculate the metrics for each .predictions files.
    - count the number of weights in the generated .model files.

⚠️ IMPORTANT: Due to the size of the produced data it was not possible to have it uploaded in GitHub. Here are the Google Drive links:
    - Preprocessing: https://drive.google.com/open?id=1bfCnXEeWZWlbAKF3fFx0tkE62LvKo0Nu
    - Trained models: https://drive.google.com/open?id=1fXEQwor4Jv8JBR8LkCB5GI4QAPyxU0dA
    
Some of the attempted trained models are residing on the tanalla.ucd.ie machine.

## How does prepro.py work?

This script provides 6 different modes of execution: 

```python
"""
    -l | --lfile    <inputfile>     Extract, save and visualize all the present labels.
                                    <inputfile> needs to be located in the root directory.
    -o | --model    <modelname>     Specify a model contained in settings.json
    -e | --efile    <inputfile>     Create a .tab file containing sequences of classes specified in the required model. 
                                    <inputfile> needs to be located in the root directory.
    -f | --ffile    <inputfile>     Given a cleaned dataset in .tab format, produce the .fasta file.
                                    <inputfile> needs to be located in the model's directory 
    -a | --ffile    <inputfile>     Analyze a .fasta file counting the number of entries per class.
                                    <inputfile> needs to be located in the model's directory
    -d | --afile    <inputfile>     Generate the three .dataset encoding files given a .fasta file
                                    Present MSAs are also taken into consideration, and sequences not preseting the MSA information are binned.
                                    Ouputs are stored into the generated "NoMSAdataset" folder.
                                    <inputfile> needs to be located in the model's directory 
    -m | --mfile                    Given the folder containing MSAs files specified in the settings.json, attach the MSA information to the .dataset files (new files are created).
                                    Ouputs are stored into the generated "MSAdataset" folder.
                                    <inputfile> needs to be located in the model's directory 
"""
```

In order for the script to work, a settings.json have to be specified, the -o flag can be omitted (but all operations will be carried out using "model_one" as configuration), but the settings.json needs to be present. The initial .tab is also to be placed in the root directory. To generate the MSAs it is also necessary to place the alignments into a folder inside the root directory.

## Setup

Example of the settings.json file.

```json
{
    "msas_dir" : "MSAs",
    "model_one" : {
        "outputdir" : "datasets",
        "labels": {
            "Cytoplasm": 1,
            "Nucleus": 2,
            "Secreted": 3,
            "Cell membrane": 4,
            "Endoplasmic reticulum membrane": 5,
            "Plastid": 6,
            "Mitochondrion": 7,
            "Other": 8
        }
    }
```

Example of the required root folder:

.

├───MSAs

│   ├───A0A023PZB3.dataset

│   └─── . . .

├───uniprot.tab

├───prepro.py

└───settings.json
