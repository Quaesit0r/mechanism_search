# Constructing A Large-Scale Reaction Mechanism Dataset Using A Template-Based Search Algorithm
This repository contains the code and resources for Constructing *A Large-Scale Reaction Mechanism Dataset Using A Template-Based Search Algorithm*. The project creates a reaction mechanism dataset using reaction SMILES and faciliate reaction mechanism predictions using different models.

# Requirements
For the dataset creation, the following packages are required:

python==3.8.19
numpy==1.24.3
pandas==2.0.3
rdkit==2024.3.3
scipy==1.10.1
tqdm==4.65.0

For model training and prediction, the packages required is recorded in the respective folders.

# Dataset creation
For dataset creation, please prepare your database in a file which each line contains a reaction SMILES. Then go to main.py and change the path to your database file and run the main.py file to create the dataset.

# Model training and prediction
We have examined our dataset with the following models:

    Retrosynthesis-Prediction Model: Located in the Retrosynthesis-Prediction folder.
    Transformer: Located in the Transformer folder.
    Graph2SMILES: Located in the Graph2SMILES folder.

For detailed setup and usage instructions, refer to the README.md file located in each model folder. Simply navigate to your preferred model’s folder and follow the steps outlined in its README.

# License
This project is licensed under the MIT License.