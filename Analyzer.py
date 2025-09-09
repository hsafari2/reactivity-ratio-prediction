
import pandas as pd
from  rdkit import Chem
import matplotlib.pyplot as plt
import numpy as np

#defince static class Analyzer
class StochasticAnalyzer:

    #definitaion of a function that analyze monomer distribution analyze of the system
    @staticmethod
    
    def analyzeFrequency(data):
        # Create a new DataFrame for SMILES conversion
        analysis_df = pd.DataFrame()
    
        # Convert SMILES to Canonical SMILES and store in the new DataFrame
        data['Canonical_SMILES_A'] = data['SMILES_A'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)) if pd.notnull(x) else None)
        data['Canonical_SMILES_B'] = data['SMILES_B'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)) if pd.notnull(x) else None)
    
        # Preparing a dictionary to hold reactivity ratios for each molecule
        reactivity_ratios = {}
    
        # Loop through each row in the data to populate the dictionary
        for index, row in data.iterrows():
            # For Canonical_SMILES_A
            smiles_a = row['Canonical_SMILES_A']
            if smiles_a:
                if smiles_a not in reactivity_ratios:
                    reactivity_ratios[smiles_a] = []
                reactivity_ratios[smiles_a].append(row['r1'])
    
            # For Canonical_SMILES_B
            smiles_b = row['Canonical_SMILES_B']
            if smiles_b:
                if smiles_b not in reactivity_ratios:
                    reactivity_ratios[smiles_b] = []
                reactivity_ratios[smiles_b].append(row['r2'])
    
        # Creating the final DataFrame
        smile = []
        frequencies = []
        reactivities = []
        for smiles, ratios in reactivity_ratios.items():
            smile.append(smiles)
            frequencies.append(len(ratios))
            reactivities.append(ratios)
    
        freq_dist_df = pd.DataFrame({'SMILE': smile, 'Frequency': frequencies, 'Reactivity Ratios': reactivities})
    
        return freq_dist_df


class Normalizer:
    @staticmethod
    def Normalize(data: pd.DataFrame, column_name: str, method: str = 'min-max') -> pd.DataFrame:
        # Check if the column exists in the DataFrame
        if column_name not in data.columns:
            raise ValueError(f"Column '{column_name}' not found in DataFrame")

        # Perform Min-Max normalization
        max_value = data[column_name].max()
        min_value = data[column_name].min()
        data[column_name] = (data[column_name] - min_value) / (max_value - min_value)
      
        return data



#plot data
class Plot:
    def __init__(self, data:pd.DataFrame = None):
        self.data = data


    
    def Plot2DFrequency (self):
        fig, ax = plt.subplots()

        frequency = [row['Frequency'] for index, row in self.data.iterrows()]
        reactivity_value = [np.mean(row['Reactivity Ratios']) for index, row in self.data.iterrows()]

        ax.scatter(frequency, reactivity_value, alpha=0.5, s=10)  # Adjust alpha for transparency, s for size

        ax.set_xlabel('Frequency')
        ax.set_ylabel('Average Reactivity Ratio')
        ax.set_title('Frequency vs Average Reactivity Ratio Scatter Plot')
        ax.grid(True)  # Add gridlines

        plt.show()









    

 
