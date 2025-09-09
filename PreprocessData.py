# in this part we define some class for preprocessign data like reading and doing some preocess to be ready for using

# add dependencies
import pandas as pd
import os
from rdkit import *
import requests



#define Class Reader
class Reader:
    
    #Initialize the reader
    def __init__(self, name :str = None, address=None , data=None):
        self.name = name
        self.address = address
        self.data = None  # Attribute to hold the DataFrame

    #set name for the Reader
    def setName (self , name):
        self.name = name
    
    #set Address for the Reader
    def setAddress(self, address:str):
        self.address = address
        
    #get name of the Reader
    def getName (self , name):
        if self.name == None:
            print ('There is no name for the Reader object, please define one name!\n')
        else:
            return self.name
    
    #get the address of the Reader
    def getAddress(self):
        return self.address

    #read file as excel
    def readDataExcel(self):
        try:
            self.data = pd.read_excel(self.address)
            print("Excel file loaded successfully.")
        except Exception as e:
            print(f"An error occurred while loading the Excel file: {e}")
        return self.data
    #read file as CSV
    def readDataCSV(self):
        pass

    
    # get data 
    def getCleanedData(self):

        #give the extension of the data
        _, file_extension = os.path.splitext(self.address)
        
        # if it is .csv pass it for csv handelling
        if file_extension in ['.csv']:
            self.readDataCSV()
            Cleaner.cleanNAN(self.data)
            Cleaner.cleanDuplication(self.data)
            Cleaner.orderData(self.data)
            return self.data
        # if it is .excell  pass it for escell handelling
        elif file_extension in ['.xls', '.xlsx']:
            self.readDataExcel()
            Cleaner.cleanNAN(self.data)
            Cleaner.cleanDuplication(self.data)
            Cleaner.orderData(self.data)
            return self.data
        else:
            print("Unsupported file format.")
    def getDataAsPdDataFrame(self):

        #give the extension of the data
        _, file_extension = os.path.splitext(self.address)
        
        # if it is .csv pass it for csv handelling
        if file_extension in ['.csv']:
            self.readDataCSV()
            Cleaner.cleanNAN(self.data)
            return self.data
        # if it is .excell  pass it for escell handelling
        elif file_extension in ['.xls', '.xlsx']:
            self.readDataExcel()
            Cleaner.cleanNAN(self.data)

            return self.data
        else:
            print("Unsupported file format.")
# is there any data that is not NAN, eliminate it or other type of cleaning like remove duplication the class is static




class Cleaner:

    #this clean method considers duplication when paired monomers are the same and also the reported reactivity ratio for each monomer is the same!
    @staticmethod
    def cleanNAN(data: pd.DataFrame) -> pd.DataFrame:
        data = data.dropna()
        data.reset_index(drop=True, inplace=True)

        return data
        


    @staticmethod
    def orderData(data: pd.DataFrame) -> pd.DataFrame:
        
        # Create a column with sorted monomer pairs
        data['Sorted_Pair'] = data.apply(lambda row: tuple(sorted([row['SMILES_A'], row['SMILES_B']])), axis=1)

        # Sort the DataFrame by the sorted pairs
        data.sort_values(by='Sorted_Pair', inplace=True)

        # Optionally, if you want to remove the 'Sorted_Pair' column after sorting
        data.drop(columns=['Sorted_Pair'], inplace=True)

        data.reset_index(drop=True, inplace=True)

        return data

    
    @staticmethod
    def cleanDuplication(data: pd.DataFrame) -> pd.DataFrame:
        # Ensure the headers are correctly assigned
        Cleaner.headerAssign(data)

        # Debugging: Print DataFrame after renaming columns
    
        # Create a new column for sorted monomer pairs
        data['Monomer_Pair'] = data.apply(lambda row: tuple(sorted([row['SMILES_A'], row['SMILES_B']])), axis=1)
    
        # Create a new column for sorted reactivity ratios
        data['Reactivity_Ratio_Pair'] = data.apply(lambda row: (row['r1'], row['r2']) if row['SMILES_A'] <= row['SMILES_B'] else (row['r2'], row['r1']), axis=1)
    
        # Create a combined column to check for duplicates
        data['Duplication_Check'] = data['Monomer_Pair'] + data['Reactivity_Ratio_Pair']
    

    
        # Drop duplicates based on the combined column
        data.drop_duplicates(subset='Duplication_Check', inplace=True)
    
        # Remove the temporary columns used for duplicate checking
        data.drop(columns=['Monomer_Pair', 'Reactivity_Ratio_Pair', 'Duplication_Check'], inplace=True)

        #Reset index
        data.reset_index(drop=True, inplace=True)

        return data

    @staticmethod
    def columnEliminator (data:pd.DataFrame , *columns_to_remove ) -> pd.DataFrame:

        data.drop(columns=list(columns_to_remove), errors='ignore', inplace=True)

        return data
        
    
    @staticmethod
    def headerAssign(data: pd.DataFrame):
        data.rename(columns={'cano_smiles_sm1': 'SMILES_A', 'cano_smiles_sm2': 'SMILES_B'}, inplace=True)

    @staticmethod
    def sortByColumn(data: pd.DataFrame , column_name) -> pd.DataFrame:
        if column_name in data.columns:
            data.sort_values(by=column_name, ascending=False ,inplace=True)
            data.reset_index(drop=True, inplace=True)
            return data
        else:
            print(f"Column '{column_name}' not found in DataFrame.")
            return data

    

class getName:

    @staticmethod
    def commonNameFromSmiles(smiles):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/synonyms/JSON"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            synonyms = data['InformationList']['Information'][0]['Synonym']
            common_name = synonyms[0]  # The first synonym is often the common name
            return common_name
        else:
            return 'No Name Could be Find!'

    @staticmethod
    def IUPACNameFromSmiles(smiles):
        
        # Convert the SMILES string to a molecule object
        molecule = Chem.MolFromSmiles(smiles)
    
        # Get the IUPAC name (requires RDKit to be compiled with InChI support)
        name = Chem.MolToInchiKey(molecule)

        return name 




