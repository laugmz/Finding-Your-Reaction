from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import pandas as pd
import os
current_directory = Path.cwd()
print("Current Directory:", current_directory.resolve())
example_file_path = current_directory / "1976_Sep2016_USPTOgrants_smiles.rsmi"
if example_file_path.exists():
    with example_file_path.open("r") as file:
        first_line = file.readline()
        print(first_line)
else:
    print("The file does not exist.")
dataFrame= pd.read_csv("1976_Sep2016_USPTOgrants_smiles.rsmi", delimiter='\t',low_memory=False)
#Delete columns that are not the reactions or the yield
columns_to_delete = ["PatentNumber", "ParagraphNum", "Year", "TextMinedYield"]
dataFrame.drop(columns=columns_to_delete, inplace=True)
#Separate the Reactants from the Products
new_columns = dataFrame['ReactionSmiles'].str.split('>', expand=True)
new_columns = new_columns.rename(columns={0: 'Reactant 1', 2 : "Product 1"})
dataFrame = pd.concat([new_columns, dataFrame.iloc[:, 1:]], axis=1)
dataFrame.drop(columns=1, inplace=True)
#Separate the reactants from each other
reactants_split = dataFrame['Reactant 1'].str.split("\.", expand = True)
reactants_split.columns = [f'Reactant {i+1}' for i in range(reactants_split.shape[1])]
dataFrame = pd.concat([reactants_split, dataFrame.iloc[:, 1:]], axis=1)
#Separate the Products from each other
products_split = dataFrame['Product 1'].str.split("\.", expand = True)
products_split.columns = [f'Product {i+1}' for i in range(products_split.shape[1])]
dataFrame = pd.concat([products_split, dataFrame.iloc[:, 1:]], axis=1)
dataFrame = dataFrame.loc[:, ~dataFrame.columns.duplicated(keep='first')]

null_rows = dataFrame[dataFrame["CalculatedYield"].isnull()]
dataFrame.dropna(subset=["CalculatedYield"], inplace=True)

pd.set_option('display.max_colwidth', None)
def main():
    
    while True:
        choice=input("Press Enter to get random products or type 'exit' to quit: ")
        if choice == "exit":
            print("Exiting.")
            break
        else: 
            filtered_columns = [pd.Series(dataFrame.iloc[:, i].dropna().unique()) for i in range(128)]
            selected_columns = pd.concat(filtered_columns, axis=0)
            random_products = selected_columns.sample(n=1)
            print(random_products)

        
        choice = input("Do you want to continue? (yes/no): ").strip().lower()
        if choice == "no":
            print("Exiting.")
            break
            
if __name__ == "__main__":
    main()

def test_answer() :
    for i in range(100):
        result = pick_random_entry(dataFrame)
        assert result.name in dataFrame.index
        assert result["Product 45"] in dataFrame["Product 45"].values
        assert result["Product 3"] in dataFrame["Product 3"].values
        assert result["Product 102"] in dataFrame["Product 102"].values
        assert result["Product 67"] in dataFrame["Product 67"].values
        assert result["Product 32"] in dataFrame["Product 32"].values
        assert result["Product 85"] in dataFrame["Product 85"].values
        assert result["Product 13"] in dataFrame["Product 13"].values
        assert result["Product 118"] in dataFrame["Product 118"].values



















