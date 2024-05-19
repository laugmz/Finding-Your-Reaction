from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import pandas as pd
import os
import random
from io import StringIO
import pytest 

current_directory = Path.cwd()
print("Current Directory:", current_directory.resolve())
example_file_path = current_directory / "1976_Sep2016_USPTOgrants_smiles.rsmi"
if example_file_path.exists():
    with example_file_path.open("r") as file:
        first_line = file.readline()
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
reactants_split = dataFrame['Reactant 1'].str.split(".", expand = True)
reactants_split.columns = [f'Reactant {i+1}' for i in range(reactants_split.shape[1])]
dataFrame = pd.concat([reactants_split, dataFrame.iloc[:, 1:]], axis=1)
#Separate the Products from each other
products_split = dataFrame['Product 1'].str.split(".", expand = True)
products_split.columns = [f'Product {i+1}' for i in range(products_split.shape[1])]
dataFrame = pd.concat([products_split, dataFrame.iloc[:, 1:]], axis=1)
dataFrame = dataFrame.loc[:, ~dataFrame.columns.duplicated(keep='first')]

null_rows = dataFrame[dataFrame["CalculatedYield"].isnull()]
dataFrame.dropna(subset=["CalculatedYield"], inplace=True)

pd.set_option('display.max_colwidth', None)
def get_random_product(dataFrame):
    filtered_columns = [pd.Series(dataFrame.iloc[:, i].dropna().unique()) for i in range(128)]
    selected_columns = pd.concat(filtered_columns, axis=0)
    random_products = selected_columns.sample(n=1)
    return random_products

def main():
    while True:
        print("Press Enter to get random products or type 'exit' to quit:")
        choice = input()
        if choice == "exit":
            print("Exiting.")
            break
        else:
            random_products = get_random_product(dataFrame)
            print(random_products)

        choice = input("Do you want to continue? (yes/no): ").strip().lower()
        if choice == "no":
            print("Exiting.")
            break

def test_get_random_product():
    result = get_random_product(dataFrame)
    assert result.iloc[0] in dataFrame.values

def test_main(monkeypatch, capsys):
    inputs = iter(["", "exit"])
    monkeypatch.setattr('builtins.input', lambda prompt='': next(inputs) if prompt == '' else prompt)

    main()

    captured = capsys.readouterr()
    output_lines = captured.out.split('\n')

    # Check if the initial prompt is in the output
    expected_prompt = "Press Enter to get random products or type 'exit' to quit:"
    assert any(expected_prompt in line.strip() for line in output_lines)
    # Check if the exit message is in the output
    assert any("Exiting." in line.strip() for line in output_lines)

if __name__ == "__main__":
    pytest.main()
