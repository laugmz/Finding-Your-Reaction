from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import re
from rdkit import RDLogger
import progressbar

def clean_text(text):
    """
    Remove some part of reactions that that have the form |f:..| and is useless"
    Arguments: - string (str) any
    Returns: - sring (str) without the expression |f:...|. ("..." stands for a squence of numbers
    """
    return re.sub(r'\|[^|]*\|', '', text)

def remove_atom_mapping(smiles):
    """
    Remove atom mapping numbers from a SMILES-like notation.
    Arguments:
    - smiles (str): The SMILES-like notation with atom mapping numbers.
    Returns:
    - smiles_without_mapping (str): The SMILES notation without atom mapping numbers.
    """
    pattern = r':\d+'
    smiles_without_mapping = re.sub(pattern, '', smiles)
    return smiles_without_mapping

def remove_percent_symbol(value):
    """
    Remove the '%' symbol from a percentage value.
    Arguments:
    - value (str): The percentage value with the '%' symbol.
    Returns:
    - value (str): The percentage value without the '%' symbol.
    """
    return value.replace('%', '')


# Create a Path object for the current directory
current_directory = Path.cwd()
# Creating a Path object for an example file that does not yet exist
example_file_path = current_directory / "1976_Sep2016_USPTOgrants_smiles.rsmi"

# Reading the contents of the file to test if we have well acess to our data
if example_file_path.exists():
    with example_file_path.open("r") as file:
        first_line = file.readline()
else:
    print("The file does not exist.")
print("The database is being processed. It will take a moment, please be patient.")   
with progressbar.ProgressBar(max_value=17, widgets=[progressbar.Percentage(), " ", progressbar.GranularBar()]) as bar:
    dataFrame= pd.read_csv("1976_Sep2016_USPTOgrants_smiles.rsmi", delimiter='\t',low_memory=False)
    bar.update(1)
    dataFrame = dataFrame.astype(str)
    dataFrame["ReactionSmiles"] = dataFrame["ReactionSmiles"].apply(lambda x: clean_text(x))
    bar.update(2)
#Delete columns that are not the reactions or the yield
    columns_to_delete = ["PatentNumber", "ParagraphNum", "Year", "TextMinedYield"]
    dataFrame.drop(columns=columns_to_delete, inplace=True)
    bar.update(3)
#Separate the Reactants from the Products
    new_columns = dataFrame['ReactionSmiles'].str.split('>', expand=True)
    new_columns = new_columns.rename(columns={0: 'Reactant 1', 2 : "Product 1"})
    bar.update(4)
    dataFrame = pd.concat([new_columns, dataFrame.iloc[:, 1:]], axis=1)
    dataFrame.drop(columns=1, inplace=True)
    bar.update(5)
#Separate the reactants from each other
    reactants_split = dataFrame['Reactant 1'].str.split("\.", expand = True)
    bar.update(6)
    reactants_split.columns = [f'Reactant {i+1}' for i in range(reactants_split.shape[1])]
    dataFrame = pd.concat([reactants_split, dataFrame.iloc[:, 1:]], axis=1)
    bar.update(7)
#Separate the Products from each other
    products_split = dataFrame['Product 1'].str.split("\.", expand = True)
    bar.update(9)
    products_split.columns = [f'Product {i+1}' for i in range(products_split.shape[1])]
    dataFrame = pd.concat([products_split, dataFrame.iloc[:, 1:]], axis=1)
    bar.update(10)
    dataFrame = dataFrame.loc[:, ~dataFrame.columns.duplicated(keep='first')]
    bar.update(11)

    dataFrame = dataFrame[dataFrame["CalculatedYield"] != 'nan']
    bar.update(12)

#Remove atom mapping
    dataFrame = dataFrame.astype(str)
    bar.update(14)
    columns_to_process = dataFrame.columns[:-1]
    bar.update(15)
    for column in columns_to_process:
        dataFrame[column] = dataFrame[column].apply(remove_atom_mapping)
        bar.update(16)
#Remove the percentage symbol from the yield's column
    dataFrame['CalculatedYield'] = dataFrame['CalculatedYield'].apply(remove_percent_symbol)
    bar.update(17)
print("Processing of the data complete!")

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
