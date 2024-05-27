from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import re
from rdkit import RDLogger

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

def main():
    
    while True:
        choice=input("Press Enter to get random products or type 'exit' to quit: ")
        if choice == "exit":
            print("Exiting.")
            break
        else: 
            filtered_columns = [pd.Series(dataFrame.iloc[:, i].dropna().unique()) for i in range(53)]
            selected_columns = pd.concat(filtered_columns, axis=0)
            random_products = selected_columns.sample(n=1)
            print(random_products)

        
        choice = input("Do you want another isomer ? (yes/no): ").strip().lower()
        if choice == "no":
            print("Exiting.")
            break

def is_smiles(smiles):
    """
    Check if a string represents a valid SMILES notation.
    Arguments:
    - smiles (str): The string to check.
    Returns:
    - is_valid (bool): True if the string is a valid SMILES notation, False otherwise.
    """
    if Chem.MolFromSmiles(smiles) == None:
        return False
    else:
        return True

def name_to_smiles(molecule_name):
    """
    Convert a molecule name to a SMILES notation using PubChemPy's PubChem database.
    Arguments:
    - molecule_name (str): The name of the molecule.
    Returns:
    - smiles (str): The SMILES notation of the molecule, or None if retrieval fails.
    """
    try:
        compound = pcp.get_compounds(molecule_name, 'name')
        if compound:
            return compound[0].canonical_smiles
        else:
            print("Error: Unable to retrieve molecule information. Please try with the SMILE of the molecule")
            return None
    except:
        print("Error: Unable to retrieve molecule information. Please try again")
        return None

def compare_molecule_with_data(element, string_input_mol, start_col=0, end_col=None):
    """
    Compares two strings regardless of spaces differences.
    
    Arguments:
    - two strings (str) to compare 
    - start_col (int): The starting column index for the search.
    - end_col (int): The ending column index for the search. If None, searches until the last column.
    Returns:
    - value (bool) : True if the two elements are the same, Flase otherwise.
    """
    return ''.join(element.split()).lower() == ''.join(string_input_mol.split()).lower()


def find_molecule_rows(dataFrame, string_input_mol, start_col=0, end_col=None):
    """
    Search through the specified range of columns in a DataFrame for the input molecule.
    
    Arguments:
    - dataFrame (pd.DataFrame): The DataFrame to search.
    - string_input_mol (str): The molecule to search for.
    - start_col (int): The starting column index for the search.
    - end_col (int): The ending column index for the search. If None, searches until the last column.

    Returns:
    - List[int]: A list of row indices where the molecule is found.
    """
    total_elements = (end_col - start_col) * len(dataFrame)
    with progressbar.ProgressBar(max_value=total_elements, widgets=[progressbar.Percentage(), " ", progressbar.GranularBar()]) as bar:
    # Initialize a list to store the row numbers
        rows = []
    # Set end_col to the last column index if not provided
        if end_col is None:
            end_col = dataFrame.shape[1]

    # Total number of elements to check
        current_element = 0
    # Iterate over the specified range of columns in the DataFrame
        for column_name in dataFrame.columns[start_col:end_col]:
            column = dataFrame[column_name]  # Get the actual column data
            for index, value in column.items():
                if compare_molecule_with_data(value, string_input_mol):
                    rows.append(index)
                current_element += 1
                bar.update(current_element)

    return rows

def is_isomer(smiles_given):
    """
    Checks if an input has an isomer in a given dataFrame.
    
    Arguments:
    - SMILES notation of a molecule (str)
    Returns:
    - Error message (str) if input isn't a SMILES notation
    - None if no isomer found
    If found:
    - Isomer in Smiles notation (str)
    - Indexes (int) of isomer's column and row
    """
    if smiles_given is None :
        return "Invalid SMILES notation, please check input."
    else:
        mol1 = Chem.MolFromSmiles(smiles_given)
        #Give molecular formula
        formula1 = rdMolDescriptors.CalcMolFormula(mol1)
    
    for idx, row in dataFrame.iterrows():
        for col in dataFrame.columns:
            if "Product" in col:
                smiles2 = row[col]
                if smiles2 == "None":
                    continue
                else:
                    mol2 = Chem.MolFromSmiles(smiles2)
                    if mol2 is None:
                        continue
                    formula2 = rdMolDescriptors.CalcMolFormula(mol2)
                    if formula1 == formula2:
                        return smiles2, idx, col

def get_molecule_name(smiles):
    """
    Get the common name of a molecule from its SMILES representation using ChemSpider.
    Args:
    - smiles (str): SMILES representation of the molecule.
    Returns:
    - name (str): Common name of the molecule.
    """
    cs = ChemSpider('your_api_key_here') 
    results = cs.simple_search(smiles)
    if results:
        return results[0].common_name
    else:
        return "Unable to retrieve molecule name."


# To print the molecular weight
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecular_weight(smiles):
    """
    Calculate the molecular weight of a molecule given its SMILES string.
    
    Args:
    smiles (str): SMILES string of the molecule.
    
    Returns:
    float: Molecular weight of the molecule.
    """
    # Convert the SMILES string to an RDKit molecule object
    molecule = Chem.MolFromSmiles(smiles)
    
    # Check if the molecule is valid
    if molecule is None:
        raise ValueError("Invalid SMILES string provided")
    
    # Calculate the molecular weight
    molecular_weight = Descriptors.MolWt(molecule)
    
    return molecular_weight

def plot_molecule_3D(smiles):
    """
    Plot a molecule in 3D with different colors for different types of atoms and bonds between atoms.
    Args:
    - smiles (str): SMILES representation of the molecule.
    """
    # Generate RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Unable to generate molecule from SMILES.")
        return

