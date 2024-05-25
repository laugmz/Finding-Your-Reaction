# Find your reaction
![Project Logo](https://github.com/laugmz/Amazing-Project/blob/main/reaction.png)

Repository used to create the project of Lilou, Eva, and Laura.

##  General use of the package
This package is designed to help one find the formation reaction which bears the highest yield for a wanted molecule through a database.

## Features
- Transformation of the data into dataframe.
- Cleaning the dataframe: remove the atom mapping and percent symbol.
- Generating an isomer dataframe in order to find isomers of the input mol if needed.
- Search in the dataframe the molecule entered by the user and eventually its isomer.
- If its found in it, returns the reaction of formation of the molecule.
- Informations about the product: UPAC name, molecular weight, 3D visualization. 

## Downloading The Database
The databse needs to be dowloaded on the user's computer directly before the package can be used.
First, `clone` the Github repository on your computer.
```
git clone https://github.com/laugmz/Amazing-Project.git
```

Then accessing the following link (https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873), download the whole folder.
Open the folder and unzip the file named "1976_Sep2016_USPTOgrants_smiles.7z". Move the created file named "1976_Sep2016_USPTOgrants_smiles.rsmi" to the Amazing-project folder. 
It is vital to not change the name of the file for the code to work.

## Installation
Install the following tools in order to run the code.

Install `pandas`:
```
pip install pandas==2.1.4
```
Install `matplotlib`:
```
pip install matplotlib
```
Install `jupyter lab`:
```
pip install jupyter lab
```
Install `chemspipy`:
```
pip install chemspipy
```
Install `rdkit`:
```
pip install rdkit
```
Install `pathlib`:
```
pip install pathlib:
```
Install `pubchempy`
```
pip install pubchempy
```

## Functions
The functions that have been used throughout the package can be grouped in four different categories.

1. Functions for the rearrangement the database:

|Function name | Function description | Test function |
|-----------------|-----------------|-----------------|
| `remove_atom_mapping` | Removes the atom mapping numbers from a SMILES-like notation|`test_mapping`|
| `remove_percent_symbol`|Removes the "%" symbol from a percentage value| `test_removepercentage`|
|`clean_string`|Removes brackets, paranthesis and plus/minus signs|`test_cleanstring`|
   

2. Function for the input of a molecule:

|Function name | Function description | Test function |
|-----------------|-----------------|-----------------|
| `main` | Inputs a random product|`test_randmoness`| 

3. Functions for the analysis of the database:

|Function name | Function description | Test function |
|-----------------|-----------------|-----------------|
| `is_smiles` | Checks if a string represents a valid SMILES notation|`test_issmiles.py`|
|`name_to_smiles`|Converts a molecule name to a SMILES notation using PubChemPy's PubChem database| `test_nametosmiles.py`|
|`compare_molecule_with_data`|Compares the inputted molecule to the ones in the database|`test_compare.py`|
|`find_molecule_rows`|Searches through the specified range of columns in a DataFrame for the input molecule|`test_findrow.py`|
|`generate_permutations`|Generates all permutations of the inputted string |`test_permutations.py`|

4. Functions for the visualization of the results

|Function name | Function description | Test function |
|-----------------|-----------------|-----------------|
| `get_molecule_name`|Gets the common name of a molecule from its SMILES representation using ChemSpider|`test_namefromspider.py`|
|`plot_molecule_3D`|Plots a molecule in 3D with different colors for different types of atoms and bonds between atoms| `test_plot3D.py`|
|`plot_bond`|Plots bonds between atoms with different styles for single and double bonds (the test tests both functions) |`test_plot3D.py`|
![Test Results](https://img.shields.io/badge/tests-0%25%20passed-red)

## Licence 
MIT Licence
