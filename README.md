![Project Logo](https://github.com/laugmz/Amazing-Project/blob/main/reaction.png)



# Find your reaction
Repository used to create the project of Lilou, Eva, and Laura.

##  General use of the package
This package is designed to help one find the formation reaction which bears the highest yield for a wanted molecule through a database.

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
Three types of functions have been used throughout the package in order to get the results.

1. Functions to rearrange the database:
   
`remove_atom_mapping

remove_percent_symbol

clean_string`


## Licence 
MIT Licence
