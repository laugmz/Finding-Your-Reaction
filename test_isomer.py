from rdkit import Chem
import pandas as pd
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
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

dataFrame = pd.DataFrame({"Product1": ["C1=CC=CC=C1", "C1=CC=CC=C1"],
                                       "Product2": ["C1=CC=CCC=C1", "C1=CC=CCC=C1"]})

def test_answer() :
    empty = is_isomer(None)
    assert empty == "Invalid SMILES notation, please check input."

    in_it = is_isomer("C1=CC=CCC=C1")
    assert in_it == ("C1=CC=CCC=C1", 0, "Product2")

    not_in_it = is_isomer("C1=CC=CC=CC=C1")
    assert not_in_it == None