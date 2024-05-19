from rdkit import Chem

def is_smiles(smiles):
    """
    Check if a string represents a valid SMILES notation.
    Args:
    - smiles (str): The string to check.
    Returns:
    - is_valid (bool): True if the string is a valid SMILES notation, False otherwise.
    """
    if Chem.MolFromSmiles(smiles)== None:
        return False
    else:
        return True

def test_answer() :
    mol1 = is_smiles("Br-")
    assert mol1 == True, f"Expected True, but got {mol1} instead"

    mol2 = is_smiles("ethylene diamine")
    assert mol2 == False, f"Expected False, but got {mol2} instead"

    mol3 = is_smiles("[C]1([C]2[CH]=[CH][C]([C]3[C]4[C]5=[C]6[C](=[CH][CH]=4)[CH]=[CH][CH]=[C]6[CH]=[CH][C]5=[CH][CH]=3)=[CH][CH]=2)[C]2[C](=[CH][CH]=[CH][CH]=2)[CH]=[CH][CH]=1")
    assert mol3 == True, f"Expected True, but got {mol1} instead"