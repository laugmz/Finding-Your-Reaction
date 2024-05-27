from chemspipy import ChemSpider

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


def test_answer() :
    empty = get_molecular_weight("")
    assert empty == "Invalid SMILES string provided"

    right = get_molecular_weight("[CH2:1]([S:3]([CH2:6][C:7]1[CH:16]=[N+:15]([O-:21])[C:14]2[C:9](=[CH:10][CH:11]=[CH:12][CH:13]=2)[N+:8]=1[O-:22])(=[O:5])=[O:4])[CH3:2]")
    assert right == 268.29400000000004

    wrong_smiles = get_molecular_weight("CHNHHNI")
    assert wrong_smiles== "Invalid SMILES string provided"