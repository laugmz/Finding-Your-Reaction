def remove_atom_mapping(smiles):
    """
    Remove atom mapping numbers from a SMILES-like notation.
    Args:
    - smiles (str): The SMILES-like notation with atom mapping numbers.
    Returns:
    - smiles_without_mapping (str): The SMILES notation without atom mapping numbers.
    """
    # Define a regular expression pattern to match atom mapping numbers
    pattern = r':\d+'
    smiles_without_mapping = re.sub(pattern, '', smiles)
    return smiles_without_mapping

def test_answer() :
    name = remove_atom_mapping("[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]")
    assert name == "[CH2]([S]([O][CH2][CH2][Br])(=[O])=[O])[CH3]", f"Expected [CH2]([S]([O][CH2][CH2][Br])(=[O])=[O])[CH3], but got {name} instead"
