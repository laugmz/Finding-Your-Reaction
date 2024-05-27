import pubchempy as pcp

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

def test_answer() :
    mol1 = name_to_smiles("Fluoride")
    assert mol1 == "[F-]"

    mol2 = name_to_smiles("Ethanol")
    assert mol2 == "CCO"

    mol3 = name_to_smiles("Diamine")
    assert mol3 == "NN" 

    mol4 = name_to_smiles(" ")
    assert mol4 == None 
