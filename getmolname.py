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

def test_answer() :
    mol1 = get_molecule_name("[Cl-]")
    assert mol1 == "Chloride", f"Expected Chloride, but got {mol1} instead"

    mol2 = name_to_smiles("NN")
    assert mol2 == "Unable to retrieve molecule name.", f"Expected Unable to retrieve molecule name., but got {mol2} instead"