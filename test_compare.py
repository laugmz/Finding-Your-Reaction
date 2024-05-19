import pandas as pd

def compare_molecule_with_data(element, string_input_mol,start_col=0, end_col=None):
    return ''.join(element.split()).lower() == ''.join(string_input_mol.split()).lower()

def test_answer() :
    try1 = is_smiles("[Br-]")
    assert mol1 == True, f"Expected True, but got {try1} instead"

    try2 = is_smiles("ethylene diamine")
    assert mol2 == False, f"Expected False, but got {try2} instead"

    mol3 = is_smiles("[C]1([C]2[CH]=[CH][C]([C]3[C]4[C]5=[C]6[C](=[CH][CH]=4)[CH]=[CH][CH]=[C]6[CH]=[CH][C]5=[CH][CH]=3)=[CH][CH]=2)[C]2[C](=[CH][CH]=[CH][CH]=2)[CH]=[CH][CH]=1")
    assert try3 == True, f"Expected True, but got {try3} instead"