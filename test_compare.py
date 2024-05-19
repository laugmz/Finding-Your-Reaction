import pandas as pd

def compare_molecule_with_data(element, string_input_mol,start_col=0, end_col=None):
    return ''.join(element.split()).lower() == ''.join(string_input_mol.split()).lower()

def test_answer() :
    same = compare_molecule_with_data("CCH", "CCH")
    assert same == True

    different = compare_molecule_with_data("CCH", "cch")
    assert different == True

    spaces = compare_molecule_with_data(" CCH ", "CCH")
    assert spaces == True

    mismatch = compare_molecule_with_data("CCH", "CCO")
    assert mismatch == False

    empty = compare_molecule_with_data("", "")
    assert empty == True

    empty2 = compare_molecule_with_data("CCO", "")
    assert empty2 == False