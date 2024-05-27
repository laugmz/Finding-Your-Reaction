import pandas as pd

def compare_molecule_with_data(element, string_input_mol, start_col=0, end_col=None):
    """
    Compares two strings regardless of spaces differences.
    
    Arguments:
    - two strings (str) to compare 
    - start_col (int): The starting column index for the search.
    - end_col (int): The ending column index for the search. If None, searches until the last column.
    Returns:
    - value (bool) : True if the two elements are the same, Flase otherwise.
    """
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