import pandas as pd
import pandas as pd
from io import StringIO

def compare_molecule_with_data(element, string_input_mol,start_col=0, end_col=None):
    return ''.join(element.split()).lower() == ''.join(string_input_mol.split()).lower()
def find_molecule_rows(dataFrame, string_input_mol, start_col=0, end_col=None):
    """
    Search through the specified range of columns in a DataFrame for the input molecule.
    
    Args:
    - dataFrame (pd.DataFrame): The DataFrame to search.
    - string_input_mol (str): The molecule to search for.
    - start_col (int): The starting column index for the search.
    - end_col (int): The ending column index for the search. If None, searches until the last column.

    Returns:
    - List[int]: A list of row indices where the molecule is found.
    """
    # Initialize a list to store the row numbers
    rows = []
    
    # Set end_col to the last column index if not provided
    if end_col is None:
        end_col = dataFrame.shape[1]
    
    # Iterate over the specified range of columns in the DataFrame
    for column_name in dataFrame.columns[start_col:end_col]:
        column = dataFrame[column_name]
        # Check if any element in the column matches the input molecule
        for index, value in column.items():
            if compare_molecule_with_data(value, string_input_mol):
                # Store the row number where the molecule is found
                rows.append(index)
    
    return rows

data = """
col1,col2,col3
H2O,C2H6,O2
CH4,H2O,N2
CO2,CH4,H2
"""

df = pd.read_csv(StringIO(data))

def test_compare_molecule_with_data():
    assert compare_molecule_with_data("H2O", "h2o") == True
    assert compare_molecule_with_data("C2 H6", "C2H6") == True
    assert compare_molecule_with_data("O2", "O 2") == True
    assert compare_molecule_with_data("H2O", "H2O2") == False

def test_find_molecule_rows_valid():
    result = find_molecule_rows(df, "H2O")
    assert result == [0, 1], f"Expected [0, 1], but got {result}"

def test_find_molecule_rows_invalid():
    result = find_molecule_rows(df, "NH3")
    assert result == [], f"Expected [], but got {result}"

def test_find_molecule_rows_specified_columns():
    result = find_molecule_rows(df, "CH4", start_col=1, end_col=3)
    assert result == [2], f"Expected [2], but got {result}"

def test_find_molecule_rows_end_col_none():
    result = find_molecule_rows(df, "CH4", start_col=1)
    assert result == [2], f"Expected [2], but got {result}"

def test_find_molecule_rows_edge_case_empty_dataframe():
    empty_df = pd.DataFrame(columns=['col1', 'col2', 'col3'])
    result = find_molecule_rows(empty_df, "H2O")
    assert result == [], f"Expected [], but got {result}"

def test_find_molecule_rows_edge_case_no_match_in_specified_range():
    result = find_molecule_rows(df, "O2", start_col=1, end_col=2)
    assert result == [], f"Expected [], but got {result}"

