import pytest
from chemspipy import ChemSpider

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

class MockResult:
    def __init__(self, common_name):
        self.common_name = common_name

@pytest.fixture
def mock_chemspider(mocker):
    mock_chemspider_class = mocker.patch('test_namefromspider.ChemSpider')
    mock_instance = mock_chemspider_class.return_value
    return mock_instance

def test_get_molecule_name_valid(mock_chemspider):
    mock_chemspider.simple_search.return_value = [MockResult('Aspirin')]
    result = get_molecule_name('CC(=O)OC1=CC=CC=C1C(=O)O')
    assert result == 'Aspirin'

def test_get_molecule_name_invalid(mock_chemspider):
    mock_chemspider.simple_search.return_value = []
    result = get_molecule_name('InvalidSMILES')
    assert result == 'Unable to retrieve molecule name.'

def test_get_molecule_name_no_results(mock_chemspider):
    mock_chemspider.simple_search.return_value = []
    result = get_molecule_name('CC')
    assert result == 'Unable to retrieve molecule name.'