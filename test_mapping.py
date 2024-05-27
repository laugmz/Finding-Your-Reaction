import re
def remove_atom_mapping(smiles):
    """
    Remove atom mapping numbers from a SMILES-like notation.
    Arguments:
    - smiles (str): The SMILES-like notation with atom mapping numbers.
    Returns:
    - smiles_without_mapping (str): The SMILES notation without atom mapping numbers.
    """
    pattern = r':\d+'
    smiles_without_mapping = re.sub(pattern, '', smiles)
    return smiles_without_mapping

def test_answer() :
    name1 = remove_atom_mapping("[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]")
    assert name1 == "[CH2]([S]([O][CH2][CH2][Br])(=[O])=[O])[CH3]", f"Expected [CH2]([S]([O][CH2][CH2][Br])(=[O])=[O])[CH3], but got {name1} instead"

    name2 = remove_atom_mapping("C(#N)C.O1CCCC1>[F:33][C:30]1[CH:29]=[C:28]([C:34]#[N:35])[C:27]([C:24]2[CH:25]=[CH:26][C:21]([CH2:20][C:17]3[C:18](=[O:19])[N:13]([C@H:10]4[CH2:11][CH2:12][C@H:4]([O:5][CH:6]([CH3:7])[C:2]([OH:3])([CH3:1])[CH3:43])[CH2:8][CH2:9]4)[C:14]4[N:15]([N:39]=[C:40]([CH3:42])[N:41]=4)[C:16]=3[CH2:36][CH2:37][CH3:38])=[CH:22][CH:23]=2)=[CH:32][CH:31]=1") 
    assert name2 == "C(#N)C.O1CCCC1>[F][C]1[CH]=[C]([C]#[N])[C]([C]2[CH]=[CH][C]([CH2][C]3[C](=[O])[N]([C@H]4[CH2][CH2][C@H]([O][CH]([CH3])[C]([OH])([CH3])[CH3])[CH2][CH2]4)[C]4[N]([N]=[C]([CH3])[N]=4)[C]=3[CH2][CH2][CH3])=[CH][CH]=2)=[CH][CH]=1", f"Expected C(#N)C.O1CCCC1>[F][C]1[CH]=[C]([C]#[N])[C]([C]2[CH]=[CH][C]([CH2][C]3[C](=[O])[N]([C@H]4[CH2][CH2][C@H]([O][CH]([CH3])[C]([OH])([CH3])[CH3])[CH2][CH2]4)[C]4[N]([N]=[C]([CH3])[N]=4)[C]=3[CH2][CH2][CH3])=[CH][CH]=2)=[CH][CH]=1, but got {name2} instead"

    name3 = remove_atom_mapping("[CH2:41]([O:40][C:38]([NH:5][C:18](=[O:20])[C@@H:17]1[CH2:21][CH2:22][CH2:23][NH:16]1)=[O:39])[C:44]1[CH:29]=[CH:28][CH:27]=[CH:26][CH:25]=1")
    assert name3 == "[CH2]([O][C]([NH][C](=[O])[C@@H]1[CH2][CH2][CH2][NH]1)=[O])[C]1[CH]=[CH][CH]=[CH][CH]=1", f"Expected [CH2]([O][C]([NH][C](=[O])[C@@H]1[CH2][CH2][CH2][NH]1)=[O])[C]1[CH]=[CH][CH]=[CH][CH]=1, but got {name3} instead"

    name4 = remove_atom_mapping(" ")
    assert name4 == " ", f"Expected , but got {name4} instead"
