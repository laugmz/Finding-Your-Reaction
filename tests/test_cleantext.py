def clean_text(text):
    """
    Remove some part of reactions that that have the form |f:..| and is useless"
    Arguments: - string (str) any
    Returns: - sring (str) without the expression |f:...|. ("..." stands for a squence of numbers
    """
    return re.sub(r'\|[^|]*\|', '', text)

def test_answer() :
    without = clean_text("[CH:12]1([C:16]([NH:1][C:2]2[CH:3]=[CH:4][C:5]([CH2:8][C:9]([OH:11])=[O:10])=[CH:6][CH:7]=2)=[O:17])[CH2:15][CH2:14][CH2:13]1")
    assert without == "[CH:12]1([C:16]([NH:1][C:2]2[CH:3]=[CH:4][C:5]([CH2:8][C:9]([OH:11])=[O:10])=[CH:6][CH:7]=2)=[O:17])[CH2:15][CH2:14][CH2:13]1"

    with_ = clean_text("CC1N=C(C(O)C=C)C=CC=1.S(Cl)(Cl)=O.ClC(C1C=CC=C(C)N=1)C=C.[C-]#N.[Na+].[CH3:30][C:31]1[N:36]=[C:35]([CH:37]([CH:40]=[CH2:41])[C:38]#[N:39])[CH:34]=[CH:33][CH:32]=1.[SH2:42]>>[CH3:30][C:31]1[N:36]=[C:35]([CH:37]([CH:40]=[CH2:41])[C:38]([NH2:39])=[S:42])[CH:34]=[CH:33][CH:32]=1 |f:3.4|")
    
    assert with_ == "CC1N=C(C(O)C=C)C=CC=1.S(Cl)(Cl)=O.ClC(C1C=CC=C(C)N=1)C=C.[C-]#N.[Na+].[CH3:30][C:31]1[N:36]=[C:35]([CH:37]([CH:40]=[CH2:41])[C:38]#[N:39])[CH:34]=[CH:33][CH:32]=1.[SH2:42]>>[CH3:30][C:31]1[N:36]=[C:35]([CH:37]([CH:40]=[CH2:41])[C:38]([NH2:39])=[S:42])[CH:34]=[CH:33][CH:32]=1 "

    empty = clean_text("", "")
    assert empty == ""

    multiple = clean_text("|f:123| |f:456| |f:789|")
    assert multiple == "  "
