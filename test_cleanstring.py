import re
def clean_string(s):
    return re.sub(r'[\[\]\(\)\+\-\#]', '', s)
    
def test_answer() :
    string1 = clean_string("CHOC1CH=CHCCH=O=CHCH=1CH3CH3")
    assert string1 == "CHOC1CH=CHCCH=O=CHCH=1CH3CH3", f"Expected CHOC1CH=CHCCH=O=CHCH=1CH3CH3, but got {string1} instead"

    string2 = clean_string("[P]([CH2][NH][CH2][C]([OH])=[O])([OH])([OH])=[O]")
    assert string2 == "PCH2NHCH2COH=OOHOH=O", f"Expected PCH2NHCH2COH=OOHOH=O, but got {string2} instead"

    string3 = clean_string("[Cl][C]1[CH]=[CH][CH]=[C]([NH][S]([CH3])(=[O])=[O])[C]=1[N]1[CH2][CH2][N]([CH2][CH2][CH2][N]2[C]3[CH2][CH2][N]([S]([NH2])(=[O])=[O])[CH2][C]=3[C]([C]3[CH]=[CH][C]([C]([F])([F])[F])=[CH][CH]=3)=[N]2)[CH2][CH2]1")
    assert string3 == "ClC1CH=CHCH=CNHSCH3=O=OC=1N1CH2CH2NCH2CH2CH2N2C3CH2CH2NSNH2=O=OCH2C=3CC3CH=CHCCFFF=CHCH=3=N2CH2CH21", f"Expected ClC1CH=CHCH=CNHSCH3=O=OC=1N1CH2CH2NCH2CH2CH2N2C3CH2CH2NSNH2=O=OCH2C=3CC3CH=CHCCFFF=CHCH=3=N2CH2CH21, but got {string3} instead"