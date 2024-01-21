import pandas as pd
from urllib.request import urlopen 
from urllib.parse import quote 

identifiers  = ['SO2CN', 'C(CN)=C(CN)2', 'SO2C(CF3)3', 'N(O)=NCN', 'SeO2CF3']
def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Did not work'
for ids in identifiers :
    print(ids, CIRconvert(ids))