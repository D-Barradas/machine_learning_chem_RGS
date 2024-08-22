# %%
import sys
# from datetime import date
import pubchempy as pcp 
import pandas as pd
from rdkit import Chem
import time


# %%
# today = date.today()
# print ( "Today's date:", today)

# %%
# all_cids = [ x for x in open("../data/compunds_cids.txt", "r")]
all_cids = [ x for x in open("../data/missing_cid.txt", "r")]

# %%
# def get_the_rdkit_features(df_to_edit):
#     mol = Chem.MolFromSmiles(df_to_edit['canonical_smiles'][0])
#     df_to_edit['Num_of_Rings'] = mol.GetRingInfo().NumRings()
#     df_to_edit["heavy_atoms"] = mol.GetNumAtoms()
#     return df_to_edit,mol
    
def get_the_rdkit_features(df_to_edit):
    mol = Chem.MolFromSmiles(df_to_edit['canonical_smiles'][0])
    try:
        df_to_edit['Num_of_Rings'] = mol.GetRingInfo().NumRings()
    except AttributeError:
        # Handle the case where mol has no attribute 'GetRingInfo'
        print("Warning: Mol has no attribute 'GetRingInfo' for molecule:", df_to_edit['canonical_smiles'][0])
        print ("Setting Num_of_Rings to 0")
        df_to_edit['Num_of_Rings'] = 0

    df_to_edit['heavy_atoms'] = mol.GetNumAtoms()
    return df_to_edit, mol 

# %%

def get_the_type_of_bond(mol):
    counter_ar, counter_single, counter_double, counter_triple = 0, 0, 0, 0

    try:
        for num, m in enumerate(mol.GetAtoms()):
            # Retrieve bond information for the current atom
            bond_type = str(mol.GetBonds()[num].GetBondType())

            # Update bond type counters
            if bond_type == "AROMATIC":
                counter_ar += 1
            elif bond_type == "SINGLE":
                counter_single += 1
            elif bond_type == "DOUBLE":
                counter_double += 1
            else:
                counter_triple += 1

    except IndexError:
        # Handle out-of-range error
        print(f"Error: Index out of range when accessing bond information, CIDS: {my_cids}")
        # exit()
        # return None

    # Create and return a DataFrame containing the bond type counts
    df_temp = pd.DataFrame([counter_ar, counter_single, counter_double, counter_triple],
                           index=['aromatic_bond', 'single_bond', 'double_bond', 'triple_bond'])
    return df_temp.T

# def get_the_type_of_bond(mol):
#     counter_ar , counter_single, counter_double, counter_triple = 0,0,0,0
#     for num ,m in enumerate ( mol.GetAtoms() ):
#         # print ( num, m.GetSymbol(), mol.GetAtomWithIdx(num).IsInRing(), mol.GetBonds()[num].GetBondType() )
#         if str ( mol.GetBonds()[num].GetBondType() ) == "AROMATIC":
#             counter_ar+=1
#         elif str ( mol.GetBonds()[num].GetBondType() )  == "SINGLE":
#             counter_single+=1
#         elif str ( mol.GetBonds()[num].GetBondType() )== "DOUBLE":
#             counter_double+=1
#         else :
#             counter_triple+=1
#     df_temp = pd.DataFrame([counter_ar , counter_single, counter_double, counter_triple],
#              index=['aromatic_bond','single_bond', 'double_bond', 'triple_bond' ])
#     return df_temp.T
    # print (num, m.GetAtomicNum() , m.GetSymbol() , mol.GetAtomWithIdx(num).IsInRing() ,  mol.GetAtomWithIdx(num).IsInRingSize(5) , mol.GetAtomWithIdx(num).IsInRingSize(6)

# %%
my_properties_list = [ 
    'cid',
    'molecular_formula',
    'molecular_weight',
    'exact_mass',
    'monoisotopic_mass' ,
    'charge' ,
    'heavy_atom_count',
    'h_bond_donor_count',
    'h_bond_acceptor_count',
    'xlogp',
    'tpsa',
    'canonical_smiles',
    'complexity',
    'covalent_unit_count',
    'bonds',
    'elements'
                      ] 

# %%
line_number = int(sys.argv[1])
my_cids = all_cids[line_number]

# %%
# c = pcp.Compound.from_cid(my_cids)
def get_compound_info(my_cids):
    while True:
        try:
            c = pcp.Compound.from_cid(my_cids)
            return c
        except pcp.PubChemHTTPError as e:
            if e.msg == "PUGREST.ServerBusy":
        # except urllib.error.HTTPError as e:
        #     if e.code == 503:
                print("Server busy, waiting for 10 seconds...")
                time.sleep(10)  # Wait for 10 seconds before retrying
            else:
                raise e
c = get_compound_info(my_cids)
# get the information on a dictionary 
temp_dict = c.to_dict( properties= my_properties_list) 
# transform the dict into pandas dataframe 
df_to_edit = pd.DataFrame.from_dict(data=temp_dict,orient='index').T

# %%
element_type = pd.get_dummies( pd.Series ( df_to_edit['elements'].iloc[0]) ,dtype=float ).sum()

# %%
df_to_edit['all_atoms_count']= df_to_edit['elements'].apply( lambda x:len(x))
df_to_edit['all_atoms_count_unique'] = df_to_edit['elements'].apply( lambda x: len ( list(dict.fromkeys(x)) ) )
# df_to_edit.drop(['bonds','elements'],axis=1 ,inplace=True)

# %%
# pd.concat( [df_to_edit ,element_type.to_frame().T],axis=1)

# %%
# for m in mol.GetAromaticAtoms():
#     print (m.GetSymbol())

# %%
df_to_edit_2 , mol  = get_the_rdkit_features(df_to_edit=df_to_edit)
df_to_edit_3 = get_the_type_of_bond(mol)

# %%
final_dataframe = pd.concat( [df_to_edit_2, df_to_edit_3,element_type.to_frame().T],axis=1)
final_dataframe.drop(['bonds','elements'],axis=1 ,inplace=True)

# %%
final_dataframe.to_csv(f"../data/preliminar_data_{my_cids.strip()}.csv", index=None )


