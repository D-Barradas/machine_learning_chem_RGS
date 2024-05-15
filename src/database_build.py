# %%
import pandas as pd 
import pubchempy as pcp 
import requests
from datetime import date
today = date.today()
print ( "Today's date:", today)

# %%
# path = '/Users/barradd/Documents/GitHub/machine_learning_chem_RGS/data/inital-data-19-nov-23.xlsx'
path = '/Users/barradd/Documents/GitHub/machine_learning_chem_RGS/data/inital-data-07-dec-23.xlsx'

# %%
my_df = pd.read_excel(path)

# %%

def query_pubchem_by_molecular_formula(molecular_formula):
  """Queries the PubChem database using a molecular formula query.

  Args:
    molecular_formula: A string representing the molecular formula.

  Returns:
    A list of PubChem compound IDs of the compounds that match the molecular formula query.
    This is correct way to query and get the Molecular formula (replace MolecularFormula )
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula/JSON"
  """

### this is the correct way to use the fastformula 
  url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{molecular_formula}/cids/JSON"

  # Make a GET request to the URL.
  response = requests.get(url)

  # Parse the JSON response.
  json_response = response.json()

    # Check if the response is successful.
  if response.status_code == 200:
    # Parse the JSON response.
    json_response = response.json()

    # Extract the PubChem compound IDs of the compounds that match the query.
    compound_ids = []
    for compound in json_response["IdentifierList"]["CID"]:
      compound_ids.append(compound)

    # Return the list of compound IDs.
    return compound_ids
  else:
    # Return an empty list if there was an error.
    return []

  # # Extract the PubChem compound IDs of the compounds that match the query.
  # compound_ids = []
  # for compound in json_response["IdentifierList"]["CID"]:
  # # for compound in json_response:
  #   # print (compound)
  #   # compound_id = compound["CID"]
  #   # compound_ids.append(compound_id)
  #   compound_ids.append(compound)

  # return compound_ids 

# %% [markdown]
# """ Here is the prperties exaplanation """
# https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables

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
    #   'atom_stereo_count',
    #   'defined_atom_stereo_count',
    #   'defined_bond_stereo_count' , 
    #   'undefined_atom_stereo_count',
    #   'undefined_bond_stereo_count' ,
    #   'rotatable_bond_count',
    #   'isotope_atom_count',     
    ] 
# c.to_dict( properties= my_properties_list)

# %%
# my_df['substituent'].head(n=5).to_list()
my_df['molecular formula'].head(n=5).to_list()

# %%
compound_cids = {} 
for compound_raw_name in my_df['molecular formula'].head(n=5).to_list() : 
    print (compound_raw_name)
    compound_ids  = query_pubchem_by_molecular_formula(compound_raw_name)
    compound_cids[compound_raw_name] = compound_ids

# %%
all_data_frames = [] 
for m in compound_cids['C7H5NO2S']: 
    # query again the database this time with the python library
    c = pcp.Compound.from_cid(m)
    # get the information on a dictionary 
    temp_dict = c.to_dict( properties= my_properties_list) 
    # transform the dict into pandas dataframe 
    temp_df = pd.DataFrame.from_dict(data=temp_dict,orient='index').T
    # save the info in a list to postprocess into a complete dataframe
    all_data_frames.append(temp_df)

# %%
df_to_edit = pd.concat(all_data_frames) 

# %%
# total number of atoms 
df_to_edit['elements'].apply( lambda x:len(x))

# %%
# total unique atoms 
df_to_edit['elements'].apply( lambda x: len ( list(dict.fromkeys(x)) ) )


# %%
df_to_edit['all_atoms_count']= df_to_edit['elements'].apply( lambda x:len(x))
df_to_edit['all_atoms_count_unique'] = df_to_edit['elements'].apply( lambda x: len ( list(dict.fromkeys(x)) ) )


# %%
df_to_edit.drop(['bonds','elements'],axis=1 ,inplace=True)

# %%
compound_cids = {} 
for compound_raw_name in my_df['molecular formula'].to_list() : 
    # print (compound_raw_name.strip())
    compound_ids  = query_pubchem_by_molecular_formula(compound_raw_name.strip())
    compound_cids[compound_raw_name] = compound_ids

# %%
all_data_frames = []
notfound_in_database = [] 

for keys in compound_cids.keys():
    num_of_cids = len ( compound_cids[keys] )
    print ( keys.strip(),num_of_cids  )

    if num_of_cids != 0 :
        for m in compound_cids[keys]: 
            # query again the database this time with the python library
            c = pcp.Compound.from_cid(m)
            # get the information on a dictionary 
            temp_dict = c.to_dict( properties= my_properties_list) 
            # transform the dict into pandas dataframe 
            temp_df = pd.DataFrame.from_dict(data=temp_dict,orient='index').T
            # save the info in a list to postprocess into a complete dataframe
            all_data_frames.append(temp_df)
    else :
        notfound_in_database.append(keys)


# %%
df_to_edit = pd.concat(all_data_frames) 
df_to_edit['all_atoms_count']= df_to_edit['elements'].apply( lambda x:len(x))
df_to_edit['all_atoms_count_unique'] = df_to_edit['elements'].apply( lambda x: len ( list(dict.fromkeys(x)) ) )
df_to_edit.drop(['bonds','elements'],axis=1 ,inplace=True)


# %%
df_to_edit.to_excel(f'../data/preliminar_data_set_for_ML_{today}.xlsx', index=False )

# %%
out = open(f'../data/not_found_in_database_{today}.txt','w')
for id_not_found in notfound_in_database:
    out.write(id_not_found)
