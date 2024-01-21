import rdkit.Chem as Chem

def convert_formula_to_smiles(formula):
  """Converts a chemical formula to SMILES format.

  Args:
    formula: A chemical formula in string format.

  Returns:
    A SMILES string for the given chemical formula.
  """

  mol = Chem.MolFromSmiles(formula)
  return Chem.MolToSmiles(mol)

def get_all_properties(mol):
  """Gets all possible properties for a given molecule.

  Args:
    mol: A RDKit molecule object.

  Returns:
    A dictionary of all possible properties for the given molecule.
  """

  properties = {}
  for prop in Chem.Descriptors.MolDescriptors.MolDescriptors:
    properties[prop.name] = prop.calc(mol)
  return properties

def convert_properties_to_dataframe(properties):
  """Converts a dictionary of properties to a Pandas DataFrame.

  Args:
    properties: A dictionary of properties.

  Returns:
    A Pandas DataFrame containing the properties.
  """

  import pandas as pd
  df = pd.DataFrame(properties, index=[0])
  return df

def main():
  # Convert a chemical formula to SMILES format.
  smiles = convert_formula_to_smiles("C6H6")

  # Get all possible properties for the molecule.
  mol = Chem.MolFromSmiles(smiles)
  properties = get_all_properties(mol)

  # Convert the properties to a Pandas DataFrame.
  df = convert_properties_to_dataframe(properties)

  # Print the DataFrame.
  print(df)

if __name__ == "__main__":
  main()
