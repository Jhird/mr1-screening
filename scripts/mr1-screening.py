from LINGO import LINGO
import pandas as pd

kegg_data = pd.read_csv('../data/kegg_compound_strings.tsv',sep='\t')
ref_data = pd.read_csv('../data/reference_compound_strings.tsv',sep='\t')
met_ref= "135409352"
smiles = ['C1=C(N=C2C(=O)NC(=NC2=N1)N)C=O'] 
#'https://pubchem.ncbi.nlm.nih.gov/compound/6-Formylpterin#section=InChI-Key'
l = LINGO()
l.add_smiles(smiles[0])
ref_data[ref_data['Pubchem-id']==135409352]
