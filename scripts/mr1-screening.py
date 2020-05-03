from LINGO import LINGO
import pandas as pd

data = pd.read_csv('../data/kegg_compound_strings.tsv')

smiles = ['Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl','CC(=O)O','CC(C)(C)O',\
          'C(C(CO[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-]'] 

l = LINGO()
l.add_smiles(smiles[0])
