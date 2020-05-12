from LINGO import LINGO
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt

kegg_data = pd.read_csv('../data/kegg_compound_strings.tsv',sep='\t')
ref_data = pd.read_csv('../data/reference_compound_strings.tsv',sep='\t')
old_ref_data = pd.read_csv('../data/old_reference_compound_strings.tsv',sep='\t')

#%% Drop NaN values and molecules with SMILES that have less or equal to 5 char
kegg_data.dropna(inplace=True)
kegg_data['len_SMILEs'] = kegg_data['Canonical-SMILES'].apply(lambda x: len(x))
kegg_data = kegg_data[kegg_data['len_SMILEs']>5]
kegg_data['len_SMILEs'].hist(bins=100)
kegg_data.drop(['len_SMILEs'],axis=1,inplace=True)

#%% Find chemical elements with two characters

can_smiles = ''.join(list(kegg_data['Canonical-SMILES'].values))
can_smiles = ''.join(can_smiles)
all_char = list(set(can_smiles))
elements_with_2_char = list(set(re.findall("[A-Z][a-z]",can_smiles)))

#%% Define functions
def can_smiles_to_lin_smiles(smiles):
    #Takes a canonical smiles string and gets it ready for the Lingo algorithm
    if not isinstance(smiles,str): return "smiles has to be string type"
    
    lin_smiles = re.sub("[0-9]{1}","0",smiles)
    lin_smiles = re.sub("Be","E",lin_smiles)
    lin_smiles = re.sub("Li","T",lin_smiles)
    lin_smiles = re.sub("Ni","y",lin_smiles)
    lin_smiles = re.sub("Sn","$",lin_smiles)
    lin_smiles = re.sub("Cu","Q",lin_smiles)
    lin_smiles = re.sub("Co","%",lin_smiles)
    lin_smiles = re.sub("Tc","&",lin_smiles)
    lin_smiles = re.sub("As","/",lin_smiles)
    lin_smiles = re.sub("Ge","J",lin_smiles)
    lin_smiles = re.sub("Bi","?",lin_smiles)
    lin_smiles = re.sub("Cd","k",lin_smiles)
    lin_smiles = re.sub("Se","z",lin_smiles)
    lin_smiles = re.sub("Tl",">",lin_smiles)
    lin_smiles = re.sub("Sr","<",lin_smiles)
    lin_smiles = re.sub("Cs",":",lin_smiles)
    lin_smiles = re.sub("Ti","x",lin_smiles)
    lin_smiles = re.sub("Rh","j",lin_smiles)
    lin_smiles = re.sub("Au","V",lin_smiles)
    lin_smiles = re.sub("Pr","q",lin_smiles)
    lin_smiles = re.sub("Sb","!",lin_smiles)
    lin_smiles = re.sub("Pt","p",lin_smiles)
    lin_smiles = re.sub("Ag","G",lin_smiles)
    lin_smiles = re.sub("Rb","W",lin_smiles)
    lin_smiles = re.sub("La","~",lin_smiles)
    lin_smiles = re.sub("Pb","w",lin_smiles)
    lin_smiles = re.sub("Hg","X",lin_smiles)
    lin_smiles = re.sub("Cr","R",lin_smiles)
    lin_smiles = re.sub("Pe","f",lin_smiles)
    lin_smiles = re.sub("Ba","v",lin_smiles)
    lin_smiles = re.sub("Br","B",lin_smiles)
    lin_smiles = re.sub("Na","D",lin_smiles)
    lin_smiles = re.sub("Mg","M",lin_smiles)
    lin_smiles = re.sub("Al","U",lin_smiles)
    lin_smiles = re.sub("Si","I",lin_smiles)
    lin_smiles = re.sub("Ca","A",lin_smiles)
    lin_smiles = re.sub("Cl","L",lin_smiles)
    lin_smiles = re.sub("Mn","m",lin_smiles)
    lin_smiles = re.sub("Fe","E",lin_smiles)
    lin_smiles = re.sub("Zn","Z",lin_smiles)
    lin_smiles = re.sub("I","Y",lin_smiles)
    return lin_smiles

def smiles_to_LINGO(smiles, q = 4):
    #Takes a smiles string (or lin_smiles string) and the LINGOs' length q(optional)
    #and returns a list with LINGOs
    LINGOs = [smiles[i:i+q] for i in range(len(smiles)-(q - 1))]
    return LINGOs

def tanimoto(smiles1, smiles2, q = 4):
    #Takes two smiles strings (or lin_smiles string) and the LINGOs' length q(optional)
    #and returns the Tanimoto's similarity score
    lin_smiles1 = can_smiles_to_lin_smiles(smiles1)
    lin_smiles2 = can_smiles_to_lin_smiles(smiles2)
    
    LINGOs1 = smiles_to_LINGO(lin_smiles1)
    LINGOs2 = smiles_to_LINGO(lin_smiles2)
    
    LINGOs_inter = list(set(LINGOs1+LINGOs2))
    
    SUM = 0
    
    for lingo in LINGOs_inter:
        NAi = len([l for l in LINGOs1 if l == lingo])
        NBi = len([l for l in LINGOs2 if l == lingo])
        
        SUM += 1 - ((abs(NAi-NBi))/(NAi + NBi))
        
    l = len(LINGOs_inter)
    
    Tc = SUM/l
    return Tc    

def getSchiffBase(smiles):
    SchiffBase = re.findall('CN=C',smiles)
    return SchiffBase
#%% Calculate Tanimoto scores between well known reference compunds (old_ref_data) 
# and kegg database
    
q = 4

Tc_scores = np.zeros((len(old_ref_data),len(kegg_data))) 
old_ref_molecules = old_ref_data['Canonical-SMILES'].values
kegg_molecules = kegg_data['Canonical-SMILES'].values

for row in range(len(old_ref_molecules)):
    for col in range(row,len(kegg_molecules)):
        Tc_scores[row][col] = tanimoto(old_ref_molecules[row],kegg_molecules[col])

# Check the tanimoto score's distribution

all_Tc_scores = []

for row in range(len(old_ref_molecules)):
    for col in range(row,len(kegg_molecules)):
        if Tc_scores[row][col]>0.5:
            all_Tc_scores.append(Tc_scores[row][col])

plt.hist(all_Tc_scores,bins=100)

# Get only the molecule's pubchme-ids for tanimoto scores >5

relevant_Tc_scores = []
ref_list = []
kegg_list = []

for row in range(len(old_ref_molecules)):
    for col in range(row,len(kegg_molecules)):
        if Tc_scores[row][col]>0.5:
            relevant_Tc_scores.append(Tc_scores[row][col])
            print(row,col)
            pubchem_id_ref = old_ref_data[old_ref_data['Canonical-SMILES'] == old_ref_molecules[row]]['Pubchem-id'].values[0]
            pubchem_id_kegg = kegg_data[kegg_data['Canonical-SMILES'] == kegg_molecules[col]]['Pubchem-id'].values[0]
            ref_list.append(pubchem_id_ref)
            kegg_list.append(pubchem_id_kegg)

# Create DataFrame with the similar Pubchem-id molecules and the corresponding
# Tanimoto score            

ref_kegg_Tanimoto = pd.DataFrame({'Pubchem-id_ref':ref_list,\
                                  'Pubchem-id_kegg':kegg_list,\
                                  'Tanimoto_score':relevant_Tc_scores})


ref_kegg_Tanimoto.to_csv('../data/ref_kegg_Tanimoto.tsv',sep='\t',index=False)




# Create a function that finds Schiff bases in the SMILES string





















