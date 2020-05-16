'''
This class is to manipulate SMILES with LINGOs
'''
import re

class LINGO:
    def __init__(self):
        self.smiles = {}
        
    def print_smiles(self):
        return self.smiles
    
    def add_smiles(self,smile):
        self.smiles[smile]=[]