import numpy as np
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from utils.data_utils import canonicalize_smiles
import pickle

class Node:
    def __init__(self, val, parent=None):
        self.smiles = val[0]
        self.ground_truth = None
        self.reagent = None
        self.plain_smiles = None
        self.childrens = []
        self.parent = parent
        self.termination = False
        self.depth = 0
        self.rank = 1
        self.score = 0
        self.pdt_smi = None
        self.rgt_smi = None
        self.rct_smi = None
        self.match = False
        self.get_ground_truth(val)
        self.get_reagent(val)
        self.get_rct_smi()
        self.get_score(val)
        self.get_depth()
        self.check_termination()

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def get_ground_truth(self, val):
        try:
            self.ground_truth = val[3]
        except:
            if self.parent:
                self.ground_truth = self.parent.ground_truth

    def get_reagent(self, val):
        try:
            self.reagent = val[4]
            if self.reagent == "No reagent\n":
                self.reagent = None
        except:
            if self.parent:
                self.reagent = self.parent.reagent

    def get_pdt_smi(self):
        if self.parent:
            self.pdt_smi = self.parent.pdt_smi
        self.pdt_smi = self.ground_truth.replace(' ', '')
        return self.pdt_smi

    def get_rgt_smi(self):
        if self.parent:
            self.rgt_smi = self.parent.rgt_smi
        if self.reagent:
            self.rgt_smi = self.reagent.replace(' ', '')
        return self.rgt_smi

    def get_rct_smi(self):
        self.rct_smi = self.smiles.replace(' ', '')
        return self.rct_smi

    def get_score(self, val):
        try:
            self.score = val[2]
        except:
            pass

    def get_rank(self, val):
        try:
            self.rank = val[1]
        except:
            pass

    def check_termination(self):
        rct_smi = canonicalize_smiles(self.get_rct_smi().replace('END', ''), trim=False, suppress_warning=True)
        pdt_smi = canonicalize_smiles(self.get_pdt_smi().replace('END', ''), trim=False, suppress_warning=True)
        if pdt_smi == rct_smi:
            self.termination = True
            self.set_match()
        elif "END" in self.smiles:
            self.termination = True
            
    def set_match(self):
        self.match = True
        if self.parent:
            self.parent.match = True
            self.parent.set_match()
        
    def get_match(self):
        return self.match

    def terminated(self):
        return self.termination
    
    def get_depth(self):
        if self.parent:
            self.depth = self.parent.get_depth() + 1
        else:
            self.depth = 0
        return self.depth
    
    def get_match_depth(self):
        if self.match:
            for child in self.childrens:
                if child.match:
                    return child.get_match_depth()
                
                
class NodeWrapper:
    def __init__(self, node, cumulative_score):
        self.node = node
        self.cumulative_score = cumulative_score

    def __lt__(self, other):
        return float(self.cumulative_score) < float(other.cumulative_score)

    def __le__(self, other):
        return float(self.cumulative_score) <= float(other.cumulative_score)

    def __gt__(self, other):
        return float(self.cumulative_score) > float(other.cumulative_score)

    def __ge__(self, other):
        return float(self.cumulative_score) >= float(other.cumulative_score)
    

def initialize_tree(src_path, tgt_path, rgt_path):
    tree = []
    with open(src_path, 'r') as src, open(tgt_path, 'r') as tgt, open(rgt_path, 'r') as rgt:
        src_lines = src.readlines()
        tgt_lines = tgt.readlines()
        rgt_lines = rgt.readlines()
        rank = [1 for i in range(1, len(src_lines) + 1)]
        score = [0 for i in range(1, len(src_lines) + 1)]
        lines = zip(src_lines, rank, score, tgt_lines, rgt_lines)
        for line in lines:
            tree.append(Node(line))
    
    return tree