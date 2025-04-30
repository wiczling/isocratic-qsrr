# filename: mcs_single.py
from rdkit import Chem
from rdkit.Chem import rdFMCS

# Setup parameters
mcs_params = rdFMCS.MCSParameters()
mcs_params.AtomCompareParameters.MatchValences = True
mcs_params.AtomCompareParameters.RingMatchesRingOnly = True
mcs_params.AtomCompareParameters.CompleteRingsOnly = True
mcs_params.AtomCompare = rdFMCS.AtomCompare.CompareElements
mcs_params.BondCompare = rdFMCS.BondCompare.CompareOrderExact
mcs_params.BondCompareParameters.CompleteRingsOnly = True
mcs_params.BondCompareParameters.RingMatchesRingOnly = True
mcs_params.Timeout = 20
mcs_params.Verbose = False
mcs_params.Threshold = 1.0

def compute_mcs_single(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return None
    
    mcs_result = rdFMCS.FindMCS([mol1, mol2], parameters=mcs_params)
    
    return mcs_result.numAtoms
