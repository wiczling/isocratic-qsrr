from rdkit import Chem
from rdkit.Chem import rdFMCS, rdmolops

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

def extract_diff_fragments(mol, match_atoms):
    """Remove matched atoms and return the atom counts of remaining fragments."""
    all_atoms = set(range(mol.GetNumAtoms()))
    unmatched = all_atoms - set(match_atoms)

    if not unmatched:
        return 0

    emol = Chem.EditableMol(mol)
    # Remove atoms in descending order to avoid index shifts
    for idx in sorted(match_atoms, reverse=True):
        emol.RemoveAtom(idx)

    try:
        frag_mol = emol.GetMol()
        fragments = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    except Exception:
        fragments = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)

    return max(frag.GetNumAtoms() for frag in fragments) if fragments else 0

def compute_mcs_and_fragment_diff(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return None

    mcs_result = rdFMCS.FindMCS([mol1, mol2], parameters=mcs_params)
    smarts = mcs_result.smartsString
    common_mol = Chem.MolFromSmarts(smarts)

    match1 = mol1.GetSubstructMatch(common_mol)
    match2 = mol2.GetSubstructMatch(common_mol)

    max_unmatched_atoms = max(
        extract_diff_fragments(mol1, match1),
        extract_diff_fragments(mol2, match2)
    )

    return max_unmatched_atoms
