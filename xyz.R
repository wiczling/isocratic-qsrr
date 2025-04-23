# Canonicalize a SMILES string using RDKit
canonicalize_smiles <- function(smiles) {
  mol <- Chem$MolFromSmiles(smiles)
  if (is.null(mol)) return(NA)
  canon <- Chem$MolToSmiles(mol, canonical = TRUE)
  return(canon)
}

find_duplicate_molecules <- function(smiles_vector) {
  canonical_smiles <- sapply(smiles_vector, canonicalize_smiles)
  duplicated_indices <- which(duplicated(canonical_smiles) | duplicated(canonical_smiles, fromLast = TRUE))
  
  data.frame(
    Original_SMILES = smiles_vector[duplicated_indices],
    Canonical_SMILES = canonical_smiles[duplicated_indices],
    Index = duplicated_indices
  )
}

duplicates <- find_duplicate_molecules(smiles$smiles)
duplicates %>% distinct(Canonical_SMILES, .keep_all = TRUE)