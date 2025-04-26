### MCS 
extract_common_fragments <- function(mol, match_atoms) {
  
  editable_mol <- Chem$RWMol(mol)
  mcs_atom_indices <- as.integer(match_atoms)
  
  # Remove all atoms not in MCS (do it in reverse to preserve indices)
  all_indices <- 0:(mol$GetNumAtoms() - 1)
  atoms_to_remove <- setdiff(all_indices, mcs_atom_indices)
  
  for (idx in rev(atoms_to_remove)) {
    editable_mol$RemoveAtom(as.integer(idx))
  }
  
  # Get final molecule
  common_mol <- editable_mol
  
  num_atoms <- common_mol$GetNumAtoms()
  for (i in seq(0, num_atoms - 1)) {
    atom <- common_mol$GetAtomWithIdx(i)
    if (!atom$IsInRing() && atom$GetIsAromatic()) {
      atom$SetIsAromatic(FALSE)
    }
  }
  
  res <- try({Chem$SanitizeMol(common_mol)}, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    common_smiles <- Chem$MolToSmiles(common_mol)
    common_smiles = paste0("Err:",common_smiles)
  } else {
    Chem$SanitizeMol(common_mol)
    common_smiles <- Chem$MolToSmiles(common_mol)
  }
  
  return(common_smiles)
}

extract_attachemnt_string <- function(mol, match_atoms) {
  
  all_atoms <- seq_len(mol$GetNumAtoms()) - 1L
  diff_atoms <- setdiff(all_atoms, match_atoms)
  
  attachment_strings <- c()
  for (atom_idx in match_atoms) {
    atom <- mol$GetAtomWithIdx(atom_idx)
    atom_symbol <- atom$GetSymbol()
    
    for (bond in atom$GetBonds()) {
      nbr_idx <- bond$GetOtherAtomIdx(atom_idx)
      if (!(nbr_idx %in% match_atoms)) {
        nbr_atom <- mol$GetAtomWithIdx(nbr_idx)
        
        
        
        nbr_symbol <- nbr_atom$GetSymbol()
        bond_type <- as.character(bond$GetBondType())
        
        attachment_strings <- c(
          attachment_strings,
          sprintf("%d:%s %d:%s %s",
                  atom_idx, atom_symbol,
                  nbr_idx, nbr_symbol,
                  bond_type)
        )
      }
    }
  }
  return(paste(attachment_strings, collapse = "."))
}

extract_diff_fragments <- function(mol, match_atoms) {
  
  all_atoms <- seq_len(mol$GetNumAtoms()) - 1L
  diff_atoms <- setdiff(all_atoms, match_atoms)
  
  if (length(diff_atoms) == 0) return("")
  rw_mol <- Chem$RWMol(mol)
  
  # Set atoms to remove (in reverse order to preserve indices)
  for (idx in sort(unlist(match_atoms), decreasing = TRUE)) {
    rw_mol$RemoveAtom(idx)
  }
  
  res <- try({
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE)
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE, sanitizeFrags = FALSE)
  } else {
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE)
  }
  
  frag_smiles <- sapply(mol_frags, function(frag) as.character(Chem$MolToSmiles(frag)))
  return(paste(frag_smiles, collapse = "."))
}

compare_smiles_pair <- function(smiles1, smiles2, .mcs_params=mcs_params) {
  
  mol1 <- Chem$MolFromSmiles(smiles1)
  mol2 <- Chem$MolFromSmiles(smiles2)
  
  if (is.null(mol1) || is.null(mol2)) {
    return(data.frame(
      smile1 = smiles1,
      smile2 = smiles2,
      common = NA,
      to_remove = NA,
      to_remove_str = NA,
      to_add = NA,
      to_add_str = NA,
      structure = "<em>Invalid SMILES</em>"
    ))
  }
  
  mcs_result <- rdFMCS$FindMCS(list(mol1, mol2), parameters = mcs_params)
  smarts <- mcs_result$smartsString
  common_mol <- Chem$MolFromSmarts(smarts)
  
  # Get atom matches for highlighting
  match1 <- mol1$GetSubstructMatch(common_mol)
  match2 <- mol2$GetSubstructMatch(common_mol)
  
  # Compute unmatched SMILES
  to_remove <- extract_diff_fragments(mol1, match1)
  to_add <- extract_diff_fragments(mol2, match2)
  
  to_remove_str <- extract_attachemnt_string(mol1, match1)
  to_add_str <- extract_attachemnt_string(mol2, match2)
  
  common_smile1 <- extract_common_fragments(mol1, match1)
  common_smile2 <- extract_common_fragments(mol2, match2)

  if (!str_detect(common_smile1, fixed("Err:"))) {
    common_smile = common_smile1
  } else if (!str_detect(common_smile2, fixed("Err:"))) {
    common_smile = common_smile2
  } else {
    common_smile = common_smile1
  }
  
  # Generate 2D coordinates
  AllChem$Compute2DCoords(mol1)
  AllChem$Compute2DCoords(mol2)
  
  # Draw image
  img <- Draw$MolsToGridImage(
    list(mol1, mol2),
    highlightAtomLists = list(match1, match2),
    subImgSize = tuple(300L, 300L),
    legends = list("Mol 1", "Mol 2")
  )
  
  # Save to temp file and encode as base64 HTML img tag
  img_file <- tempfile(fileext = ".png")
  img$save(img_file)
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  
  gc()
  py_run_string("import gc; gc.collect()")
   img_magick <- image_read(img_file)
   print(img_magick)  # Displays in RStudio viewer or default graphics device
  
  # Return data frame row
  data.frame(
    smile1 = smiles1,
    smile2 = smiles2,
    common = common_smile,
    to_remove = to_remove,
    to_remove_str= to_remove_str,
    to_add = to_add,
    to_add_str=to_add_str,
    structure = sprintf('<img src="%s" style="width:600px; height:auto;"/>', img_base64), 
    stringsAsFactors = FALSE)
}


make_similarity_matrix_fun <- function(lower_tri_df) {
  n <- max(lower_tri_df$row, lower_tri_df$col)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)
  for (k in 1:nrow(lower_tri_df)) {
    i <- lower_tri_df$row[k]
    j <- lower_tri_df$col[k]
    similarity <- lower_tri_df$similarity[k]
    similarity_matrix[i, j] <- similarity
    similarity_matrix[j, i] <- similarity
  }
  diag(similarity_matrix) <- 1; 
  return(similarity_matrix)}


# Convert to lower triangle with indices (Option 3)
similarity_to_ltr_fun <- function(similarity_matrix){
  
lower_tri_indices <- which(lower.tri(similarity_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
  row = lower_tri_indices[, 1],
  col = lower_tri_indices[, 2],
  similarity = similarity_matrix[lower_tri_indices] )

 return(lower_tri_df)}


make_filename_safe <- function(name) {
  name <- iconv(name, from = "latin1", to = "UTF-8")
  name |>
    tolower() |>
    (\(x) gsub("[^a-z0-9]+", "_", x))() |>
    (\(x) gsub("^_|_$", "", x))() |>
    (\(x) substr(x, 1, 100))()
}

to_pylist <- function(vec) {
  if (length(vec) == 0) {
    list()
  } else {
    as.list(vec)
  }
}
