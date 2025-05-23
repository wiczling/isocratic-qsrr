### MCS 
extract_diff_fragments_combined <- function(idx1, idx2, smiles, .mcs_params = mcs_params, .fg_hierarchy = fg_hierarchy) {
  
  # Parse SMILES strings
  mol1 <- Chem$MolFromSmiles(smiles[idx1])
  mol2 <- Chem$MolFromSmiles(smiles[idx2])
  
  # Find Maximum Common Substructure (MCS)
  mcs_result <- rdFMCS$FindMCS(list(mol1, mol2), parameters = .mcs_params)
  smarts <- mcs_result$smartsString
  common_mol <- Chem$MolFromSmarts(smarts)
  
  # Get matching atoms in mol1
  match1 <- mol1$GetSubstructMatch(common_mol)
  
  # Set molecule and match atoms
  mol <- mol1
  match_atoms <- match1
  
  # Get all atoms and differing atoms
  all_atoms <- seq_len(mol$GetNumAtoms()) - 1L
  diff_atoms <- setdiff(all_atoms, match_atoms)
  
  # If no differing atoms, return empty data frame
  if (length(diff_atoms) == 0) {
    return(data.frame(
      atom_idx = integer(0),
      atom_symbol = character(0),
      matched_fg = character(0),
      frag_smiles = character(0),
      frag_id = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Initialize lists for results
  atom_list <- c()
  symbol_list <- c()
  fg_list <- c()
  smiles_list <- c()
  frag_id_list <- c()
  
  # --- Functional Group Matching ---
  atom_matches <- list()
  for (fg in .fg_hierarchy) {
    fg_name <- fg$name
    fg_smarts <- fg$smarts
    fg_mol <- Chem$MolFromSmarts(fg_smarts)
    
    if (is.null(fg_mol)) next
    
    # Get all matches of the functional group in the molecule
    matches <- mol$GetSubstructMatches(fg_mol)
    
    if (length(matches) == 0) next
    
    for (match in matches) {
      match <- as.integer(match)
      for (atom_idx in match) {
        if (atom_idx %in% diff_atoms) {
          key <- as.character(atom_idx)
          if (is.null(atom_matches[[key]])) {
            atom_matches[[key]] <- list(fg_name)
          } else {
            atom_matches[[key]] <- unique(c(atom_matches[[key]], fg_name))
          }
        }
      }
    }
  }
  
  # --- Fragment Extraction ---
  # Create a copy of the molecule
  rw_mol <- Chem$RWMol(mol)
  
  # Store original atom indices as properties
  for (i in all_atoms) {
    atom <- rw_mol$GetAtomWithIdx(i)
    atom$SetProp("orig_idx", as.character(i))
  }
  
  # Remove match_atoms (in reverse order to preserve indices)
  for (idx in sort(unlist(match_atoms), decreasing = TRUE)) {
    rw_mol$RemoveAtom(idx)
  }
  
  # Get fragments
  res <- try({
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE)
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE, sanitizeFrags = FALSE)
  } else {
    mol_frags <- rdmolops$GetMolFrags(rw_mol, asMols = TRUE)
  }
  
  # Process each fragment with a unique fragment ID
  frag_idx <- 0L
  for (frag in mol_frags) {
    frag_idx <- frag_idx + 1L
    frag_smile <- as.character(Chem$MolToSmiles(frag))
    frag_atoms <- c()
    
    # Get original atom indices from the fragment
    for (i in seq_len(frag$GetNumAtoms()) - 1L) {
      atom <- frag$GetAtomWithIdx(i)
      if (atom$HasProp("orig_idx")) {
        orig_idx <- as.integer(atom$GetProp("orig_idx"))
        if (orig_idx %in% diff_atoms) {
          frag_atoms <- c(frag_atoms, orig_idx)
        }
      }
    }
    
    # Add each atom's details to the lists
    if (length(frag_atoms) > 0) {
      for (atom_idx in frag_atoms) {
        atom_list <- c(atom_list, atom_idx)
        symbol_list <- c(symbol_list, mol$GetAtomWithIdx(atom_idx)$GetSymbol())
        
        # Get functional groups for this atom
        fg_names <- atom_matches[[as.character(atom_idx)]]
        fg_str <- if (is.null(fg_names)) NA_character_ else paste(fg_names, collapse = ", ")
        fg_list <- c(fg_list, fg_str)
        
        smiles_list <- c(smiles_list, frag_smile)
        frag_id_list <- c(frag_id_list, frag_idx)
      }
    }
  }
  
  # Create data frame
  result <- data.frame(
    atom_idx = atom_list,
    atom_symbol = symbol_list,
    matched_fg = fg_list,
    frag_smiles = smiles_list,
    frag_id = frag_id_list,
    stringsAsFactors = FALSE
  )
  
  # Sort by atom index for consistency
  if (nrow(result) > 0) {
    result <- result[order(result$atom_idx), ]
  }
  
  result$idx1 = idx1
  result$idx2 = idx2
  
  return(result)
}

extend_fg_hierarchy <- function(fg_hierarchy, new_groups) {
  for (group in new_groups) {
    if (!("name" %in% names(group)) || !("smarts" %in% names(group))) {
      warning("Invalid custom functional group structure.")
      next
    }
    fg_hierarchy <- append(fg_hierarchy, list(group))
  }
  return(fg_hierarchy)
}


FixBrokenAromaticity <- function(editable_mol) {
  # Create editable molecule
  ed_mol <- Chem$RWMol(editable_mol)

  # Try simple Kekulize
  success <- TRUE
  tryCatch({
    Chem$SanitizeMol(ed_mol)
    Chem$Kekulize(ed_mol, clearAromaticFlags=FALSE)
  }, error = function(e) {
    success <<- FALSE
  })
  
  if (!success) {
    # If Kekulize fails, fix problematic rings
    ri <- ed_mol$GetRingInfo()
    atom_rings <- ri$AtomRings()
    
    for (ring in atom_rings) {

      ring_py <- reticulate::r_to_py(ring)
      
      submol <- Chem$PathToSubmol(ed_mol, ring_py)
      
      ring_success <- TRUE
      tryCatch({
        Chem$Kekulize(submol, clearAromaticFlags=FALSE)
      }, error = function(e) {
        ring_success <<- FALSE
      })
      
      if (!ring_success) {
        # Clear aromaticity flags on broken ring
        for (idx in ring) {
          atom <- ed_mol$GetAtomWithIdx(as.integer(idx))
          atom$SetIsAromatic(FALSE)
        }
        if (length(ring) >= 2) {  # Only if the ring has 2+ atoms
          for (i in 1:(length(ring) - 1)) {
            for (j in (i + 1):length(ring)) {
              bond <- ed_mol$GetBondBetweenAtoms(as.integer(ring[i]), as.integer(ring[j]))
              if (!is.null(bond)) {
                bond$SetIsAromatic(FALSE)
              }
            }
          }
        }
      }
    }
    
    # Try Kekulize again after cleanup
    tryCatch({
      Chem$Kekulize(ed_mol, clearAromaticFlags=TRUE)
    }, error = function(e) {})
  }
  
  # Sanitize molecule
  Chem$SanitizeMol(ed_mol)
  
  return(ed_mol)
}

extract_common_fragments <- function(mol, match_atoms) {
  
  editable_mol <- Chem$RWMol(mol)
  mcs_atom_indices <- as.integer(match_atoms)
  
  # Remove all atoms not in MCS (do it in reverse to preserve indices)
  all_indices <- 0:(mol$GetNumAtoms() - 1)
  atoms_to_remove <- setdiff(all_indices, mcs_atom_indices)
  
  for (idx in rev(atoms_to_remove)) {
    editable_mol$RemoveAtom(as.integer(idx))
  }
  
  res <- try({
    common_mol <- FixBrokenAromaticity(editable_mol)
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    common_smiles <- "NA"
  } else {
    common_mol <- FixBrokenAromaticity(editable_mol)
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
      common1 = NA,
      common2 = NA,
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
  
  
   # gc()
   # py_run_string("import gc; gc.collect()")
   # img_magick <- image_read(img_file)
   # print(img_magick)  # Displays in RStudio viewer or default graphics device
    # 
  # Return data frame row
  data.frame(
    smile1 = smiles1,
    smile2 = smiles2,
    common1 = common_smile1,
    common2 = common_smile2,
    to_remove = to_remove,
    to_remove_str= to_remove_str,
    to_add = to_add,
    to_add_str=to_add_str,
    structure = sprintf('<img src="%s" style="width:600px; height:auto;"/>', img_base64), 
    stringsAsFactors = FALSE)
}

# Function to get number of atoms from SMILES
get_atom_count <- function(smiles) {
  mol <- Chem$MolFromSmiles(smiles)
  return(mol$GetNumAtoms())
}

are_smiles_identical <- function(smiles1, smiles2) {
  
  mol1 <- Chem$MolFromSmiles(smiles1)
  mol2 <- Chem$MolFromSmiles(smiles2)
  
  if (is.null(mol1) || is.null(mol2)) {
    return(NA)
  }

  can1 <- Chem$MolToSmiles(mol1, canonical = TRUE)
  can2 <- Chem$MolToSmiles(mol2, canonical = TRUE)
  return(identical(can1, can2))

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

make_distance_matrix_fun <- function(lower_tri_df) {
  n <- max(lower_tri_df$row, lower_tri_df$col)
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  for (k in 1:nrow(lower_tri_df)) {
    i <- lower_tri_df$row[k]
    j <- lower_tri_df$col[k]
    distance <- lower_tri_df$distance[k]
    distance_matrix[i, j] <- distance
    distance_matrix[j, i] <- distance
  }
  diag(distance_matrix) <- 0; 
  return(distance_matrix)}

# Convert to lower triangle with indices
similarity_to_ltr_fun <- function(similarity_matrix){
  
lower_tri_indices <- which(lower.tri(similarity_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
  row = lower_tri_indices[, 1],
  col = lower_tri_indices[, 2],
  similarity = similarity_matrix[lower_tri_indices] )

 return(lower_tri_df)}

# Convert to lower triangle with indices
distance_to_ltr_fun <- function(distance_matrix){
  
  lower_tri_indices <- which(lower.tri(distance_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
    row = lower_tri_indices[, 1],
    col = lower_tri_indices[, 2],
    distance = distance_matrix[lower_tri_indices] )
  
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
