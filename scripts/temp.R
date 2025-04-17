custom_fg  <- list(
  "Allenic Carbon" = "[CX2](=C)=C",
  "Vinylic Carbon" = "[CX3]=[CX3]",
  # Alkynes
  "Acetylenic Carbon" = "[CX2]#C",
  # Carbonyl-related groups
  "Carbonyl (Low Specificity)" = "[CX3]=[OX1]",
  "Carbonyl (Both Resonance Forms)" = "[$([CX3]=[OX1]),$([CX3+]-[OX1-])]",
  "Carbonyl with Carbon" = "[CX3](=[OX1])C",
  "Carbonyl with Nitrogen" = "[OX1]=CN",
  "Carbonyl with Oxygen" = "[CX3](=[OX1])O",
  "Acyl Halide" = "[CX3](=[OX1])[F,Cl,Br,I]",
  "Aldehyde" = "[CX3H1](=O)[#6]",
  "Anhydride" = "[CX3](=[OX1])[OX2][CX3](=[OX1])",
  "Amide" = "[NX3][CX3](=[OX1])[#6]",
  "Amidinium" = "[NX3][CX3]=[NX3+]",
  "Carbamate" = "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
  "Carbamic Ester" = "[NX3][CX3](=[OX1])[OX2H0]",
  "Carbamic Acid" = "[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]",
  "Carboxylate Ion" = "[CX3](=O)[O-]",
  "Carbonic Acid or Ester" = "[CX3](=[OX1])(O)O",
  "Carbonic Acid or Monoester" = "[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]",
  "Carbonic Diester" = "C[OX2][CX3](=[OX1])[OX2]C",
  "Carboxylic Acid" = "[CX3](=O)[OX2H1]",
  "Carboxylic Acid or Conjugate Base" = "[CX3](=O)[OX1H0-,OX2H1]",
  "Cyanamide" = "[NX3][CX2]#[NX1]",
  "Ester" = "[#6][CX3](=O)[OX2H0][#6]",
  "Ketone" = "[#6][CX3](=O)[#6]",
  
  # Ethers
  "Ether" = "[OD2]([#6])[#6]",
  
  # Nitrogen-containing groups
  "Primary or Secondary Amine" = "[NX3;H2,H1;!$(NC=O)]",
  "Enamine" = "[NX3][CX3]=[CX3]",
  "Primary Amine (Not Amide)" = "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",
  "Two Primary or Secondary Amines" = "[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]",
  "Enamine or Aniline Nitrogen" = "[NX3][$(C=C),$(cc)]",
  
  # Azides
  "Azide Group" = "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]",
  "Azide Ion" = "[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]",
  
  # Azo compounds
  "Azo Nitrogen (Low Specificity)" = "[NX2]=N",
  "Diazene" = "[NX2]=[NX2]",
  "Azoxy Nitrogen" = "[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]",
  "Diazo Nitrogen" = "[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]",
  "Azole" = "[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]",
  
  # Hydrazines and hydrazones
  "Hydrazine" = "[NX3][NX3]",
  "Hydrazone" = "[NX3][NX2]=[*]",
  
  # Imines
  "Substituted Imine" = "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]",
  "Imine (Schiff Base)" = "[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]",
  "Iminium" = "[NX3+]=[CX3]",
  
  # Imides
  "Unsubstituted Dicarboximide" = "[CX3](=[OX1])[NX3H][CX3](=[OX1])",
  "Substituted Dicarboximide" = "[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])",
  "Dicarboxdiimide" = "[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])",
  
  # Nitrates
  "Nitrate Group" = "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
  "Nitrate Anion" = "[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]",
  
  # Nitriles
  "Nitrile" = "[NX1]#[CX2]",
  "Isonitrile" = "[CX1-]#[NX2+]",
  
  # Nitro compounds
  "Nitro Group" = "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
  "Two Nitro Groups" = "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8].[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
  
  # Nitroso
  "Nitroso Group" = "[NX2]=[OX1]",
  
  # N-Oxides
  "N-Oxide" = "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
  
  # Hydroxyl groups
  "Hydroxyl" = "[OX2H]",
  "Hydroxyl in Alcohol" = "[#6][OX2H]",
  "Hydroxyl in Carboxylic Acid" = "[OX2H][CX3]=[OX1]",
  "Hydroxyl in H-O-P-" = "[OX2H]P",
  "Enol" = "[OX2H][#6X3]=[#6]",
  "Phenol" = "[OX2H][cX3]:[c]",
  "Enol or Phenol" = "[OX2H][$(C=C),$(cc)]",
  "Acidic Hydroxyl" = "[$([OH]-*=[!#6])]",
  
  # Peroxides
  "Peroxide" = "[OX2,OX1-][OX2,OX1-]",
  
  # Phosphoric compounds
  "Phosphoric Acid Groups" = "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
  "Phosphoric Ester Groups" = "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
  
  # Sulfur-containing groups
  "Carbo-Thiocarboxylate" = "[S-][CX3](=S)[#6]",
  "Carbo-Thioester" = "S([#6])[CX3](=O)[#6]",
  "Thio Analog of Carbonyl" = "[#6X3](=[SX1])([!N])[!N]",
  "Thiol, Sulfide, or Disulfide Sulfur" = "[SX2]",
  "Thiol" = "[#16X2H]",
  "Sulfur with At Least One Hydrogen" = "[#16!H0]",
  
  # Sulfides
  "Sulfide" = "[#16X2H0]",
  "Mono-Sulfide" = "[#16X2H0][!#16]",
  "Di-Sulfide" = "[#16X2H0][#16X2H0]",
  "Two Sulfides" = "[#16X2H0][!#16].[#16X2H0][!#16]",
  
  # Sulfinates
  "Sulfinate" = "[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
  "Sulfinic Acid" = "[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
  
  # Sulfones
  "Sulfone (Low Specificity)" = "[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]",
  "Sulfone (High Specificity)" = "[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]",
  "Sulfonic Acid" = "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
  "Sulfonate" = "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
  "Sulfonamide" = "[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]",
  "Carbo-Azosulfone" = "[SX4](C)(C)(=O)=N",
  
  # Sulfoxides
  "Sulfoxide (Low Specificity)" = "[$([#16X3]=[OX1]),$([#16X3+][OX1-])]",
  "Sulfoxide (High Specificity)" = "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
  
  # Sulfates
  "Sulfate" = "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
  "Sulfuric Acid Ester (Low Specificity)" = "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]",
  "Sulfuric Acid Diester" = "[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]",
  
  # Sulfamates
  "Sulfamate" = "[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]",
  "Sulfamic Acid" = "[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]",
  
  # Sulfenes
  "Sulfenic Acid" = "[#16X2][OX2H,OX1H0-]",
  "Sulfenate" = "[#16X2][OX2H0]",
  
  # Halides
  "Any Carbon with Halogen" = "[#6][F,Cl,Br,I]",
  "Halogen" = "[F,Cl,Br,I]",
  "Three Halides" = "[F,Cl,Br,I].[F,Cl,Br,I].[F,Cl,Br,I]"
)