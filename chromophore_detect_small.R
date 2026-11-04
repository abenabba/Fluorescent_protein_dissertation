# --------------------------------------------------
# 1. Enter amino acid sequences
# --------------------------------------------------

proteins <- c(
  "2099C" = paste0(
    "MCSCGGHYRC",
    "DFKTIYKAKK",
    "VVQLPDYHSI",
    "DHRIEILSHD",
    "KDYNKVELHE",
    "HASSRHGSRH",
    "GMDELYKALR",
    "VAASLLNIKQ",
    "WMVMRLSVV"
  ),
  
  "5967S" = paste0(
    "MSASAATDTN",
    "KSTKGMDELY",
    "KGSQELADYI",
    "QERTVSKGPV",
    "LIPFTHYLSF",
    "QTAISKDRNE",
    "TRDHMVFLEF",
    "FSACGHTHGM",
    "DELYHCGQAS",
    "RLTFGSDCLV",
    "YKIRFDGVNF",
    "PPKSPVIQKT",
    "TMKWEPSTEI",
    "LYVRDGVAGW",
    "RH"
  )
)

# Verify lengths, including initiating methionine
print(nchar(proteins))

stopifnot(
  nchar(proteins["2099C"]) == 89,
  nchar(proteins["5967S"]) == 142
)


# --------------------------------------------------
# 2. Scan for potential chromophore tripeptides
# --------------------------------------------------

find_chromophore_motifs <- function(sequence, protein_name) {
  
  # Examine every overlapping 3-residue window
  positions <- seq_len(nchar(sequence) - 2)
  
  triplets <- substring(
    sequence,
    positions,
    positions + 2
  )
  
  # Chromophore-like motifs:
  # X = any amino acid
  # Middle residue = Y, H, W, or F
  # Third residue = G
  
  middle <- substr(triplets, 2, 2)
  third  <- substr(triplets, 3, 3)
  
  matches <- middle %in% c("Y", "H", "W", "F") &
    third == "G"
  
  matched_positions <- positions[matches]
  
  # Extract sequence context: 5 residues on either side
  context <- substring(
    sequence,
    pmax(1, matched_positions - 5),
    pmin(nchar(sequence), matched_positions + 7)
  )
  
  data.frame(
    protein = protein_name,
    start = matched_positions,
    end = matched_positions + 2,
    motif = triplets[matches],
    category = ifelse(
      middle[matches] %in% c("Y", "H"),
      "Y/H candidate",
      "W/F exploratory"
    ),
    context = context
  )
}


# --------------------------------------------------
# 3. Analyze both proteins
# --------------------------------------------------

results <- do.call(
  rbind,
  Map(
    find_chromophore_motifs,
    unname(proteins),
    names(proteins)
  )
)

rownames(results) <- NULL

print(results)