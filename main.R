library(stringr)  # Load the stringr package

# Step 1: Preparing the data

# Specify the path to the .txt file
file_path <- "dane/train.txt"

# Read the contents of the file
lines <- readLines(file_path)

# Initialize variables
protein_data <- list()
current_protein <- NULL

# Process the lines
for (i in 1:length(lines)) {
  line <- lines[i]

  if (startsWith(line, ">")) {
    # Start of a new protein
    if (!is.null(current_protein)) {
      # Store the previous protein
      protein_data[[current_protein$id]] <- current_protein
    }

    # Create a new protein object
    current_protein <- list()
    current_protein$id <- substring(line, 2)  # Remove the ">" symbol
    current_protein$sequence <- ""
  } else if (i %% 3 == 2) {
    # Store the protein sequence
    current_protein$sequence <- line
  } else if (i %% 3 == 0) {
    # Store the orthogonally encoded sequence
    current_protein$encoded_sequence <- line
  }
}

# Store the last protein
protein_data[[current_protein$id]] <- current_protein

# Print the stored protein data
for (id in names(protein_data)) {
  protein <- protein_data[[id]]
  print(paste("ID:", protein$id))
  print(paste("Sequence:", protein$sequence))
  print(paste("Encoded Sequence:", protein$encoded_sequence))
}

#--------Part with orthogonal encoding-------------------------------------------------------------------------------------

# Define a function to orthogonally encode protein sequences
orthogonal_encode <- function(sequence) {
  # Define a list of amino acids
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  # Initialize an empty matrix to store the encoded sequence
  encoded_sequence <- matrix(0, nrow = length(sequence), ncol = length(amino_acids))
  
  # Loop through the sequence and encode each amino acid
  for (i in 1:length(sequence)) {
    aa <- substring(sequence, i, i)
    aa_index <- match(aa, amino_acids)
    if (!is.na(aa_index)) {
      encoded_sequence[i, aa_index] <- 1
    }
  }
  
  return(encoded_sequence)
}

# Add the orthogonal encoding to each protein
for (id in names(protein_data)) {
  protein <- protein_data[[id]]
  protein$encoded_sequence <- orthogonal_encode(protein$sequence)
}

# Print the stored protein data with orthogonal encoding
for (id in names(protein_data)) {
  protein <- protein_data[[id]]
  print(paste("ID:", protein$id))
  print(paste("Sequence:", protein$sequence))
  print("Encoded Sequence:")
  print(protein$encoded_sequence)
}
#--------End of Part with orthogonal encoding-------------------------------------------------------------------------------------

  

# Step 2: Sliding window encoding

# window_size <- 20  # Specify the desired window size
#
# # Initialize variables
# encoded_data <- list()
#
# # Process each protein entry
# for (id in names(protein_data)) {
#   protein <- protein_data[[id]]
#   sequence <- protein$sequence
#   encoded_sequence <- protein$encoded_sequence
#
#   # Apply sliding window encoding
#   for (i in 1:(nchar(sequence) - window_size + 1)) {
#     window <- substring(sequence, i, i + window_size - 1)
#     encoded_window <- substring(encoded_sequence, i, i + window_size - 1)
#
#     # Store the encoded window with its corresponding secondary structure class label
#     encoded_data[[length(encoded_data) + 1]] <- list(window = window, encoded_window = encoded_window, label = "H/E/C")
#   }
# }




