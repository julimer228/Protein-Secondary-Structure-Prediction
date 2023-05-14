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




