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


#--------Sliding window encoding-------------------------------------------------------------------------------------  

library(caret)  # Load the caret package for model training and evaluation

# Step 3: Grid search for optimal window length

window_lengths <- seq(5, 25, by = 5)  # Set up a sequence of window lengths to test
cv_folds <- 5  # Number of cross-validation folds

# Function to create sliding window encoded data for a given window size
create_encoded_data <- function(window_size, protein_data) {
  encoded_data <- list()
  
  for (id in names(protein_data)) {
    protein <- protein_data[[id]]
    sequence <- protein$sequence
    encoded_sequence <- protein$encoded_sequence
    
    for (i in 1:(nchar(sequence) - window_size + 1)) {
      window <- substring(sequence, i, i + window_size - 1)
      encoded_window <- substring(encoded_sequence, i, i + window_size - 1)
      
      encoded_data[[length(encoded_data) + 1]] <- list(window = window, encoded_window = encoded_window, label = "H/E/C")
    }
  }
  
  return(encoded_data)
}

# Function to train and evaluate a model with a specific window length using cross-validation
evaluate_window_length <- function(window_length, protein_data, cv_folds) {
  # Create encoded data for the given window length
  encoded_data <- create_encoded_data(window_length, protein_data)
  
  # Convert the encoded data into a data frame for use with caret
  encoded_data_df <- do.call(rbind, lapply(encoded_data, data.frame))
  
  # Set up cross-validation
  control <- trainControl(method = "cv", number = cv_folds, classProbs = TRUE)
  
  # Train and evaluate the model
  model <- train(label ~ ., data = encoded_data_df, method = "rf", trControl = control)
  accuracy <- max(model$results$Accuracy)
  
  return(accuracy)
}

# Perform grid search to find the optimal window length
results <- data.frame()
for (window_length in window_lengths) {
  accuracy <- evaluate_window_length(window_length, protein_data, cv_folds)
  results <- rbind(results, data.frame(window_length = window_length, accuracy = accuracy))
}

# Print the results
print(results)

# Find the optimal window length
optimal_window_length <- results[which.max(results$accuracy), "window_length"]
print(paste("Optimal window length:", optimal_window_length))

#--------End of Sliding window encoding-------------------------------------------------------------------------------------  

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




