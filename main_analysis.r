################################################################################
# PROJECT: PCA by Efficient Power Method Implementation on Sparse Movie Ratings Data
# AUTHOR: ROHAN WAGLE
################################################################################

##### Setting Up Libraries #####
if (!require("Matrix")) install.packages("Matrix")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("plotly")) install.packages("plotly")

library(Matrix); library(ggplot2); library(ggrepel); 
library(stringr); library(dplyr); library(gridExtra); library(plotly)


####### LOADING AND STORING DATA #####

# Load Ratings (u.data)
ratings <- read.table("ml-100k/u.data", sep = "\t", 
                      col.names = c("user_id", "movie_id", "rating", "timestamp"))

# Load Users (u.user)
users <- read.table("ml-100k/u.user", sep = "|", 
                    col.names = c("user_id", "age", "gender", "occupation", "zip"))

# Load Movie Metadata (u.item)
movies_raw <- read.table("ml-100k/u.item", sep = "|", quote = "", 
                         encoding = "Latin-1", header = FALSE)
movie_titles <- movies_raw[, 1:2]
colnames(movie_titles) <- c("movie_id", "title")

# Mapping Genres to Ratings Data (Splitting strings like "Action|1")
genre_meta <- read.table("ml-100k/u.genre")
genre_names <- sapply(strsplit(as.character(genre_meta$V1), "\\|"), `[`, 1)
colnames(movies_raw)[6:24] <- genre_names

primary_genre <- apply(movies_raw[, 6:24], 1, function(x) {
  match_idx <- which(x == 1)
  if (length(match_idx) > 0) return(genre_names[match_idx[1]])
  return("None")
})

####### Dataset Features and Demographics ###########

### Creating Freqency Distribution of Movies by Genre

# Sum the binary genre wise columns 6-24 from raw data (mapped using u.genre)
genre_counts <- colSums(movies_raw[, 6:24])

# Create a clean data frame for plotting
genre_dist <- data.frame(
  Genre = names(genre_counts),
  MovieCount = as.numeric(genre_counts)
)

# Sort the data frame so the plot is ordered from most to least frequent
genre_dist <- genre_dist[order(-genre_dist$MovieCount), ]

# Generate Plot
ggplot(genre_dist, aes(x = reorder(Genre, MovieCount), y = MovieCount)) +
  geom_bar(stat = "identity", fill = "midnightblue", alpha = 0.8) +
  geom_text(aes(label = MovieCount), hjust = -0.2, size = 3, color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "MovieLens 100K Dataset: Number of Movies per Genre",
       x = "Genre",
       y = "Total Number of Movies") +
  # Adjusting y-axis limit to accommodate the labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

########## Power Method Implementation ##########

# Prepare Matrix
X <- sparseMatrix(i = ratings$user_id, j = ratings$movie_id, x = ratings$rating)
n <- nrow(X); p <- ncol(X)
column_means <- as.numeric(colMeans(X)) 

# Analysing Data Matrix Sparsity

# Calculate Sparsity Metrics
total_possible <- nrow(X) * ncol(X)
actual_ratings <- nnzero(X)  # Number of non-zero entries in sparse matrix
missing_ratings <- total_possible - actual_ratings
sparsity_pct <- (1 - (actual_ratings / total_possible)) * 100
density_pct <- (actual_ratings / total_possible) * 100

# Create Summary Table
sparsity_table <- data.frame(
  Metric = c("Total Number of Cells", "Actual Ratings", 
             "Missing Ratings", "Sparsity Percentage", "Density Percentage"),
  Value = c(format(total_possible, big.mark=","), 
            format(actual_ratings, big.mark=","), 
            format(missing_ratings, big.mark=","),
            paste0(round(sparsity_pct, 2), "%"),
            paste0(round(density_pct, 2), "%"))
)

print("--- Data Matrix Sparsity Table ---")
print(sparsity_table)

# Create Matrix Spy Plot to Visualise Sparsity

image(X, 
      main = "Data Matrix Spy Plot: MovieLens 100K", sub = "Black dots indicate observed ratings; White space indicates missing data",
      xlab = "Movie ID (1-1682)", ylab = "User ID (1-943)",
      colorkey = FALSE, col.regions = "black")


# Function for Gram Schmidt Algorithm
standard_gram_schmidt <- function(A) {
  k_dim <- ncol(A); p_dim <- nrow(A)
  Q <- matrix(0, nrow = p_dim, ncol = k_dim)
  for (j in 1:k_dim) {
    v_j <- A[, j]
    if (j > 1) {
      for (i in 1:(j - 1)) {
        v_j <- v_j - sum(A[, j] * Q[, i]) * Q[, i]
      }
    }
    norm_v_j <- sqrt(sum(v_j^2))
    Q[, j] <- if (norm_v_j < 1e-15) 0 else v_j / norm_v_j
  }
  return(Q)
}

# Power Method Loop
k <- 10; epsilon <- 1e-6; max_iter <- 1500; set.seed(42)

V_raw_previous <- matrix(0, nrow = p, ncol = k)
Q_current <- standard_gram_schmidt(matrix(rnorm(p * k), p, k))
convergence_log <- c()

cat("Executing Power Method...\n")

for (m in 1:max_iter) {
  # Forward Pass
  Z <- X %*% Q_current
  mu_T_Q <- as.numeric(t(column_means) %*% Q_current)
  U_tmp <- sweep(Z, 2, mu_T_Q, "-")
  
  # Backward Pass
  W <- t(X) %*% U_tmp
  V_raw_current <- as.matrix(W - (column_means %o% colSums(U_tmp)))
  
  # Convergence
  dist_sum <- sum(sqrt(colSums((V_raw_current - V_raw_previous)^2)))
  convergence_log <- c(convergence_log, dist_sum)
  if (dist_sum < epsilon) {
    cat("Converged at iteration:", m, "\n"); break
  }
  
  V_raw_previous <- V_raw_current
  Q_current <- standard_gram_schmidt(V_raw_current)
}

####### Recovering Outputs ####

# Recover Squared Singular Values 
eigenvalues <- sqrt(colSums(V_raw_current^2))
sigma <- sqrt(eigenvalues)

# Total Market Variance (Squared Frobenius Norm of Centered X)
total_var <- sum(X^2) - n * sum(column_means^2)

# Calculate Proportion of Variance Captured 
prop_variance <- eigenvalues / total_var

# Extract U and V Orthogonal Matrices
V_matrix <- Q_current 
U_matrix <- sweep(U_tmp, 2, sigma, "/")

cat(sprintf("\nTotal Variance Explained by k=10: %.4f (%.2f%%)\n", 
            sum(prop_variance), sum(prop_variance) * 100))

######### Convergence Related Plots and Metrics #####

# Final 10 Steps
n_iters <- length(convergence_log)
last_10_idx <- (n_iters - 9):n_iters
last_10_idx <- last_10_idx[last_10_idx > 0]

# Calculate the convergence ratio
reduction_rates <- c(NA, convergence_log[2:n_iters] / convergence_log[1:(n_iters-1)])

convergence_summary <- data.frame(
  Iteration = last_10_idx,
  Error_Magnitude = format(convergence_log[last_10_idx], scientific = TRUE, digits = 4),
  Convergence_Ratio = round(reduction_rates[last_10_idx], 5)
)

cat("\n--- Final Subspace Convergence: Last 10 Iterations ---\n")
print(convergence_summary)

# Log-linear Convergence Plot
conv_df <- data.frame(
  Iteration = 1:n_iters, 
  Error = convergence_log
)

ggplot(conv_df, aes(x = Iteration, y = Error)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue", size = 0.7, shape = 16) +
  
  geom_hline(yintercept = epsilon, linetype = "dotted", color = "red", size = 0.8) +
  
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  coord_cartesian(ylim = c(1e-8, max(convergence_log) * 1.2), expand = FALSE) +
  
  theme_bw() + 
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(size = 11)
  ) +
  labs(
    title = expression(paste("Convergence of power method for ", lambda, " = 1, 2, ..., 10")),
    x = "number of iterations",
    y = "error in eigenvector"
  ) +
  annotate("text", x = n_iters * 0.15, y = epsilon * 0.4, 
           label = "Error Threshold", color = "red", size = 3, fontface = "italic")

# Verifying Orthogonalisation

ortho_check <- t(V_matrix) %*% V_matrix
ortho_summary <- round(as.matrix(ortho_check), 10)

cat("\n--- Orthogonality Check (V^T * V) ---\n")
print(ortho_summary)

# Scree Plot

# Prepare Data
scree_df <- data.frame(
  PC = 1:length(prop_variance),
  Variance = prop_variance
)

# Plot
ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_line(color = "blue", size = 0.6) +
  geom_point(color = "blue", size = 1.5, shape = 16) +
  
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01),
    limits = c(0, max(prop_variance) * 1.1),
    expand = c(0, 0) 
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    axis.title = element_text(size = 11)
  ) +
  
  labs(
    title = "Scree Plot",
    x = "Principal Component Index",
    y = "Proportion of Variance Explained"
  )

####### Analysis Of Loadings #####

# Prepare Concise Genre Labelling
top_genres <- c("Drama", "Comedy", "Thriller", "Action", "Romance")

# Clean titles: Remove (Year) AND shorten Dr. Strangelove
clean_titles <- movie_titles$title %>%
  str_remove("\\s\\(\\d{4}\\)$") %>%
  str_replace("Dr. Strangelove or: How I Learned to Stop Worrying and Love the Bomb", "Dr. Strangelove")

movie_analysis <- data.frame(
  title = clean_titles, 
  PC1 = V_matrix[, 1], 
  PC2 = V_matrix[, 2], 
  OriginalGenre = as.character(primary_genre)
)

# Genre Aggregation
movie_analysis$PlotGenre <- ifelse(
  movie_analysis$OriginalGenre %in% top_genres, 
  movie_analysis$OriginalGenre, 
  "Others"
)

movie_analysis$PlotGenre <- factor(
  movie_analysis$PlotGenre, 
  levels = c(top_genres, "Others")
)

# Selecting Movies With Extreme Loadings for Labelling
movie_analysis$magnitude <- sqrt(movie_analysis$PC1^2 + movie_analysis$PC2^2)

extreme_indices <- unique(c(
  order(movie_analysis$PC1)[1:6],
  order(movie_analysis$PC1, decreasing=T)[1:6],
  order(movie_analysis$PC2)[1:6],
  order(movie_analysis$PC2, decreasing=T)[1:6],
  order(movie_analysis$magnitude, decreasing=T)[1:15]
))

top_labels <- movie_analysis[extreme_indices, ]

# Generate loadings Plot
ggplot(movie_analysis, aes(x = PC1, y = PC2, color = PlotGenre)) +
  geom_point(aes(alpha = ifelse(PlotGenre == "Others", 0.4, 0.7),
                 size = ifelse(PlotGenre == "Others", 1.4, 1.2))) +
  
  geom_text_repel(
    data = top_labels, 
    aes(label = title), 
    size = 3.2, 
    color = "black", 
    fontface = "bold",
    box.padding = 1.0,         # Space around text
    point.padding = 0.5,       # Space from the dot
    force = 40,                
    force_pull = 0.15,         # Allow labels to float away from crowded spots
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.4,
    min.segment.length = 0,
    seed = 42,
    direction = "both"
  ) +
  
  # Colour Palette
  scale_color_manual(values = c(
    "Drama" = "#E41A1C", "Comedy" = "#377EB8", "Thriller" = "#4DAF4A", 
    "Action" = "#984EA3", "Romance" = "#FF7F00", "Others" = "#2F4F4F"
  )) +
  
  scale_alpha_identity() +
  scale_size_identity() +
  
  # Expand canvas 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  
  # Visual Guides
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.4) +
  
  theme_minimal() +
  labs(
    title = "PCA Loadings Plot: MovieLens 100K Data",
    x = "PC1",
    y = "PC2",
    color = "Genre"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank()
  )

# Generate interactive loadings plots

# Extract Release Year from Titles using Regex
# Standard format is "Movie Title (1995)"
movie_analysis$Year <- str_extract(movie_titles$title, "\\(\\d{4}\\)") %>%
  str_replace_all("[\\(\\)]", "") # Remove parentheses

# Construct the Full Genre Profile
# This uses the binary columns (6-24) from the raw data
movie_analysis$FullProfile <- apply(movies_raw[, 6:24], 1, function(x) {
  match_idxs <- which(x == 1)
  if (length(match_idxs) > 0) {
    return(paste(genre_names[match_idxs], collapse = " | "))
  }
  return("None")
})

# 3. Build the plot
# Note: the 'text' aesthetic is what plotly uses for the hover tool
p <- ggplot(movie_analysis, aes(x = PC1, y = PC2, color = PlotGenre,
                                text = paste0("<b>", title, "</b><br>",
                                              "Year: ", Year, "<br>",
                                              "Genres: ", FullProfile))) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = market_colors) + 
  theme_minimal() +
  labs(title = "Interactive PC Loadings Plot: MovieLens 100K",
       x = "PC1 (Market Engagement)",
       y = "PC2 (Critically Acclaimed vs. Blockbusters)")

# Convert to Plotly
# 'tooltip' argument tells plotly to only show the 'text' aesthetic we built
interactive_plot <- ggplotly(p, tooltip = "text")

# 5. Save as HTML (This file can be hosted or shared)
htmlwidgets::saveWidget(interactive_plot, "PCA_Interactive_Loadings.html")

# Display in RStudio Viewer
interactive_plot

# Compact Loadings Plots (PC1-PC4)

# Consolidate Title and Genre Data
movie_titles_clean <- clean_titles 
genre_mapping <- as.character(primary_genre)

# Re-define Colour Palette to ensure the plot has its palette
market_colors <- c(
  "Drama" = "#E41A1C", "Comedy" = "#377EB8", "Thriller" = "#4DAF4A", 
  "Action" = "#984EA3", "Romance" = "#FF7F00", "Others" = "#2F4F4F"
)

combined_pc_data <- lapply(1:4, function(k_idx) {
  # Get indices for the most positive and most negative loadings
  top_5 <- order(V_matrix[, k_idx], decreasing = TRUE)[1:5]
  bot_5 <- order(V_matrix[, k_idx], decreasing = FALSE)[1:5]
  
  data.frame(
    PC = paste0("PC ", k_idx),
    title = movie_titles_clean[c(bot_5, top_5)],
    genre = genre_mapping[c(bot_5, top_5)],
    score = V_matrix[c(bot_5, top_5), k_idx]
  )
}) %>% bind_rows()

# Create Plot
ggplot(combined_pc_data, aes(x = reorder(title, score), y = score, fill = genre)) +
  geom_bar(stat = "identity", width = 0.75) +
  
  # Logic to place title labels: inside the bar, adjusted by sign
  geom_text(aes(label = title, 
                y = ifelse(score > 0, -0.0005, 0.0005), 
                hjust = ifelse(score > 0, 1.05, -0.05)), 
            size = 2.5, fontface = "bold") +
  
  scale_fill_manual(values = market_colors) +
  scale_x_discrete(expand = expansion(add = c(1, 1))) +
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  
  # Faceting by PC to show all 4 PC plots at once
  facet_wrap(~PC, scales = "free_y", ncol = 2) +
  
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "bottom"
  ) +
  labs(
    title = "Top/Bottom Movie Loadings (PC 1-4)",
    x = NULL, y = "Loadings",
    fill = "Genre"
  )

########## USER PCA SCORES ANALYSIS ##########

# Scale U Matrix
U_scaled <- sweep(U_matrix, 2, sigma, "*")

# Convert to Data Frame and Merge
user_scores <- data.frame(user_id = 1:nrow(U_scaled), as.matrix(U_scaled))
colnames(user_scores)[2:11] <- paste0("PC", 1:10)

analysis_df <- merge(user_scores, users, by = "user_id")

# T-tests for Gender on the PC1 and PC2 axis
t.test(PC1 ~ gender, data = analysis_df)
t.test(PC2 ~ gender, data = analysis_df)

# Filter scores by gender
female_scores <- analysis_df$PC1[analysis_df$gender == "F"]
male_scores   <- analysis_df$PC1[analysis_df$gender == "M"]

# Conduct one-sided T-test (Testing if Male mean is GREATER)
pc1_onesided <- t.test(male_scores, female_scores, 
                       alternative = "greater", 
                       var.equal = FALSE) # Welch's T-test (safer for unequal N)

print(pc1_onesided)


# Calculate Centroids of scores along PC1 and PC2 by Gender
gender_centroids <- analysis_df %>%
  group_by(gender) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# Create PC Scores Plot
ggplot(analysis_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = gender), alpha = 0.5, size = 1.5) +
  
  geom_point(data = gender_centroids, aes(fill = gender), 
             size = 2, shape = 23, color = "black", stroke = 1.5) +
  
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  
  scale_color_manual(values = c("M" = "#377EB8", "F" = "#E41A1C")) +
  scale_fill_manual(values = c("M" = "#377EB8", "F" = "#E41A1C")) +
  
  theme_minimal() +
  labs(
    title = "PC Scores Plot: MovieLens 100K Data",
    x = "PC1",
    y = "PC2"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Age-related Analysis

# Correlation Tests
cor_pc1 <- cor.test(analysis_df$age, analysis_df$PC1)
cor_pc2 <- cor.test(analysis_df$age, analysis_df$PC2)

cat("Correlation Age vs PC1 Scores: ", cor_pc1$estimate, " (p =", cor_pc1$p.value, ")\n")
cat("Correlation Age vs PC2 Scores: ", cor_pc2$estimate, " (p =", cor_pc2$p.value, ")\n")

# Visualizations
plot1 <- ggplot(analysis_df, aes(x = age, y = PC1)) +
  geom_point(alpha = 0.5, color = "#2980b9") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_minimal() +
  labs(title = "Age vs. PC Score (PC1)",
       subtitle = paste("Correlation:", round(cor_pc1$estimate, 3)),
       x = "Age", y = "PC1 Score")

plot2 <- ggplot(analysis_df, aes(x = age, y = PC2)) +
  geom_point(alpha = 0.5, color = "#c0392b") +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_minimal() +
  labs(title = "Age vs. PC Score (PC2)",
       subtitle = paste("Correlation:", round(cor_pc2$estimate, 3)),
       x = "Age", y = "PC2 Score")

grid.arrange(plot1, plot2, ncol = 2)

################################################################################
# ADDITIONAL ANALYSIS: Gender-Specific PCA using Simultaneous Power Method
################################################################################

# Subset the existing full data matrix 'X' using gender IDs from 'users'
female_user_ids <- users$user_id[users$gender == "F"]
male_user_ids   <- users$user_id[users$gender == "M"]

X_female <- X[female_user_ids, , drop = FALSE]
X_male   <- X[male_user_ids, , drop = FALSE]

# Simultaneous Power Method Function
run_gender_pca <- function(X_sub, k = 10, epsilon = 1e-6, max_iter = 1000) {
  n_sub <- nrow(X_sub)
  p_sub <- ncol(X_sub)
  col_means_sub <- as.numeric(colMeans(X_sub))
  
  set.seed(42)
  Q_curr <- as.matrix(standard_gram_schmidt(matrix(rnorm(p_sub * k), p_sub, k)))
  V_raw_prev <- matrix(0, nrow = p_sub, ncol = k)
  conv_log <- c()
  
  for (m in 1:max_iter) {
    # Forward Pass
    Z_sub <- X_sub %*% Q_curr
    mu_T_Q_sub <- as.numeric(t(col_means_sub) %*% Q_curr)
    U_tmp_sub <- sweep(as.matrix(Z_sub), 2, mu_T_Q_sub, "-")
    
    # Backward Pass
    W_sub <- t(X_sub) %*% U_tmp_sub
    V_raw_curr <- as.matrix(W_sub - (col_means_sub %o% colSums(U_tmp_sub)))
    
    # Convergence Tracking
    dist <- sum(sqrt(colSums((V_raw_curr - V_raw_prev)^2)))
    conv_log <- c(conv_log, dist)
    if (dist < epsilon) break
    
    V_raw_prev <- V_raw_curr
    Q_curr <- standard_gram_schmidt(V_raw_curr)
  }
  
  vals <- sqrt(colSums(V_raw_curr^2))
  tot_v <- sum(X_sub^2) - n_sub * sum(col_means_sub^2)
  
  return(list(V = Q_curr, eigenvalues = vals, log = conv_log, total_var = tot_v))
}

# --- 3. Execution ---
cat("Executing PCA for Female Submatrix...\n")
pca_female <- run_gender_pca(X_female)

cat("Executing PCA for Male Submatrix...\n")
pca_male   <- run_gender_pca(X_male)

# Make Scree and Convergence Plots
prepare_diag_data <- function(pca_res, label) {
  list(
    scree = data.frame(PC = 1:10, Var = pca_res$eigenvalues / pca_res$total_var, Group = label),
    conv = data.frame(Iteration = 1:length(pca_res$log), Error = pca_res$log, Group = label)
  )
}

f_diag <- prepare_diag_data(pca_female, "Female")
m_diag <- prepare_diag_data(pca_male, "Male")

# Calculate Cumulative Variance captured for both submatrices
cum_var_female <- sum(f_diag$scree$Var) * 100
cum_var_male   <- sum(m_diag$scree$Var) * 100
var_subtitle   <- sprintf("Cumulative Variance Explained (First 10 PCs) - Female: %.2f%% | Male: %.2f%%", 
                          cum_var_female, cum_var_male)

# Combined Scree Plot
p_scree <- ggplot(rbind(f_diag$scree, m_diag$scree), aes(x = PC, y = Var, color = Group)) +
  geom_line(size = 0.8) + geom_point(size = 2) + 
  scale_x_continuous(breaks = 1:10) +
  theme_bw() +
  labs(title = "Scree Plot Comparison: Female vs Male", 
       subtitle = var_subtitle,
       y = "Prop. Variance Explained", x = "Principal Component")

# Combined Log-Linear Convergence Plot (with threshold line)
p_conv <- ggplot(rbind(f_diag$conv, m_diag$conv), aes(x = Iteration, y = Error, color = Group)) +
  geom_line(size = 0.8) + scale_y_log10() + 
  geom_hline(yintercept = 1e-6, linetype = "dashed", color = "red", size = 0.8) +
  theme_bw() +
  labs(title = "Log-Linear Convergence Comparison", 
       subtitle = "Dashed red line indicates convergence threshold (1e-6)",
       y = "Error (Log Scale)", x = "Iteration") +
  annotate("text", x = max(f_diag$conv$Iteration)*0.2, y = 1.5e-6, 
           label = "Threshold (10^-6)", color = "red", vjust = -0.5, fontface = "italic")

# Display diagnostic plots
grid.arrange(p_scree, p_conv, ncol = 2)

# Structured Loading Plots (Top 5 / Bottom 5)
get_structured_loadings <- function(pca_res, gender_label) {
  lapply(1:2, function(i) {
    top_5 <- order(pca_res$V[, i], decreasing = TRUE)[1:5]
    bot_5 <- order(pca_res$V[, i], decreasing = FALSE)[1:5]
    data.frame(
      PC = paste0("PC", i),
      Title = clean_titles[c(bot_5, top_5)],
      Score = pca_res$V[c(bot_5, top_5), i],
      Type = rep(c("Bottom 5", "Top 5"), each = 5),
      Gender = gender_label,
      Genre = as.character(primary_genre[c(bot_5, top_5)])
    )
  }) %>% bind_rows()
}

loadings_female <- get_structured_loadings(pca_female, "Female")
loadings_male   <- get_structured_loadings(pca_male, "Male")

# Plot 1: Female Loadings (PC1 and PC2)
p_load_female <- ggplot(loadings_female, aes(x = reorder(Title, Score), y = Score, fill = Genre)) +
  geom_bar(stat = "identity", width = 0.75) + 
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = market_colors) +
  theme_bw() + 
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Female Viewership: Top/Bottom Movie Loadings", 
       subtitle = "PC1 and PC2",
       x = NULL, y = "Loadings")

# Plot 2: Male Loadings (PC1 and PC2)
p_load_male <- ggplot(loadings_male, aes(x = reorder(Title, Score), y = Score, fill = Genre)) +
  geom_bar(stat = "identity", width = 0.75) + 
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = market_colors) +
  theme_bw() + 
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(title = "Male Viewership: Top/Bottom Movie Loadings", 
       subtitle = "PC1 and PC2",
       x = NULL, y = "Loadings")

# Display the loading plots sequentially
print(p_load_female)
print(p_load_male)
