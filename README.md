1. check_significance
– Input:
score_list (3D array): List of scores for species pairs
thres_noise (scalar): Noise threshold
thres_S (scalar): Score threshold
thres_R (scalar): Region coverage threshold
T (scalar): Number of time points
– Output:
is_significant (logical): Whether the pair is significant
–Description:
Determines if a species pair interaction is significant based on score and region coverage thresholds.


2. analyze_links_and_nugent
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots and summary statistics 
– Description:
Analyzes the correlation between the percentage of linked species pairs and Nugent scores for each patient, producing scatter plots and summary statistics.


3. analyze_link_pairs_nugent
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots and summary statistics
– Description:
Computes the number of significant species pairs per patient and their correlation with Nugent scores, displaying results in correlation plots.


4. analyze_link_sensitivity
–Input:
vaginal_data (table): Vaginal microbiome data
–Output:
plots and summary statistics
–Description:
Performs a sensitivity analysis of detected links across a range of thresholds, correlating link counts with Nugent scores and displaying results in tables and plots.


5. plot_pairs_by_nugent_category
–Input:
vaginal_data (table): Vaginal microbiome data
–Output:
box plots and statistics
–Description:
Creates a box plot of the number of significant species pairs per patient, grouped by Nugent score categories (Normal, Intermediate, BV).


6. analyze_link_percent_sensitivity
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots and summary tables 
–Description:
Sensitivity analysis of the percentage of linked species pairs, correlating these percentages with Nugent scores across multiple thresholds.


7. plot_link_percentage_by_nugent_category
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
box plots and summary statistics
–Description:
Box plot showing the percentage of linked species pairs per patient, grouped by Nugent score categories.


8. plot_interaction_types_boxplot
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
box plots and summary statistics 
– Description:
Box plot comparing the percentage of linked and non-linked species pairs by Nugent score category.


9. analyze_normalized_link_percentages
– Input:
vaginal_data (table): Vaginal microbiome data
Output:
plots/tables 
– Description:
Computes normalized link percentages (scaled to sum to 1) and their correlation with Nugent scores, visualized across thresholds.


10. plot_link_metrics_comparison
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots
– Description:
Compares raw link counts, link percentages, and normalized link percentages across thresholds using box plots.


11. analyze_link_entropy_correlations
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots
– Description:
Analyzes the correlation between link entropy (Shannon entropy of link/no-link probabilities) and Nugent scores, across thresholds.


12. analyze_entropy_link_percentage_correlation
– Input:
vaginal_data (table): Vaginal microbiome data
-- Output:
plots 
– Description:
Correlates link entropy and link percentage across thresholds


13. analyze_entropy_by_threshold
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
plots/correlation tables 
– Description:
Comprehensive analysis of Shannon entropy, link percentages, and Nugent scores across thresholds, with correlation visualizations


14. analyze_species_pair_probabilities
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
link_prob_result (struct): Contains link probability matrices, species names, and summary table
– Description:
Calculates and visualizes the probability that each species pair is linked or not linked, including a heatmap.


15. analyze_species_link_percentages
– Input:
vaginal_data (table): Vaginal microbiome data
– Output:
summary tables
– Description:
Analyzes the link percentages at the single-species level, tracking possible and actual connections for each species
