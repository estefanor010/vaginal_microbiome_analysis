__Keywords:__
1) General ODE-Based Inference (GOBI)

__Status:__ 
Under active development and continuous updates. 

__Overview:__
This project performs quantitative analysis of the vaginal microbiome, focusing on pairwise species interactions, species-level link percentages, and their correlations with clinical health metrics such as the Nugent score. The goal is to better understand microbial community dynamics and to aid in diagnosing vaginal health states including healthy, mild, and bacterial vaginosis (B.V.) conditions. 

__References (Inference Methods and Dataset):__
1) Park, S. H., Ha, S., & Kim, J. K. (2023). A general model-based causal inference method overcomes the curse of synchrony and indirect effect. Nature Communications, 14, 39983.
https://www.nature.com/articles/s41467-023-39983-4
   Method download source: https://github.com/Mathbiomed/GOBI
  
3) Gajer, Pawel, et al. "Temporal dynamics of the human vaginal microbiota." Science translational medicine 4.132 (2012): 132ra52-132ra52.

__Structure:__ 
```
vaginal_microbiome_analysis/
│
├── main.m                     
├── vaginal.xlsx               
│
├── utils/                    
│   ├── check_significance.m   % statistical significance to test links using GOBI 
│   ├── parameters.m           % thresholds and figure dimensions 
│   ├── process_data.m         % data preprocessing per patient (relative abundance and filtering)
│
├── analysis/                  
│   ├── analyze_links_and_nugent.m            % link analysis correlated with Nugent scores
│   ├── analyze_species_pair_probabilities.m  % pairwise species links 
│   ├── analyze_species_link_percentages.m    % species-level links 
│   ├── analyze_link_entropy_correlations.m   % (In development) entropy analysis correlated with Nugent scores 
│
├── helpers/                   
│   └── calculate_links.m      % calculates microbial links using RDS from GOBI 
│
└── plots/                    
    ├── plot_link_correlations.m               % Correlations between links and Nugent
    ├── plot_species_pair_heatmap.m             % Heatmaps of species pair interactions
    ├── plot_individual_species_percentage.m    % Species-level link bar charts
    ├── plot_pairs_by_nugent_category.m         % (In development) boxplots by Nugent groups
    ├── plot_link_percentage_by_nugent_category.m % (In development) category comparisons
    ├── plot_interaction_types_boxplot.m        % (In development) link vs no-link distributions
    └── plot_link_metrics_comparison.m           % (In development) threshold metric comparisons

```

__Getting Started:__

1) Clone the vaginal_microbiome_analysis repository.
2) Open utils/parameters.m and adjust parameters as needed for your study.
3) Run main.m to execute the full analysis pipeline.


