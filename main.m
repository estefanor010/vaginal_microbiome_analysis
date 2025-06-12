% vaginal microbiome analysis pipeline 
clear; clc; close all;
addpath('utils');          
addpath('helpers');        
addpath('analysis');       
addpath('plots');    

params = parameters();
vaginal_data = readtable('vaginal.xlsx');

display(vaginal_data);


vaginal_data.Properties.VariableNames

% correlation with pairs % and nugent 
results = analyze_links_and_nugent(vaginal_data, params);

% species pair probabilities analysis 
pair_results = analyze_species_pair_probabilities(vaginal_data, params);

%Individual species link analysis 
species_results = analyze_individual_species_link_percentages(vaginal_data, params);

%Threshold sensitivity analysis (pending)
%sensitivity_results = analyze_link_sensitivity_func(vaginal_data, params);

%(pending) 
% plots.plot_interaction_types_boxplot(results, params);
% plots.plot_link_percentage_by_nugent_category(results, params);
% plots.plot_link_metrics_comparison(results, params); 

fprintf('\Analyses completed\n');


