function plot_species_pair_heatmap(prob_matrix, species_names, parameters)
% Visualizes species interaction probabilities
% Inputs:
%   prob_matrix:     pair-wise prob matrix 
%   species_names:   species names array 
%   parameters:     parameters 

    
    clean_names = regexprep(species_names, '_', ' ');
    figure('Position', [100, 100, parameters.plot_width, parameters.plot_height]);
    % Generate heatmap
    h = heatmap(clean_names, clean_names, prob_matrix);
    h.Title = 'Species Link Probabilities';
    h.Colormap = parula;
    h.ColorLimits = [0 1];  % For probability visualization
    h.FontSize = 12;
    h.ColorbarVisible = 'on';
    h.Title  = 'Species link %';

    h.MissingDataLabel = 'NA';  
    h.MissingDataColor = [1 1 1];

    cb = colorbar('westoutside'); 
    
    %grid lines play with the labels 
    s = struct(h);
    s.Axes.XAxisLocation = 'top';
    s.Axes.YAxisLocation = 'right';
    grid on;
    
   
    set(gca, 'FontSize', 12, 'FontName', 'Arial');
end
