function plot_individual_species_percentage(species_results, parameters)
%  Visualizes species level % 
% Inputs:
%   species_results (struct): species_names, link_percentages,
%   presence_names, possible_links, actual_links
%   parameters      


    %can change species to show 

    total_number = length(species_results.link_percentages);
    num_to_show =  max(20,total_number); 

    


    figure('Position', [100, 100, parameters.plot_width, parameters.plot_height]);
 
    bar_data_horizontal = species_results.link_percentages(1:num_to_show);
    species_names = regexprep(species_results.species_names(1:num_to_show), '_', ' ');
    barh(bar_data_horizontal, 'FaceColor', [0.2 0.4 0.6]);
    %label num_to_show with species_names 
    set(gca, 'YTick', 1:num_to_show, 'YTickLabel', species_names);
    
    % Add presence count annotations
    for species_quantity = 1:num_to_show
       
        text(bar_data_horizontal(species_quantity) + 0.5, species_quantity, ...
            sprintf('presence (Global) = %d', species_results.presence_counts(species_quantity)), ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 10);  
    end

    % Formatting with fixed font sizes
    title('Individual species link %', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Link %', 'FontSize', 12);
    ylabel('Species', 'FontSize', 12);
    
    grid off;
    
    x_value = max(bar_data_horizontal); 
    if x_value == 0
        x_limits = [0 1];
    else 
        x_limits = [0 x_value * 1.15];
    end

    set(gca, 'FontSize', 11, 'YDir', 'reverse', 'Xlim', x_limits);
end
