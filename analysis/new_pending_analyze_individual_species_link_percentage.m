%similar idea of analyze_species_pair_probabilities.m but storage is
%species node instead of pair-wise. 

function species_results = analyze_individual_species_link_percentage(vaginal_data, parameters)
% Computes species level % (Global)
% Inputs:
%   vaginal_data    
%   parameters     
% Outputs:
%   species_results (struct):


    all_species_names = vaginal_data.Properties.VariableNames(11:end);
    
    disp('all_species_names');
    disp(all_species_names);

    num_species_total = length(all_species_names);

    presence_counts = zeros(num_species_total, 1);
    possible_links = zeros(num_species_total, 1);
    actual_links = zeros(num_species_total, 1);

    subject_ids = unique(vaginal_data.SubjectID);
    num_patients = length(subject_ids);

    for patient = 1:num_patients
        % Process patient data using centralized function
        [rel_abun, ~, idx_species] = process_data(vaginal_data, subject_ids(patient), parameters);

        %num_present_species = length(idx_species); consider fixing the
        %idx_species

        num_present_species = size(rel_abun,2);
        if num_present_species < 2
            continue;  
        end

        % track presence globally
        presence_counts(idx_species) = presence_counts(idx_species) + 1;

        for idx1 = 1:num_present_species
            current_species = idx_species(idx1);
            possible_links(current_species) = possible_links(current_species) + (num_present_species - 1);

            for idx2 = (idx1+1):num_present_species
                partner_species = idx_species(idx2);
                
                has_link = compute_links(rel_abun(:, idx1), rel_abun(:, idx2), parameters);
                actual_links(current_species) = actual_links(current_species) + has_link;
                actual_links(partner_species) = actual_links(partner_species) + has_link;
            end
        end
    end

    link_percentages = zeros(num_species_total, 1);
    valid_species = possible_links > 0;
    link_percentages(valid_species) = (actual_links(valid_species) ./ possible_links(valid_species)) * 100;

    disp('actual links num');
    disp(actual_links(valid_species));

    [sorted_species_names, sort_idx] = sort(all_species_names);
    species_results = struct( ...
        'species_names',    {sorted_species_names}, ...
        'link_percentages', link_percentages(sort_idx), ...
        'presence_counts',  presence_counts(sort_idx), ...
        'possible_links',   possible_links(sort_idx), ...
        'actual_links',     actual_links(sort_idx) ...
    );

    has_data = species_results.presence_counts > 0;
    species_results = structfun(@(temp) temp(has_data), species_results, 'UniformOutput', false);


    disp(species_results.link_percentages);

    plot_individual_species_percentage(species_results, parameters);
end
