function link_prob_result = analyze_species_pair_probabilities(vaginal_data, parameters)
    % Computes the probability matrix and presence matrix (pair-wise) 

    % Input
    % vaginal data: tempperal time series
    % parameters: 

    % Output 
    %   link_prob_result (struct): pair-wise prob links, pair-wise prescence count, species names, 


    subject_ids = unique(vaginal_data.SubjectID);
    num_patients = length(subject_ids);
    all_species_names = {};
    pair_counts = containers.Map(); 
    pair_presence = containers.Map(); 
    
    for patient = 1:num_patients
        [rel_abun, ~, ~] = process_data(vaginal_data, subject_ids(patient), parameters);
        species_names = vaginal_data.Properties.VariableNames(11:end); 
        
        num_species = size(rel_abun, 2);
        if num_species < 2
            continue;
        end
        
        % species amount after filtering 
        current_species = species_names(1:num_species); 
        all_species_names = union(all_species_names, current_species);
        
        % obtain pair links on a patient 
        for idx1 = 1:num_species
            for idx2 = (idx1+1):num_species
                pair_species_indices = sort({current_species{idx1}, current_species{idx2}});
                pair_indices = sprintf('%s-%s', pair_species_indices{1}, pair_species_indices{2});
                
                %presence pair y/n
                if isKey(pair_presence, pair_indices)
                    pair_presence(pair_indices) = pair_presence(pair_indices) + 1;
                else
                    pair_presence(pair_indices) = 1;
                end
                
                % presence pair has a link y/n
                has_link = compute_links(rel_abun(:,idx1), rel_abun(:,idx2), parameters);
                if isKey(pair_counts, pair_indices)
                    pair_counts(pair_indices) = pair_counts(pair_indices) + has_link;
                else
                    pair_counts(pair_indices) = has_link;
                end
            end
        end
    end

    
    num_species = length(all_species_names);
    [sorted_species, sort_idx] = sort(all_species_names);
    
    link_prob_matrix = zeros(num_species);
    presence_count_matrix = zeros(num_species); 

    for idx1 = 1:num_species
        for idx2 = (idx1+1):num_species
            pair_species = sprintf('%s-%s', sorted_species{idx1}, sorted_species{idx2});

            %update presence or set to 0 
            if isKey(pair_presence, pair_species)
                presence = pair_presence(pair_species);
                presence_count_matrix(idx1, idx2) = presence;
            else
                presence = 0;
                presence_count_matrix(idx1, idx2) = 0;
            end

            %normalize prob matrix indices or set to NAN
            if presence > 0 && isKey(pair_counts, pair_species)
                link_prob_matrix(idx1, idx2) = pair_counts(pair_species) / presence;
            elseif presence == 0
                link_prob_matrix(idx1, idx2) = NaN;
            end
        end
    end

    %symmetric and set NAN 
    link_prob_matrix = triu(link_prob_matrix) + triu(link_prob_matrix, 1)';
    presence_count_matrix = triu(presence_count_matrix) + triu(presence_count_matrix, 1)';
    
    nan_mask = isnan(triu(link_prob_matrix)) | isnan(triu(link_prob_matrix,1)');
    link_prob_matrix(nan_mask) = NaN;

    
    
    % Generate visualization
    plot_species_pair_heatmap(link_prob_matrix, sorted_species, parameters);

    display(all_species_names);
    
    % Store results
    link_prob_result = struct();
    link_prob_result.link_prob_matrix = link_prob_matrix;
    link_prob_result.presence_count_matrix = presence_count_matrix;
    link_prob_result.species_names = sorted_species;
end
