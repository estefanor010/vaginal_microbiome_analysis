
function link_prob_result = analyze_species_pair_probabilities(vaginal_data, parameters)
    %Computes the probability matrix and presence matrix (pair-wise) 

    % Input
    % vaginal data: tempperal time series
    % parameters: 

    % Output 
    %   link_prob_result (struct): pair-wise prob links, pair-wise prescence count, species names, 

    subject_ids = unique(vaginal_data.SubjectID);
    num_patients = length(subject_ids);
    
    
    all_species_names = {};
    pair_counts = containers.Map(); %store how many times a pair is linked
    pair_presence = containers.Map(); %stores how many times a pair is present
    
    for patient = 1:num_patients
        [rel_abun, ~,~] = process_data(vaginal_data, subject_ids(patient), parameters);
        species_names = vaginal_data.Properties.VariableNames(11:end); %strings... 
        
        num_species = size(rel_abun, 2);
        if num_species < 2
            continue;
        end
        
        % species amount after filtering 
        current_species = species_names(1:num_species); 

        %modify: I might be overcounting species... it should be more like a set...
        all_species_names = union(all_species_names, current_species);
        
        % obtain pair links on a patient 
        for idx1 = 1:num_species
            for idx2 = (idx1+1):num_species
                %print pair indices on the console for debugging
    
                
                pair_species_indices = sort({current_species{idx1}, current_species{idx2}});
                pair_indices = sprintf('%s-%s', pair_species_indices{1}, pair_species_indices{2});
                
                
                % store pair presence in a hash map
                pair_presence(pair_indices) = pair_presence.get(pair_indices, 0) + 1;
                
                % store pair links in a hash map
                has_link = helpers.calculate_links(rel_abun(:,idx1), rel_abun(:,idx2), parameters);
                pair_counts(pair_indices) = pair_counts.get(pair_indices, 0) + has_link;
            end
        end
    end

    %once store the presence and links quantity hashmaps
    %iterates on both to create:
        % a pair-wise probability link matrix (below)
        % a pair-wise presence count matrix 
    
        
    
    num_species = length(all_species_names); %check overcounting species
    link_prob_matrix = zeros(num_species);
    presence_count_matrix = zeros(num_species); 
    
    % Sort species alphabetically
    [sorted_species, idx] = sort(all_species_names);
    
    for idx1 = 1:num_species
        for idx2 = (idx1+1):num_species
            pair_species = sprintf('%s-%s', sorted_species{idx1}, sorted_species{idx2});
            presence = pair_presence.get(pair_species, 0);
            links = pair_counts.get(pair_species, 0);
            
            presence_count_matrix(idx1,idx2) = presence;
            if presence > 0
                link_prob_matrix(idx1,idx2) = links / presence;
            end
        end
    end
    
    % symmetric prob and presence matrix 
    link_prob_matrix = link_prob_matrix + link_prob_matrix';
    presence_count_matrix = presence_count_matrix + presence_count_matrix';
    
    % plot a heatmap (pending function)
    plots.plot_species_pair_heatmap(link_prob_matrix, sorted_species, parameters);
    
    % Store results
    link_prob_result = struct();
    link_prob_result.link_prob_matrix = link_prob_matrix;
    link_prob_result.presence_count_matrix = presence_count_matrix;
    link_prob_result.species_names = sorted_species;
end
