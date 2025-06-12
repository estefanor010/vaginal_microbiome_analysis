function [rel_abun_filtered, nugent_scores, keep_species] = process_data(vaginal_data, subject_id, parameters)
    % Initial data processing (normalization and filtering) 
    % Input:
    %   method: 'any' or 'max' (default at max)
    %   threshold: default = 0.03
    %Output:

    rows = vaginal_data.SubjectID == subject_id;


    
    
    abundance = table2array(vaginal_data(rows, 11:end));
    

    nugent_scores = vaginal_data.("NugentScore")(rows);
    relative_abundance = abundance ./ sum(abundance, 2);

    % Apply species filtering
    if strcmpi(parameters.filter_method, 'any')
        keep_species = any(relative_abundance > parameters.abundance_threshold, 1);
    else
        keep_species = max(relative_abundance, [], 1) > parameters.abundance_threshold;
    end

    rel_abun_filtered = relative_abundance(:, keep_species);
end
