
function results = analyze_links_and_nugent(vaginal_data, parameters)
    

    subject_ids = unique(vaginal_data.SubjectID);
    num_patients = length(subject_ids);

    
   
    link_counts = zeros(num_patients, 1);
    total_possible_links = zeros(num_patients, 1);
    link_percentages = zeros(num_patients, 1);
    patient_mean_nugent = zeros(num_patients, 1);
    patient_median_nugent = zeros(num_patients, 1); 

    for patient = 1:num_patients
        subject_id = subject_ids(patient);
        [rel_abun_1_percent, nugent_scores,~] = process_data(vaginal_data, subject_id, parameters);
        
        % display('rel_abun_1_percent');
        % display(rel_abun_1_percent);

        display('nugent_scores');
        display(nugent_scores);
       
        patient_mean_nugent(patient) = mean(nugent_scores, 'omitnan');
        patient_median_nugent(patient) = median(nugent_scores, 'omitnan');  % Median calculation

        num_species = size(rel_abun_1_percent, 2);
        if num_species < 2
            link_percentages(patient) = NaN;
            continue;
        end
        
       
        total_links = (num_species * (num_species - 1)) / 2;
        link_count = 0;
        
        for idx1 = 1:num_species
            for idx2 = (idx1+1):num_species
                link_count = link_count + compute_links(...
                    rel_abun_1_percent(:,idx1), rel_abun_1_percent(:,idx2), parameters);
            end
        end
        
        
        link_counts(patient) = link_count;
        total_possible_links(patient) = total_links;
        link_percentages(patient) = (link_count / total_links) * 100;
    end

    display('link_counts');
    display(link_counts);

    display('total_links');
    display(total_links);
    
    %keep links that are nonzero 
    valid_idx = total_possible_links > 0;
    results.link_percent = link_percentages(valid_idx);
    results.mean_nugent = patient_mean_nugent(valid_idx);
    results.median_nugent = patient_median_nugent(valid_idx);  % Store median
    
    display(results.link_percent);
    display(results.mean_nugent);
    display(results.median_nugent);

    display(results);
    
    % plot correlations (pairs %, nugent type)
    plot_link_correlations(results, parameters);
end
