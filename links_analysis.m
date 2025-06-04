vaginal_data = readtable("vaginal.xlsx");
an
function is_significant = check_significance(score_list, thres_noise, thres_S, thres_R, T)
    is_significant = false;
    for type = 1:2
        score = reshape(score_list(:,:,type), [T, T]);
        loca = find(abs(score) > thres_noise);

        if ~isempty(loca)
            s = sum(score(loca)) / sum(abs(score(loca)));
            r = length(loca) / (T * T / 2);

            disp('scores in check_significance')
            disp(s);

            if abs(s) >= thres_S && r >= thres_R
                is_significant = true;
                break;
            end
        end
    end
end

% Initial task percent correlation with pairs one threshold  

% function analyze_links_and_nugent(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.7;  
%     thres_R = 0.9;  
%     time_interval = 1;
% 
%     
%     link_counts = zeros(1, num_patients);
%     total_possible_links = zeros(1, num_patients);
%     link_percentages = zeros(1, num_patients);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
% 
%     
%     for patient = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         
%         patient_mean_nugent(patient) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(patient) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % normalized by 1% 
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         link_counter = 0;
% 
%         % test pairs direction 
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % test i --> j direction 
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % test j --> i direction 
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 if has_link_ij || has_link_ji
%                     link_counter = link_counter + 1;
%                 end
%             end
%         end
% 
%         
%         total_possible = (num_species * (num_species - 1))/2;
%         
%         link_counts(patient) = link_counter;
%         total_possible_links(patient) = total_possible;
%         if total_possible > 0
%             link_percentages(patient) = (link_counter/total_possible) * 100;
%         else
%             link_percentages(patient) = 0;
%         end
%     end
% 
%     % Remove any patients with no data
%     nonempty_idx = total_possible_links > 0;
%     link_percentages = link_percentages(nonempty_idx);
%     patient_mean_nugent = patient_mean_nugent(nonempty_idx);
%     patient_median_nugent = patient_median_nugent(nonempty_idx);
%     subject_ids = subject_ids(nonempty_idx);
% 
%     % Create correlation plots
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Top-left: Mean Nugent vs Link % (Spearman)
%     subplot(2, 2, 1);
%     scatter(link_percentages, patient_mean_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     [sorted_links, idx] = sort(link_percentages);
%     sorted_nugent = patient_mean_nugent(idx);
%     coef = polyfit(sorted_links, sorted_nugent, 1);
%     xfit = linspace(min(link_percentages), max(link_percentages), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%   
%     [rho, pval] = corr(link_percentages', patient_mean_nugent', 'Type', 'Spearman');
%     text(min(link_percentages) + 5, max(patient_mean_nugent) - 1, ...
%         sprintf('Spearman rho = %.2f\np = %.4f', rho, pval), ...
%         'FontSize', 12);
% 
%     title('Link % vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link %', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Top-right: Mean Nugent vs Link % (Pearson)
%     subplot(2, 2, 2);
%     scatter(link_percentages, patient_mean_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     
%     coef = polyfit(link_percentages, patient_mean_nugent, 1);
%     xfit = linspace(min(link_percentages), max(link_percentages), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     
%     [r, p] = corrcoef(link_percentages, patient_mean_nugent);
%     pearson_corr = r(1,2);
%     pearson_p = p(1,2);
%     text(min(link_percentages) + 5, max(patient_mean_nugent) - 1, ...
%         sprintf('Pearson r = %.2f\np = %.4f', pearson_corr, pearson_p), ...
%         'FontSize', 12);
% 
%     title('Link % vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link %', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Bottom-left: Median Nugent vs Link % (Spearman)
%     subplot(2, 2, 3);
%     scatter(link_percentages, patient_median_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     [sorted_links, idx] = sort(link_percentages);
%     sorted_nugent = patient_median_nugent(idx);
%     coef = polyfit(sorted_links, sorted_nugent, 1);
%     xfit = linspace(min(link_percentages), max(link_percentages), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     [rho, pval] = corr(link_percentages', patient_median_nugent', 'Type', 'Spearman');
%     text(min(link_percentages) + 5, max(patient_median_nugent) - 1, ...
%         sprintf('Spearman rho = %.2f\np = %.4f', rho, pval), ...
%         'FontSize', 12);
% 
%     title('Link % vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link %', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Bottom-right: Median Nugent vs Link % (Pearson)
%     subplot(2, 2, 4);
%     scatter(link_percentages, patient_median_nugent, 70, 'filled', 'b');
%     hold on;

%     coef = polyfit(link_percentages, patient_median_nugent, 1);
%     xfit = linspace(min(link_percentages), max(link_percentages), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     [r, p] = corrcoef(link_percentages, patient_median_nugent);
%     pearson_corr = r(1,2);
%     pearson_p = p(1,2);
%     text(min(link_percentages) + 5, max(patient_median_nugent) - 1, ...
%         sprintf('Pearson r = %.2f\np = %.4f', pearson_corr, pearson_p), ...
%         'FontSize', 12);
% 
%     title('Link % vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link %', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     
%     figure('Position', [100, 100, 800, 600]);
%     low_nugent = [0, 3];    % Healthy
%     mid_nugent = [4, 6];    % Intermediate
%     high_nugent = [7, 10];  % BV
% 
%     % patients grouped by nugent type 
%     low_indices = patient_mean_nugent >= low_nugent(1) & patient_mean_nugent <= low_nugent(2);
%     mid_indices = patient_mean_nugent >= mid_nugent(1) & patient_mean_nugent <= mid_nugent(2);
%     high_indices = patient_mean_nugent >= high_nugent(1) & patient_mean_nugent <= high_nugent(2);
% 
%     % Calculate means for each category
%     low_link_mean = mean(link_percentages(low_indices), 'omitnan');
%     mid_link_mean = mean(link_percentages(mid_indices), 'omitnan');
%     high_link_mean = mean(link_percentages(high_indices), 'omitnan');
% 
%     bar(categorical({'Low Nugent (0-3)', 'Mid Nugent (4-6)', 'High Nugent (7-10)'}), ...
%         [low_link_mean, mid_link_mean, high_link_mean]);
%     title('Mean Link Percentage by Nugent Category', 'FontSize', 14);
%     ylabel('Link Percentage (%)', 'FontSize', 12);
%     grid on;
% 
%     % patient pairs by nugent type 
%     text(1, low_link_mean + 2, sprintf('n=%d', sum(low_indices)), 'HorizontalAlignment', 'center');
%     text(2, mid_link_mean + 2, sprintf('n=%d', sum(mid_indices)), 'HorizontalAlignment', 'center');
%     text(3, high_link_mean + 2, sprintf('n=%d', sum(high_indices)), 'HorizontalAlignment', 'center');
% 
%     
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% analyze_links_and_nugent(vaginal_data);



%% counts pairs with one threshold 
% 
% function analyze_link_pairs_nugent(vaginal_data)
%     
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%
%     thres_noise = 0;
%     thres_S = 0.2;  % Threshold for score
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Initialize arrays for storing results
%     link_counts = zeros(1, num_patients);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
% 
%     for patient = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(patient) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(patient) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         link_count = 0;
% 
%         % Check for links between species pairs
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % If either direction has a link, count it
%                 if has_link_ij || has_link_ji
%                     link_count = link_count + 1;
%                 end
%             end
%         end
% 
%         % Store results
%         link_counts(patient) = link_count;
%     end
% 
%     % Remove any patients with no data
%     valid_indices = ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     link_counts = link_counts(valid_indices);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
%     subject_ids = subject_ids(valid_indices);
% 
%     % Print summary statistics
%     fprintf('\nSummary Statistics:\n');
%     fprintf('Total Patients: %d\n', sum(valid_indices));
%     fprintf('Mean Link Count: %.2f\n', mean(link_counts, 'omitnan'));
%     fprintf('Mean Nugent Score: %.2f\n', mean(patient_mean_nugent, 'omitnan'));
% 
%     % Create 4x4 correlation plots
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % 1. Mean Nugent vs Link Count (Spearman)
%     subplot(2, 2, 1);
%     scatter(link_counts, patient_mean_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     [sorted_links, idx] = sort(link_counts);
%     sorted_nugent = patient_mean_nugent(idx);
%     coef = polyfit(sorted_links, sorted_nugent, 1);
%     xfit = linspace(min(link_counts), max(link_counts), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Calculate and display Spearman correlation
%     [rho, pval] = corr(link_counts', patient_mean_nugent', 'Type', 'Spearman');
%     text(min(link_counts) + 1, max(patient_mean_nugent) - 1, ...
%         sprintf('Spearman rho = %.2f\np = %.4f', rho, pval), ...
%         'FontSize', 12);
% 
%     title('Link Count vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Link Count (Pearson)
%     subplot(2, 2, 2);
%     scatter(link_counts, patient_mean_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     % Add regression line
%     coef = polyfit(link_counts, patient_mean_nugent, 1);
%     xfit = linspace(min(link_counts), max(link_counts), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Calculate and display Pearson correlation
%     [r, p] = corrcoef(link_counts, patient_mean_nugent);
%     pearson_corr = r(1,2);
%     pearson_p = p(1,2);
%     text(min(link_counts) + 1, max(patient_mean_nugent) - 1, ...
%         sprintf('Pearson r = %.2f\np = %.4f', pearson_corr, pearson_p), ...
%         'FontSize', 12);
% 
%     title('Link Count vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Link Count (Spearman)
%     subplot(2, 2, 3);
%     scatter(link_counts, patient_median_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     [sorted_links, idx] = sort(link_counts);
%     sorted_nugent = patient_median_nugent(idx);
%     coef = polyfit(sorted_links, sorted_nugent, 1);
%     xfit = linspace(min(link_counts), max(link_counts), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Calculate and display Spearman correlation
%     [rho, pval] = corr(link_counts', patient_median_nugent', 'Type', 'Spearman');
%     text(min(link_counts) + 1, max(patient_median_nugent) - 1, ...
%         sprintf('Spearman rho = %.2f\np = %.4f', rho, pval), ...
%         'FontSize', 12);
% 
%     title('Link Count vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Link Count (Pearson)
%     subplot(2, 2, 4);
%     scatter(link_counts, patient_median_nugent, 70, 'filled', 'b');
%     hold on;
% 
%     % Add regression line
%     coef = polyfit(link_counts, patient_median_nugent, 1);
%     xfit = linspace(min(link_counts), max(link_counts), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Calculate and display Pearson correlation
%     [r, p] = corrcoef(link_counts, patient_median_nugent);
%     pearson_corr = r(1,2);
%     pearson_p = p(1,2);
%     text(min(link_counts) + 1, max(patient_median_nugent) - 1, ...
%         sprintf('Pearson r = %.2f\np = %.4f', pearson_corr, pearson_p), ...
%         'FontSize', 12);
% 
%     title('Link Count vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% end

% To run the function:
% vaginal_data = readtable("vaginal.xlsx");
% analyze_link_pairs_nugent(vaginal_data);


%% 2nd task: pairs correlated with nugent score 

% function analyze_link_sensitivity(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     all_link_counts = zeros(num_patients, num_thresholds);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
% 
%     % Initialize correlation results table
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanMean', 'SpearmanMeanP', 'PearsonMean', 'PearsonMeanP'});
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(pat) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Store results
%             all_link_counts(pat, t) = link_count;
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     all_link_counts = all_link_counts(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
%     subject_ids = subject_ids(valid_indices);
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         link_counts = all_link_counts(:, t);
% 
%         % Spearman correlation with mean Nugent
%         [rho, pval] = corr(link_counts, patient_mean_nugent', 'Type', 'Spearman');
% 
%         % Pearson correlation with mean Nugent
%         [r, p] = corrcoef(link_counts, patient_mean_nugent');
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanMean(t) = rho;
%         corr_results.SpearmanMeanP(t) = pval;
%         corr_results.PearsonMean(t) = pearson_corr;
%         corr_results.PearsonMeanP(t) = pearson_p;
%     end
% 
%     % Print summary statistics for default threshold (0.2)
%     fprintf('\nSummary Statistics (Threshold = 0.2):\n');
%     fprintf('Total Patients: %d\n', sum(valid_indices));
%     fprintf('Mean Link Count: %.2f\n', mean(all_link_counts(:, 2), 'omitnan'));
%     fprintf('Mean Nugent Score: %.2f\n', mean(patient_mean_nugent, 'omitnan'));
% 
%     % Display correlation results table
%     disp('Correlation Results by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results by Threshold', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create 4x4 correlation plots with multiple regression lines
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Define colors for different thresholds
%     colors = jet(num_thresholds);
% 
%     % 1. Mean Nugent vs Link Count (Spearman)
%     subplot(2, 2, 1);
%     hold on;
% 
%     % Plot scatter for default threshold (0.2)
%     default_idx = 2; % 0.2 is the second threshold
%     scatter(all_link_counts(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_counts = all_link_counts(:, t);
%         [sorted_links, idx] = sort(link_counts);
%         sorted_nugent = patient_mean_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_counts), max(link_counts), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Count vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Link Count (Pearson)
%     subplot(2, 2, 2);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_counts(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_counts = all_link_counts(:, t);
%         coef = polyfit(link_counts, patient_mean_nugent, 1);
%         xfit = linspace(min(link_counts), max(link_counts), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Count vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Link Count (Spearman)
%     subplot(2, 2, 3);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_counts(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_counts = all_link_counts(:, t);
%         [sorted_links, idx] = sort(link_counts);
%         sorted_nugent = patient_median_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_counts), max(link_counts), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Count vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Link Count (Pearson)
%     subplot(2, 2, 4);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_counts(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_counts = all_link_counts(:, t);
%         coef = polyfit(link_counts, patient_median_nugent, 1);
%         xfit = linspace(min(link_counts), max(link_counts), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Count vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Number of Linked Pairs', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Add a common colorbar for all subplots
%     colormap(jet);
%     cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
%     cb.Ticks = linspace(0, 1, num_thresholds);
%     cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false);
%     ylabel(cb, 'Threshold Value', 'FontSize', 12);
% 
%     % Create a figure showing average link counts by threshold
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Calculate mean link counts for each threshold
%     mean_link_counts = mean(all_link_counts, 1);
% 
%     % Create bar chart
%     bar(thresholds, mean_link_counts);
%     title('Mean Number of Linked Pairs by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Mean Number of Linked Pairs', 'FontSize', 12);
%     grid on;
% 
%     % Add text labels with counts
%     for t = 1:num_thresholds
%         text(thresholds(t), mean_link_counts(t) + 0.5, sprintf('%.1f', mean_link_counts(t)), ...
%              'HorizontalAlignment', 'center', 'FontSize', 10);
%     end
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% analyze_link_sensitivity(vaginal_data);
% 
% %% 2.1 box plot pairs as patients categorized (low, mid, high)
% function plot_pairs_by_nugent_category(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.3;      % Threshold for score
%     thres_R = 0.9;      % Threshold for region coverage
%     time_interval = 1;
% 
%     % Initialize arrays for storing results
%     link_counts = zeros(1, num_patients);
%     patient_mean_nugent = zeros(1, num_patients);
%     nugent_categories = cell(1, num_patients);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         % Categorize this patient's Nugent score
%         if patient_mean_nugent(pat) >= 0 && patient_mean_nugent(pat) <= 3
%             nugent_categories{pat} = 'Normal (0-3)';
%         elseif patient_mean_nugent(pat) > 3 && patient_mean_nugent(pat) <= 6
%             nugent_categories{pat} = 'Intermediate (4-6)';
%         elseif patient_mean_nugent(pat) > 6 && patient_mean_nugent(pat) <= 10
%             nugent_categories{pat} = 'BV (7-10)';
%         else
%             nugent_categories{pat} = 'Unknown';
%         end
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         link_count = 0;
% 
%         % Check for links between species pairs
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % If either direction has a link, count it
%                 if has_link_ij || has_link_ji
%                     link_count = link_count + 1;
%                 end
%             end
%         end
% 
%         % Store results
%         link_counts(pat) = link_count;
%     end
% 
%     % Remove any patients with no data
%     valid_indices = ~isnan(patient_mean_nugent) & link_counts > 0;
%     link_counts = link_counts(valid_indices);
%     nugent_categories = nugent_categories(valid_indices);
% 
%     % Convert cell array to categorical for boxplot
%     nugent_cat = categorical(nugent_categories, {'Normal (0-3)', 'Intermediate (4-6)', 'BV (7-10)'});
% 
%     % Create the box plot figure
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Create boxplot
%     boxplot(link_counts', nugent_cat);
% 
%     % Add title and labels
%     title('Number of Pairs by Mean Nugent Score Categories', 'FontSize', 14);
%     xlabel('Mean Nugent Score Category', 'FontSize', 12);
%     ylabel('Number of Pairs', 'FontSize', 12);
% 
%     % Add legend text explaining the categories
%     % annotation('textbox', [0.15, 0.8, 0.3, 0.1], ...
%     %     'String', {
%     %         'Mean Nugent Score Categories:', ...
%     %         '• Normal (0-3): Healthy vaginal microbiome, typically dominated by Lactobacillus species', ...
%     %         '• Intermediate (4-6): Transitional state with mixed microbiota', ...
%     %         '• BV (7-10): Bacterial vaginosis, characterized by diverse anaerobic bacteria'
%     %     }, ...
%     %     'FitBoxToText', 'on', ...
%     %     'BackgroundColor', [0.95 0.95 0.95], ...
%     %     'EdgeColor', [0.7 0.7 0.7]);
% 
%     % Count patients in each category
%     normal_count = sum(strcmp(nugent_categories, 'Normal (0-3)'));
%     intermediate_count = sum(strcmp(nugent_categories, 'Intermediate (4-6)'));
%     bv_count = sum(strcmp(nugent_categories, 'BV (7-10)'));
% 
%     % Add patient counts to the plot
%     text(1, max(link_counts) * 0.9, sprintf('n=%d', normal_count), 'HorizontalAlignment', 'center');
%     text(2, max(link_counts) * 0.9, sprintf('n=%d', intermediate_count), 'HorizontalAlignment', 'center');
%     text(3, max(link_counts) * 0.9, sprintf('n=%d', bv_count), 'HorizontalAlignment', 'center');
% 
%     % Print summary statistics
%     fprintf('\n----- PAIR COUNT BY NUGENT CATEGORY -----\n');
% 
%     normal_links = link_counts(strcmp(nugent_categories, 'Normal (0-3)'));
%     fprintf('Normal (0-3): Mean = %.2f, Median = %.2f (n=%d)\n', ...
%         mean(normal_links), median(normal_links), normal_count);
% 
%     intermediate_links = link_counts(strcmp(nugent_categories, 'Intermediate (4-6)'));
%     fprintf('Intermediate (4-6): Mean = %.2f, Median = %.2f (n=%d)\n', ...
%         mean(intermediate_links), median(intermediate_links), intermediate_count);
% 
%     bv_links = link_counts(strcmp(nugent_categories, 'BV (7-10)'));
%     fprintf('BV (7-10): Mean = %.2f, Median = %.2f (n=%d)\n', ...
%         mean(bv_links), median(bv_links), bv_count);
% 
%     % Save the figure
%     saveas(gcf, 'pairs_by_nugent_boxplot.png');
%     saveas(gcf, 'pairs_by_nugent_boxplot.fig');
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% plot_pairs_by_nugent_category(vaginal_data);

% %% 3rd task: pairs percent correlated with nugent score
% function analyze_link_percent_sensitivity(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     all_link_percentages = zeros(num_patients, num_thresholds);
%     total_possible_pairs = zeros(num_patients, 1);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
%     valid_data = false(1, num_patients);  % Track which patients have valid data
% 
%     % Initialize correlation results table
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanMean', 'SpearmanMeanP', 'PearsonMean', 'PearsonMeanP'});
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(pat) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
%         total_possible_pairs(pat) = total_possible;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate and store percentage
%             if total_possible > 0
%                 all_link_percentages(pat, t) = (link_count / total_possible) * 100;
%             else
%                 all_link_percentages(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data & ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     all_link_percentages = all_link_percentages(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
%     subject_ids = subject_ids(valid_indices);
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
% 
%         % Spearman correlation with mean Nugent
%         [rho, pval] = corr(link_percentages, patient_mean_nugent', 'Type', 'Spearman');
% 
%         % Pearson correlation with mean Nugent
%         [r, p] = corrcoef(link_percentages, patient_mean_nugent');
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanMean(t) = rho;
%         corr_results.SpearmanMeanP(t) = pval;
%         corr_results.PearsonMean(t) = pearson_corr;
%         corr_results.PearsonMeanP(t) = pearson_p;
%     end
% 
%     % Print summary statistics for default threshold (0.2)
%     fprintf('\nSummary Statistics (Threshold = 0.2):\n');
%     fprintf('Total Patients: %d\n', sum(valid_indices));
%     fprintf('Mean Link Percentage: %.2f%%\n', mean(all_link_percentages(:, 2), 'omitnan'));
%     fprintf('Mean Nugent Score: %.2f\n', mean(patient_mean_nugent, 'omitnan'));
% 
%     % Display correlation results table
%     disp('Correlation Results by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results by Threshold (Link Percentage)', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create 4x4 correlation plots with multiple regression lines
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Define colors for different thresholds
%     colors = jet(num_thresholds);
% 
%     % 1. Mean Nugent vs Link Percentage (Spearman)
%     subplot(2, 2, 1);
%     hold on;
% 
%     % Plot scatter for default threshold (0.2)
%     default_idx = 2; % 0.2 is the second threshold
%     scatter(all_link_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         [sorted_links, idx] = sort(link_percentages);
%         sorted_nugent = patient_mean_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Link Percentage (Pearson)
%     subplot(2, 2, 2);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         coef = polyfit(link_percentages, patient_mean_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Link Percentage (Spearman)
%     subplot(2, 2, 3);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         [sorted_links, idx] = sort(link_percentages);
%         sorted_nugent = patient_median_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Link Percentage (Pearson)
%     subplot(2, 2, 4);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         coef = polyfit(link_percentages, patient_median_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Add a common colorbar for all subplots
%     colormap(jet);
%     cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
%     cb.Ticks = linspace(0, 1, num_thresholds);
%     cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false);
%     ylabel(cb, 'Threshold Value', 'FontSize', 12);
% 
%     % Create a figure showing average link percentages by threshold
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Calculate mean link percentages for each threshold
%     mean_link_percentages = mean(all_link_percentages, 1);
% 
%     % Create bar chart
%     bar(thresholds, mean_link_percentages);
%     title('Mean Link Percentage by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Mean Link Percentage (%)', 'FontSize', 12);
%     grid on;
% 
%     % Add text labels with percentages
%     for t = 1:num_thresholds
%         text(thresholds(t), mean_link_percentages(t) + 2, sprintf('%.1f%%', mean_link_percentages(t)), ...
%              'HorizontalAlignment', 'center', 'FontSize', 10);
%     end
% 
%     % Create a figure showing link percentages by Nugent category
%     figure('Position', [100, 100, 1200, 600]);
% 
%     % Define Nugent score categories
%     low_nugent = [0, 3];    % Healthy
%     mid_nugent = [4, 6];    % Intermediate
%     high_nugent = [7, 10];  % BV
% 
%     % Categorize patients
%     low_indices = patient_mean_nugent >= low_nugent(1) & patient_mean_nugent <= low_nugent(2);
%     mid_indices = patient_mean_nugent >= mid_nugent(1) & patient_mean_nugent <= mid_nugent(2);
%     high_indices = patient_mean_nugent >= high_nugent(1) & patient_mean_nugent <= high_nugent(2);
% 
%     % Calculate means for each category and threshold
%     nugent_categories = {'Low (0-3)', 'Mid (4-6)', 'High (7-10)'};
%     mean_by_category = zeros(3, num_thresholds);
% 
%     for t = 1:num_thresholds
%         mean_by_category(1, t) = mean(all_link_percentages(low_indices, t), 'omitnan');
%         mean_by_category(2, t) = mean(all_link_percentages(mid_indices, t), 'omitnan');
%         mean_by_category(3, t) = mean(all_link_percentages(high_indices, t), 'omitnan');
%     end
% 
%     % Create grouped bar chart
%     bar(mean_by_category);
%     title('Mean Link Percentage by Nugent Category and Threshold', 'FontSize', 14);
%     xlabel('Nugent Category', 'FontSize', 12);
%     ylabel('Mean Link Percentage (%)', 'FontSize', 12);
%     set(gca, 'XTick', 1:3, 'XTickLabel', nugent_categories);
%     legend(arrayfun(@(x) sprintf('Threshold = %.1f', x), thresholds, 'UniformOutput', false), 'Location', 'best');
%     grid on;
% 
%     % Add patient counts
%     text(1, max(mean_by_category(1,:)) + 5, sprintf('n=%d', sum(low_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
%     text(2, max(mean_by_category(2,:)) + 5, sprintf('n=%d', sum(mid_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
%     text(3, max(mean_by_category(3,:)) + 5, sprintf('n=%d', sum(high_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
% end
% 
% 
% 
% 
% %To run the function:
% vaginal_data = readtable("vaginal.xlsx");
% analyze_link_percent_sensitivity(vaginal_data);

%% 3.1 box plot percentage by nugent type 

% function plot_link_percentage_by_nugent_category(vaginal_data)
%     % This function creates a box plot showing the percentage of linked species pairs
%     % for each Nugent score category (low, intermediate, high)
% 
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.3;      % Threshold for score (using 0.3 as requested)
%     thres_R = 0.9;      % Threshold for region coverage
%     time_interval = 1;
% 
%     % Initialize arrays for storing results
%     link_counts = zeros(1, num_patients);
%     total_possible_links = zeros(1, num_patients);
%     link_percentages = zeros(1, num_patients);
%     patient_mean_nugent = zeros(1, num_patients);
%     nugent_categories = cell(1, num_patients);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         % Categorize this patient's Nugent score
%         if patient_mean_nugent(pat) >= 0 && patient_mean_nugent(pat) <= 3
%             nugent_categories{pat} = 'Normal (0-3)';
%         elseif patient_mean_nugent(pat) > 3 && patient_mean_nugent(pat) <= 6
%             nugent_categories{pat} = 'Intermediate (4-6)';
%         elseif patient_mean_nugent(pat) > 6 && patient_mean_nugent(pat) <= 10
%             nugent_categories{pat} = 'BV (7-10)';
%         else
%             nugent_categories{pat} = 'Unknown';
%         end
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         link_count = 0;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
%         total_possible_links(pat) = total_possible;
% 
%         % Check for links between species pairs
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % If either direction has a link, count it
%                 if has_link_ij || has_link_ji
%                     link_count = link_count + 1;
%                 end
%             end
%         end
% 
%         % Store results
%         link_counts(pat) = link_count;
% 
%         % Calculate percentage of linked pairs
%         if total_possible > 0
%             link_percentages(pat) = (link_count/total_possible) * 100;
%         else
%             link_percentages(pat) = 0;
%         end
%     end
% 
%     % Remove any patients with no data or no possible links
%     valid_indices = ~isnan(patient_mean_nugent) & total_possible_links > 0;
%     link_percentages = link_percentages(valid_indices);
%     nugent_categories = nugent_categories(valid_indices);
% 
%     % Convert cell array to categorical for boxplot
%     nugent_cat = categorical(nugent_categories, {'Normal (0-3)', 'Intermediate (4-6)', 'BV (7-10)'});
% 
%     % Create the box plot figure
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Create boxplot
%     boxplot(link_percentages', nugent_cat);
% 
%     % Add title and labels
%     title('Percentage of Linked Pairs by Mean Nugent Score Categories (Threshold = 0.3)', 'FontSize', 14);
%     xlabel('Mean Nugent Score Category', 'FontSize', 12);
%     ylabel('Percentage of Linked Pairs (%)', 'FontSize', 12);
% 
%     % Count patients in each category
%     normal_count = sum(strcmp(nugent_categories, 'Normal (0-3)'));
%     intermediate_count = sum(strcmp(nugent_categories, 'Intermediate (4-6)'));
%     bv_count = sum(strcmp(nugent_categories, 'BV (7-10)'));
% 
%     % Add patient counts to the plot
%     text(1, max(link_percentages) * 0.9, sprintf('n=%d', normal_count), 'HorizontalAlignment', 'center');
%     text(2, max(link_percentages) * 0.9, sprintf('n=%d', intermediate_count), 'HorizontalAlignment', 'center');
%     text(3, max(link_percentages) * 0.9, sprintf('n=%d', bv_count), 'HorizontalAlignment', 'center');
% 
%     % Print summary statistics
%     fprintf('\n----- LINK PERCENTAGE BY NUGENT CATEGORY (Threshold = 0.3) -----\n');
% 
%     normal_percentages = link_percentages(strcmp(nugent_categories, 'Normal (0-3)'));
%     fprintf('Normal (0-3): Mean = %.2f%%, Median = %.2f%% (n=%d)\n', ...
%         mean(normal_percentages), median(normal_percentages), normal_count);
% 
%     intermediate_percentages = link_percentages(strcmp(nugent_categories, 'Intermediate (4-6)'));
%     fprintf('Intermediate (4-6): Mean = %.2f%%, Median = %.2f%% (n=%d)\n', ...
%         mean(intermediate_percentages), median(intermediate_percentages), intermediate_count);
% 
%     bv_percentages = link_percentages(strcmp(nugent_categories, 'BV (7-10)'));
%     fprintf('BV (7-10): Mean = %.2f%%, Median = %.2f%% (n=%d)\n', ...
%         mean(bv_percentages), median(bv_percentages), bv_count);
% 
%     % Save the figure
%     saveas(gcf, 'link_percentage_by_nugent_boxplot.png');
%     saveas(gcf, 'link_percentage_by_nugent_boxplot.fig');
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% plot_link_percentage_by_nugent_category(vaginal_data);

%% 3.2 box plot percent link and no link by nugent category 
% 
% function plot_interaction_types_boxplot(vaginal_data)
%     % This function creates a boxplot showing the distribution of linked and non-linked
%     % species pairs percentages for each Nugent score category
% 
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.3;      % Threshold for score
%     thres_R = 0.9;      % Threshold for region coverage
%     time_interval = 1;
% 
%     % Initialize arrays for storing results
%     link_percentages = zeros(1, num_patients);
%     nolink_percentages = zeros(1, num_patients);
%     patient_mean_nugent = zeros(1, num_patients);
%     nugent_categories = cell(1, num_patients);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         % Categorize this patient's Nugent score
%         if patient_mean_nugent(pat) >= 0 && patient_mean_nugent(pat) <= 3
%             nugent_categories{pat} = 'Low';
%         elseif patient_mean_nugent(pat) > 3 && patient_mean_nugent(pat) <= 6
%             nugent_categories{pat} = 'Mid';
%         elseif patient_mean_nugent(pat) > 6 && patient_mean_nugent(pat) <= 10
%             nugent_categories{pat} = 'High';
%         else
%             nugent_categories{pat} = 'Unknown';
%         end
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         link_count = 0;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
% 
%         % Check for links between species pairs
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % If either direction has a link, count it
%                 if has_link_ij || has_link_ji
%                     link_count = link_count + 1;
%                 end
%             end
%         end
% 
%         % Calculate percentages
%         if total_possible > 0
%             link_percentages(pat) = (link_count/total_possible) * 100;
%             nolink_percentages(pat) = ((total_possible - link_count)/total_possible) * 100;
%         else
%             link_percentages(pat) = 0;
%             nolink_percentages(pat) = 0;
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = ~isnan(patient_mean_nugent) & (link_percentages + nolink_percentages > 0);
%     link_percentages = link_percentages(valid_indices);
%     nolink_percentages = nolink_percentages(valid_indices);
%     nugent_categories = nugent_categories(valid_indices);
% 
%     % Group data by Nugent category
%     low_indices = strcmp(nugent_categories, 'Low');
%     mid_indices = strcmp(nugent_categories, 'Mid');
%     high_indices = strcmp(nugent_categories, 'High');
% 
%     % Create data for boxplot
%     data = [];
%     group = [];
%     colors = [];
% 
%     % Define colors for Nugent categories
%     low_color = [0, 0, 1];    % Blue for Low Nugent
%     mid_color = [1, 0.5, 0];  % Orange for Mid Nugent
%     high_color = [1, 0.8, 0]; % Yellow for High Nugent
% 
%     % Add linked data
%     if sum(low_indices) > 0
%         data = [data; link_percentages(low_indices)'];
%         group = [group; ones(sum(low_indices), 1)];
%         colors = [colors; repmat(low_color, sum(low_indices), 1)];
%     end
% 
%     if sum(mid_indices) > 0
%         data = [data; link_percentages(mid_indices)'];
%         group = [group; 2*ones(sum(mid_indices), 1)];
%         colors = [colors; repmat(mid_color, sum(mid_indices), 1)];
%     end
% 
%     if sum(high_indices) > 0
%         data = [data; link_percentages(high_indices)'];
%         group = [group; 3*ones(sum(high_indices), 1)];
%         colors = [colors; repmat(high_color, sum(high_indices), 1)];
%     end
% 
%     % Add non-linked data
%     if sum(low_indices) > 0
%         data = [data; nolink_percentages(low_indices)'];
%         group = [group; 4*ones(sum(low_indices), 1)];
%         colors = [colors; repmat(low_color, sum(low_indices), 1)];
%     end
% 
%     if sum(mid_indices) > 0
%         data = [data; nolink_percentages(mid_indices)'];
%         group = [group; 5*ones(sum(mid_indices), 1)];
%         colors = [colors; repmat(mid_color, sum(mid_indices), 1)];
%     end
% 
%     if sum(high_indices) > 0
%         data = [data; nolink_percentages(high_indices)'];
%         group = [group; 6*ones(sum(high_indices), 1)];
%         colors = [colors; repmat(high_color, sum(high_indices), 1)];
%     end
% 
%     % Create figure
%     figure('Position', [100, 100, 1000, 600]);
% 
%     % Create custom labels for x-axis
%     %x_labels = {'Link-Low', 'Link-Mid', 'Link-High', 'NoLink-Low', 'NoLink-Mid', 'NoLink-High'};
% 
%     % Create boxplot
%     % h = boxplot(data, group, 'Colors', [low_color; mid_color; high_color; low_color; mid_color; high_color], ...
%     %     'Labels', x_labels, 'LabelOrientation', 'horizontal');
% 
%     h = boxplot(data, group, 'Colors', [low_color; mid_color; high_color; low_color; mid_color; high_color]);
% 
%     % Customize boxplot appearance
%     set(gca, 'FontSize', 12);
% 
%     % Add title and labels
%     title('Distribution of Interaction Types by Nugent Category', 'FontSize', 16);
%     ylabel('Percentage (%)', 'FontSize', 14);
%     ylim([0 100]);
% 
%     % Create a custom legend
%     hold on;
%     h1 = plot(NaN, NaN, 'Color', low_color, 'LineWidth', 2);
%     h2 = plot(NaN, NaN, 'Color', mid_color, 'LineWidth', 2);
%     h3 = plot(NaN, NaN, 'Color', high_color, 'LineWidth', 2);
%     legend([h1, h2, h3], {'Low Nugent (0-3)', 'Mid Nugent (4-6)', 'High Nugent (7-10)'}, ...
%         'Location', 'northeast', 'FontSize', 12);
% 
% 
% 
%     % Add text labels for the two main groups
%     text(2, -10, 'Link', 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
%     text(5, -10, 'No Link', 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% 
%     % Add patient counts
%     text(1, 95, sprintf('n=%d', sum(low_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
%     text(2, 95, sprintf('n=%d', sum(mid_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
%     text(3, 95, sprintf('n=%d', sum(high_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
%     text(4, 95, sprintf('n=%d', sum(low_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
%     text(5, 95, sprintf('n=%d', sum(mid_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
%     text(6, 95, sprintf('n=%d', sum(high_indices)), 'HorizontalAlignment', 'center', 'FontSize', 10);
% 
%     % Print summary statistics
%     fprintf('\n----- INTERACTION TYPES BY NUGENT CATEGORY -----\n');
%     fprintf('Low Nugent (0-3): Linked = %.2f%%, Non-linked = %.2f%% (n=%d)\n', ...
%         mean(link_percentages(low_indices), 'omitnan'), mean(nolink_percentages(low_indices), 'omitnan'), sum(low_indices));
%     fprintf('Mid Nugent (4-6): Linked = %.2f%%, Non-linked = %.2f%% (n=%d)\n', ...
%         mean(link_percentages(mid_indices), 'omitnan'), mean(nolink_percentages(mid_indices), 'omitnan'), sum(mid_indices));
%     fprintf('High Nugent (7-10): Linked = %.2f%%, Non-linked = %.2f%% (n=%d)\n', ...
%         mean(link_percentages(high_indices), 'omitnan'), mean(nolink_percentages(high_indices), 'omitnan'), sum(high_indices));
% 
%     % Save the figure
%     saveas(gcf, 'interaction_types_boxplot.png');
%     saveas(gcf, 'interaction_types_boxplot.fig');
% end
% 
% plot_interaction_types_boxplot(vaginal_data);


%% 4th normalized percentages of pairs correlated with nugent score 

% function analyze_normalized_link_percentages(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     all_link_percentages = zeros(num_patients, num_thresholds);
%     total_possible_pairs = zeros(num_patients, 1);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
%     valid_data = false(1, num_patients);  % Track which patients have valid data
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(pat) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
%         total_possible_pairs(pat) = total_possible;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate and store percentage
%             if total_possible > 0
%                 all_link_percentages(pat, t) = (link_count / total_possible) * 100;
%             else
%                 all_link_percentages(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data & ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     all_link_percentages = all_link_percentages(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
%     subject_ids = subject_ids(valid_indices);
% 
%     % Initialize array for normalized percentages
%     all_normalized_percentages = zeros(size(all_link_percentages));
% 
%     % Normalize link percentages for each threshold so they sum to 1 (0-1 scale)
%     for t = 1:num_thresholds
%         % Get the link percentages for this threshold
%         link_percentages = all_link_percentages(:, t);
% 
%         % Calculate the sum of all percentages
%         total_percentage = sum(link_percentages);
% 
%         % Normalize to make them sum to 1
%         if total_percentage > 0
%             normalized_percentages = link_percentages / total_percentage;  % Sum to 1
%             all_normalized_percentages(:, t) = normalized_percentages;
%         else
%             all_normalized_percentages(:, t) = zeros(size(link_percentages));
%         end
%     end
% 
%     % Print summary statistics for default threshold (0.2)
%     fprintf('\nSummary Statistics (Threshold = 0.2):\n');
%     fprintf('Total Patients: %d\n', sum(valid_indices));
%     fprintf('Mean Normalized Link Percentage: %.4f\n', mean(all_normalized_percentages(:, 2), 'omitnan'));
%     fprintf('Mean Nugent Score: %.2f\n', mean(patient_mean_nugent, 'omitnan'));
% 
%     % Calculate correlations for normalized percentages
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanMean', 'SpearmanMeanP', 'PearsonMean', 'PearsonMeanP'});
% 
%     for t = 1:num_thresholds
%         norm_percentages = all_normalized_percentages(:, t);
% 
%         % Spearman correlation with mean Nugent
%         [rho, pval] = corr(norm_percentages, patient_mean_nugent', 'Type', 'Spearman');
% 
%         % Pearson correlation with mean Nugent
%         [r, p] = corrcoef(norm_percentages, patient_mean_nugent');
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanMean(t) = rho;
%         corr_results.SpearmanMeanP(t) = pval;
%         corr_results.PearsonMean(t) = pearson_corr;
%         corr_results.PearsonMeanP(t) = pearson_p;
%     end
% 
%     % Display correlation results table
%     disp('Correlation Results for Normalized Percentages by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results by Threshold (Normalized Link Percentage)', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create 4x4 correlation plots with multiple regression lines for normalized percentages
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Define colors for different thresholds
%     colors = jet(num_thresholds);
% 
%     % 1. Mean Nugent vs Normalized Link Percentage (Spearman)
%     subplot(2, 2, 1);
%     hold on;
% 
%     % Plot scatter for default threshold (0.2)
%     default_idx = 2; % 0.2 is the second threshold
%     scatter(all_normalized_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         norm_percentages = all_normalized_percentages(:, t);
%         [sorted_links, idx] = sort(norm_percentages);
%         sorted_nugent = patient_mean_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(norm_percentages), max(norm_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Normalized Link Percentage vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Normalized Link Percentage (Pearson)
%     subplot(2, 2, 2);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_normalized_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         norm_percentages = all_normalized_percentages(:, t);
%         coef = polyfit(norm_percentages, patient_mean_nugent, 1);
%         xfit = linspace(min(norm_percentages), max(norm_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Normalized Link Percentage vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Normalized Link Percentage (Spearman)
%     subplot(2, 2, 3);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_normalized_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         norm_percentages = all_normalized_percentages(:, t);
%         [sorted_links, idx] = sort(norm_percentages);
%         sorted_nugent = patient_median_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(norm_percentages), max(norm_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Normalized Link Percentage vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Normalized Link Percentage (Pearson)
%     subplot(2, 2, 4);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_normalized_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         norm_percentages = all_normalized_percentages(:, t);
%         coef = polyfit(norm_percentages, patient_median_nugent, 1);
%         xfit = linspace(min(norm_percentages), max(norm_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Normalized Link Percentage vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Add a common colorbar for all subplots
%     colormap(jet);
%     cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
%     cb.Ticks = linspace(0, 1, num_thresholds);
%     cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false);
%     ylabel(cb, 'Threshold Value', 'FontSize', 12);
% 
%     % Create a figure showing average normalized link percentages by threshold
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Calculate mean normalized link percentages for each threshold
%     mean_norm_percentages = mean(all_normalized_percentages, 1);
% 
%     % Create bar chart
%     bar(thresholds, mean_norm_percentages);
%     title('Mean Normalized Link Percentage by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Mean Normalized Link Percentage (0-1)', 'FontSize', 12);
%     grid on;
% 
%     % Add text labels with percentages
%     for t = 1:num_thresholds
%         text(thresholds(t), mean_norm_percentages(t) + 0.02, sprintf('%.3f', mean_norm_percentages(t)), ...
%              'HorizontalAlignment', 'center', 'FontSize', 10);
%     end
% 
%     % Create a comparison figure showing original vs normalized percentages
%     figure('Position', [100, 100, 1200, 500]);
% 
%     % Subplot for original percentages
%     subplot(1, 2, 1);
%     boxplot(all_link_percentages, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     title('Original Link Percentages by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Link Percentage (%)', 'FontSize', 12);
%     grid on;
% 
%     % Subplot for normalized percentages
%     subplot(1, 2, 2);
%     boxplot(all_normalized_percentages, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     title('Normalized Link Percentages by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
%     grid on;
% end

% To run the function:
% vaginal_data = readtable("vaginal.xlsx");
% analyze_normalized_link_percentages(vaginal_data);



%% 5th box plots of thresholds and pairs, percents, normalized links 
%% Threshold sensitivity analysis 
% 
 
% 
% function plot_link_metrics_comparison(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     all_link_counts = zeros(num_patients, num_thresholds);
%     all_link_percentages = zeros(num_patients, num_thresholds);
%     total_possible_pairs = zeros(num_patients, 1);
%     valid_data = false(1, num_patients);  % Track which patients have valid data
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
%         total_possible_pairs(pat) = total_possible;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Store raw count
%             all_link_counts(pat, t) = link_count;
% 
%             % Calculate and store percentage
%             if total_possible > 0
%                 all_link_percentages(pat, t) = (link_count / total_possible) * 100;
%             else
%                 all_link_percentages(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data;
%     all_link_counts = all_link_counts(valid_indices, :);
%     all_link_percentages = all_link_percentages(valid_indices, :);
%     total_possible_pairs = total_possible_pairs(valid_indices);
% 
%     % Initialize array for normalized percentages
%     all_normalized_percentages = zeros(size(all_link_percentages));
% 
%     % Normalize link percentages for each threshold so they sum to 1 (0-1 scale)
%     for t = 1:num_thresholds
%         % Get the link percentages for this threshold
%         link_percentages = all_link_percentages(:, t);
% 
%         % Calculate the sum of all percentages
%         total_percentage = sum(link_percentages);
% 
%         % Normalize to make them sum to 1
%         if total_percentage > 0
%             normalized_percentages = link_percentages / total_percentage;
%             all_normalized_percentages(:, t) = normalized_percentages;
%         else
%             all_normalized_percentages(:, t) = zeros(size(link_percentages));
%         end
%     end
% 
%     % Create a figure with 3 subplots showing boxplots for each metric
%     figure('Position', [100, 100, 1800, 600]);
% 
%     % 1. Raw Link Counts
%     subplot(1, 3, 1);
%     boxplot(all_link_counts, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     title('Raw Number of Linked Pairs by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Number of Linked Pairs', 'FontSize', 12);
% 
% 
%     % 2. Link Percentages
%     subplot(1, 3, 2);
%     boxplot(all_link_percentages, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     title('Link Percentage by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Link Percentage (%)', 'FontSize', 12);
% 
% 
%     % 3. Normalized Link Percentages
%     subplot(1, 3, 3);
%     boxplot(all_normalized_percentages, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     title('Normalized Link Percentage by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Normalized Link Percentage (0-1)', 'FontSize', 12);
% 
% end

% To run the function:
% vaginal_data = readtable("vaginal.xlsx");
% plot_link_metrics_comparison(vaginal_data);



% 6th entropy correlated with nugent 

% function analyze_link_entropy_correlations(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     link_entropy = zeros(num_patients, num_thresholds);
%     patient_mean_nugent = zeros(num_patients, 1);
%     patient_median_nugent = zeros(num_patients, 1);
%     valid_data = false(num_patients, 1);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
%         nugent_scores = vaginal_data.NugentScore(rows);
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(nugent_scores, 'omitnan');
%         patient_median_nugent(pat) = median(nugent_scores, 'omitnan');
% 
%         % Skip patients with insufficient data
%         if size(abun, 1) < 3 || sum(~isnan(nugent_scores)) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Calculate total possible pairs
%             total_possible = (num_species * (num_species - 1))/2;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate probability of link and no link
%             if total_possible > 0
%                 p_link = link_count / total_possible;
%                 p_no_link = 1 - p_link;
% 
%                 % Calculate Shannon entropy: -sum(p_i * log(p_i))
%                 % Only include non-zero probabilities in calculation
%                 entropy = 0;
%                 if p_link > 0
%                     entropy = entropy - p_link * log2(p_link);
%                 end
%                 if p_no_link > 0
%                     entropy = entropy - p_no_link * log2(p_no_link);
%                 end
% 
%                 link_entropy(pat, t) = entropy;
%             else
%                 link_entropy(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data & ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     link_entropy = link_entropy(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
% 
%     % Initialize correlation results table
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanMean', 'SpearmanMeanP', 'PearsonMean', 'PearsonMeanP'});
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
% 
%         % Spearman correlation with mean Nugent
%         [rho, pval] = corr(entropy_values, patient_mean_nugent, 'Type', 'Spearman');
% 
%         % Pearson correlation with mean Nugent
%         [r, p] = corrcoef(entropy_values, patient_mean_nugent);
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanMean(t) = rho;
%         corr_results.SpearmanMeanP(t) = pval;
%         corr_results.PearsonMean(t) = pearson_corr;
%         corr_results.PearsonMeanP(t) = pearson_p;
%     end
% 
%     % Find the threshold with the highest absolute Spearman correlation
%     [~, max_spearman_idx] = max(abs(corr_results.SpearmanMean));
%     max_spearman_threshold = thresholds(max_spearman_idx);
%     max_spearman_rho = corr_results.SpearmanMean(max_spearman_idx);
%     max_spearman_p = corr_results.SpearmanMeanP(max_spearman_idx);
% 
%     % Display correlation results table
%     disp('Correlation Results for Link Entropy by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results by Threshold (Link Entropy)', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create 4x4 correlation plots with multiple regression lines for entropy
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Define colors for different thresholds
%     colors = jet(num_thresholds);
% 
%     % 1. Mean Nugent vs Link Entropy (Spearman)
%     subplot(2, 2, 1);
%     hold on;
% 
%     % Plot scatter for default threshold (0.2)
%     default_idx = 2; % 0.2 is the second threshold
%     scatter(link_entropy(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
%         [sorted_entropy, idx] = sort(entropy_values);
%         sorted_nugent = patient_mean_nugent(idx);
%         coef = polyfit(sorted_entropy, sorted_nugent, 1);
%         xfit = linspace(min(entropy_values), max(entropy_values), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Entropy vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Entropy', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Link Entropy (Pearson)
%     subplot(2, 2, 2);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(link_entropy(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
%         coef = polyfit(entropy_values, patient_mean_nugent, 1);
%         xfit = linspace(min(entropy_values), max(entropy_values), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Entropy vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Entropy (bits)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Link Entropy (Spearman)
%     subplot(2, 2, 3);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(link_entropy(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
%         [sorted_entropy, idx] = sort(entropy_values);
%         sorted_nugent = patient_median_nugent(idx);
%         coef = polyfit(sorted_entropy, sorted_nugent, 1);
%         xfit = linspace(min(entropy_values), max(entropy_values), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Entropy vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Entropy', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Link Entropy (Pearson)
%     subplot(2, 2, 4);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(link_entropy(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
%         coef = polyfit(entropy_values, patient_median_nugent, 1);
%         xfit = linspace(min(entropy_values), max(entropy_values), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Entropy vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Entropy (bits)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Add a common colorbar for all subplots
%     colormap(jet);
%     cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
%     cb.Ticks = linspace(0, 1, num_thresholds);
%     cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false);
%     ylabel(cb, 'Threshold Value', 'FontSize', 12);
% 
%     % Create a figure showing average link entropy by threshold
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Calculate mean link entropy for each threshold
%     mean_link_entropy = mean(link_entropy, 1);
% 
%     % Create bar chart
%     bar(thresholds, mean_link_entropy);
%     title('Mean Link Entropy by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Mean Link Entropy', 'FontSize', 12);
%     grid on;
% 
%     % Add text labels with entropy values
%     for t = 1:num_thresholds
%         text(thresholds(t), mean_link_entropy(t) + 0.02, sprintf('%.3f', mean_link_entropy(t)), ...
%              'HorizontalAlignment', 'center', 'FontSize', 10);
%     end
% 
%     % Create a figure for the highest Spearman correlation
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Get entropy values for the threshold with highest correlation
%     entropy_values = link_entropy(:, max_spearman_idx);
% 
%     % Create scatter plot
%     scatter(entropy_values, patient_mean_nugent, 100, 'filled', 'b');
%     hold on;
% 
%     % Add trend line (Spearman)
%     [sorted_entropy, idx] = sort(entropy_values);
%     sorted_nugent = patient_mean_nugent(idx);
%     coef = polyfit(sorted_entropy, sorted_nugent, 1);
%     xfit = linspace(min(entropy_values), max(entropy_values), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 3);
% 
%     % Add correlation statistics
%     text(min(entropy_values) + 0.1*(max(entropy_values)-min(entropy_values)), ...
%          max(patient_mean_nugent) - 0.2*(max(patient_mean_nugent)-min(patient_mean_nugent)), ...
%          sprintf('Spearman rho = %.3f\np-value = %.4f', max_spearman_rho, max_spearman_p), ...
%          'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
% 
%     title('Highest Correlation: Link Entropy vs Mean Nugent Score', 'FontSize', 16, 'FontWeight', 'bold');
%     xlabel('Link Entropy', 'FontSize', 14);
%     ylabel('Mean Nugent Score', 'FontSize', 14);
%     grid on;
% 
%     % Print summary statistics
%     fprintf('\nSummary Statistics:\n');
%     fprintf('Total Patients with Valid Data: %d\n', sum(valid_indices));
% 
%     % Print entropy statistics for threshold = 0.2
%     fprintf('\nLink Entropy Statistics (Threshold = 0.2):\n');
%     fprintf('Mean Link Entropy: %.3f bits\n', mean(link_entropy(:, 2)));
%     fprintf('Min Link Entropy: %.3f bits\n', min(link_entropy(:, 2)));
%     fprintf('Max Link Entropy: %.3f bits\n', max(link_entropy(:, 2)));
% 
%     % Print correlation statistics for threshold = 0.2
%     fprintf('\nCorrelations for Link Entropy (Threshold = 0.2):\n');
%     [r, p] = corrcoef(link_entropy(:, 2), patient_mean_nugent);
%     fprintf('With Mean Nugent - Pearson: r = %.3f (p = %.4f)\n', r(1,2), p(1,2));
% 
%     [rho, pval] = corr(link_entropy(:, 2), patient_mean_nugent, 'Type', 'Spearman');
%     fprintf('With Mean Nugent - Spearman: rho = %.3f (p = %.4f)\n', rho, pval);
% 
%     % Print highest correlation statistics
%     fprintf('\nHighest Correlation:\n');
%     fprintf('Threshold: %.1f\n', max_spearman_threshold);
%     fprintf('Spearman rho: %.3f (p = %.4f)\n', max_spearman_rho, max_spearman_p);
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% analyze_link_entropy_correlations(vaginal_data);
% 



%% 7th entropy correlated with links
% 
% function analyze_entropy_link_percentage_correlation(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     link_entropy = zeros(num_patients, num_thresholds);
%     link_percentages = zeros(num_patients, num_thresholds);
%     valid_data = false(num_patients, 1);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Skip patients with insufficient data
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Calculate total possible pairs
%             total_possible = (num_species * (num_species - 1))/2;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate link percentage
%             if total_possible > 0
%                 link_percentages(pat, t) = (link_count / total_possible) * 100;
% 
%                 % Calculate probability of link and no link
%                 p_link = link_count / total_possible;
%                 p_no_link = 1 - p_link;
% 
%                 % Calculate Shannon entropy: -sum(p_i * log(p_i))
%                 % Only include non-zero probabilities in calculation
%                 entropy = 0;
%                 if p_link > 0
%                     entropy = entropy - p_link * log2(p_link);
%                 end
%                 if p_no_link > 0
%                     entropy = entropy - p_no_link * log2(p_no_link);
%                 end
% 
%                 link_entropy(pat, t) = entropy;
%             else
%                 link_percentages(pat, t) = 0;
%                 link_entropy(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data;
%     link_entropy = link_entropy(valid_indices, :);
%     link_percentages = link_percentages(valid_indices, :);
% 
%     % Initialize correlation results table
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanCorr', 'SpearmanP', 'PearsonCorr', 'PearsonP'});
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         entropy_values = link_entropy(:, t);
%         percentage_values = link_percentages(:, t);
% 
%         % Spearman correlation
%         [rho, pval] = corr(entropy_values, percentage_values, 'Type', 'Spearman');
% 
%         % Pearson correlation
%         [r, p] = corrcoef(entropy_values, percentage_values);
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanCorr(t) = rho;
%         corr_results.SpearmanP(t) = pval;
%         corr_results.PearsonCorr(t) = pearson_corr;
%         corr_results.PearsonP(t) = pearson_p;
%     end
% 
%     % Find the threshold with the highest absolute Spearman correlation
%     [~, max_spearman_idx] = max(abs(corr_results.SpearmanCorr));
%     max_spearman_threshold = thresholds(max_spearman_idx);
%     max_spearman_rho = corr_results.SpearmanCorr(max_spearman_idx);
%     max_spearman_p = corr_results.SpearmanP(max_spearman_idx);
% 
%     % Display correlation results table
%     disp('Correlation Results between Link Entropy and Link Percentage by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results: Link Entropy vs Link Percentage', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create scatter plots for each threshold
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Calculate subplot dimensions
%     subplot_rows = ceil(sqrt(num_thresholds));
%     subplot_cols = ceil(num_thresholds / subplot_rows);
% 
%     for t = 1:num_thresholds
%         subplot(subplot_rows, subplot_cols, t);
% 
%         entropy_values = link_entropy(:, t);
%         percentage_values = link_percentages(:, t);
% 
%         % INVERTED AXES: percentage on x-axis, entropy on y-axis
%         scatter(percentage_values, entropy_values, 70, 'filled', 'b');
%         hold on;
% 
%         % Add regression line
%         coef = polyfit(percentage_values, entropy_values, 1);
%         xfit = linspace(min(percentage_values), max(percentage_values), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%         % Calculate and display correlations
%         [r, p] = corrcoef(percentage_values, entropy_values);
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         [rho, pval] = corr(percentage_values, entropy_values, 'Type', 'Spearman');
% 
%         title(sprintf('Threshold = %.1f', thresholds(t)), 'FontSize', 12);
%         xlabel('Link Percentage (%)', 'FontSize', 10);  % Inverted axis label
%         ylabel('Link Entropy', 'FontSize', 10);  % Inverted axis label
% 
%         text(min(percentage_values) + 0.1*(max(percentage_values)-min(percentage_values)), ...
%              min(entropy_values) + 0.8*(max(entropy_values)-min(entropy_values)), ...
%              sprintf('r = %.2f (p = %.4f)\nrho = %.2f (p = %.4f)', ...
%                     pearson_corr, pearson_p, rho, pval), ...
%              'FontSize', 8);
% 
%         grid on;
%     end
% 
%     % Add a super title
%     sgtitle('Link Percentage vs Link Entropy by Threshold', 'FontSize', 16);
% 
%     % Create a theoretical curve figure
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Generate theoretical curve
%     p_values = linspace(0, 1, 1000);
%     entropy_values = zeros(size(p_values));
% 
%     for i = 1:length(p_values)
%         p = p_values(i);
%         q = 1 - p;
% 
%         if p > 0 && q > 0
%             entropy_values(i) = -p*log2(p) - q*log2(q);
%         end
%     end
% 
%     % Plot theoretical curve with INVERTED AXES
%     plot(p_values*100, entropy_values, 'k-', 'LineWidth', 2);
%     hold on;
% 
%     % Add actual data points for threshold = 0.2 with INVERTED AXES
%     scatter(link_percentages(:, 2), link_entropy(:, 2), 70, 'filled', 'r');
% 
%     title('Theoretical Relationship: Link Percentage vs Link Entropy', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);  % Inverted axis label
%     ylabel('Link Entropy', 'FontSize', 12);  % Inverted axis label
%     grid on;
%     legend({'Theoretical Curve', 'Actual Data (Threshold = 0.2)'}, 'Location', 'best');
% 
%     % Create a figure for the highest Spearman correlation
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Get values for the threshold with highest correlation
%     entropy_values = link_entropy(:, max_spearman_idx);
%     percentage_values = link_percentages(:, max_spearman_idx);
% 
%     % Create scatter plot with INVERTED AXES
%     scatter(percentage_values, entropy_values, 100, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     coef = polyfit(percentage_values, entropy_values, 1);
%     xfit = linspace(min(percentage_values), max(percentage_values), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 3);
% 
%     % Add correlation statistics
%     text(min(percentage_values) + 0.1*(max(percentage_values)-min(percentage_values)), ...
%          max(entropy_values) - 0.2*(max(entropy_values)-min(entropy_values)), ...
%          sprintf('Spearman rho = %.3f\np-value = %.4f', max_spearman_rho, max_spearman_p), ...
%          'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
% 
%     title('Highest Correlation: Link Percentage vs Link Entropy', 'FontSize', 16, 'FontWeight', 'bold');
%     xlabel('Link Percentage (%)', 'FontSize', 14);  % Inverted axis label
%     ylabel('Link Entropy', 'FontSize', 14);  % Inverted axis label
%     grid on;
% 
%     % Print summary statistics
%     fprintf('\nSummary Statistics:\n');
%     fprintf('Total Patients with Valid Data: %d\n', sum(valid_indices));
% 
%     % Print entropy statistics for threshold = 0.2
%     fprintf('\nLink Entropy Statistics (Threshold = 0.2):\n');
%     fprintf('Mean Link Entropy: %.3f \n', mean(link_entropy(:, 2)));
%     fprintf('Min Link Entropy: %.3f \n', min(link_entropy(:, 2)));
%     fprintf('Max Link Entropy: %.3f \n', max(link_entropy(:, 2)));
% 
%     % Print link percentage statistics for threshold = 0.2
%     fprintf('\nLink Percentage Statistics (Threshold = 0.2):\n');
%     fprintf('Mean Link Percentage: %.2f%%\n', mean(link_percentages(:, 2)));
%     fprintf('Min Link Percentage: %.2f%%\n', min(link_percentages(:, 2)));
%     fprintf('Max Link Percentage: %.2f%%\n', max(link_percentages(:, 2)));
% 
%     % Print correlation statistics for threshold = 0.2
%     fprintf('\nCorrelation between Link Percentage and Link Entropy (Threshold = 0.2):\n');
%     [r, p] = corrcoef(link_percentages(:, 2), link_entropy(:, 2));
%     fprintf('Pearson: r = %.3f (p = %.4f)\n', r(1,2), p(1,2));
% 
%     [rho, pval] = corr(link_percentages(:, 2), link_entropy(:, 2), 'Type', 'Spearman');
%     fprintf('Spearman: rho = %.3f (p = %.4f)\n', rho, pval);
% 
%     % Print highest correlation statistics
%     fprintf('\nHighest Correlation:\n');
%     fprintf('Threshold: %.1f\n', max_spearman_threshold);
%     fprintf('Spearman rho: %.3f (p = %.4f)\n', max_spearman_rho, max_spearman_p);
% end
% 
% vaginal_data = readtable("vaginal.xlsx");
% analyze_entropy_link_percentage_correlation(vaginal_data);


% function analyze_entropy_by_threshold(vaginal_data)
%     % This function analyzes the relationship between Shannon entropy, link percentages,
%     % and Nugent scores across different thresholds in vaginal microbiome data.
%     %
%     % Parameters:
%     % vaginal_data - Table containing vaginal microbiome data with SubjectID and NugentScore columns
% 
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     shannon_entropy = zeros(num_patients, num_thresholds);
%     link_percentages = zeros(num_patients, num_thresholds);
%     patient_mean_nugent = zeros(num_patients, 1);
%     valid_data = false(num_patients, 1);
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         % Skip patients with insufficient data
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Calculate total possible pairs
%             total_possible = (num_species * (num_species - 1))/2;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate link percentage
%             if total_possible > 0
%                 link_percentages(pat, t) = (link_count / total_possible) * 100;
% 
%                 % Calculate probability of link and no link
%                 p_link = link_count / total_possible;
%                 p_no_link = 1 - p_link;
% 
%                 % Calculate Shannon entropy: -sum(p_i * log(p_i))
%                 % Only include non-zero probabilities in calculation
%                 entropy = 0;
%                 if p_link > 0
%                     entropy = entropy - p_link * log2(p_link);
%                 end
%                 if p_no_link > 0
%                     entropy = entropy - p_no_link * log2(p_no_link);
%                 end
% 
%                 shannon_entropy(pat, t) = entropy;
%             else
%                 shannon_entropy(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data & ~isnan(patient_mean_nugent);
%     shannon_entropy = shannon_entropy(valid_indices, :);
%     link_percentages = link_percentages(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
% 
%     % Create a table for correlation results between Shannon entropy and Nugent scores
%     entropy_nugent_corr = table('Size', [num_thresholds, 4], ...
%                       'VariableTypes', {'double', 'double', 'double', 'double'}, ...
%                       'VariableNames', {'Threshold', 'SpearmanRho', 'SpearmanP', 'PearsonR'});
% 
%     % Create a table for correlation results between Shannon entropy and link percentages
%     entropy_link_corr = table('Size', [num_thresholds, 4], ...
%                       'VariableTypes', {'double', 'double', 'double', 'double'}, ...
%                       'VariableNames', {'Threshold', 'SpearmanRho', 'SpearmanP', 'PearsonR'});
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         % Entropy vs Nugent correlations
%         [rho, pval] = corr(shannon_entropy(:, t), patient_mean_nugent, 'Type', 'Spearman');
%         [r, ~] = corrcoef(shannon_entropy(:, t), patient_mean_nugent);
% 
%         entropy_nugent_corr.Threshold(t) = thresholds(t);
%         entropy_nugent_corr.SpearmanRho(t) = rho;
%         entropy_nugent_corr.SpearmanP(t) = pval;
%         entropy_nugent_corr.PearsonR(t) = r(1,2);
% 
%         % Entropy vs Link Percentage correlations
%         [rho, pval] = corr(shannon_entropy(:, t), link_percentages(:, t), 'Type', 'Spearman');
%         [r, ~] = corrcoef(shannon_entropy(:, t), link_percentages(:, t));
% 
%         entropy_link_corr.Threshold(t) = thresholds(t);
%         entropy_link_corr.SpearmanRho(t) = rho;
%         entropy_link_corr.SpearmanP(t) = pval;
%         entropy_link_corr.PearsonR(t) = r(1,2);
%     end
% 
%     % Display correlation tables
%     disp('Correlation between Shannon Entropy and Nugent Scores:');
%     disp(entropy_nugent_corr);
% 
%     disp('Correlation between Shannon Entropy and Link Percentages:');
%     disp(entropy_link_corr);
% 
%     % Create box plot for Shannon entropy by threshold
%     figure('Position', [100, 100, 1200, 600]);
%     boxplot(shannon_entropy, 'Labels', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false));
%     hold on;
% 
%     % Add mean markers to the boxplot
%     % mean_entropy = mean(shannon_entropy, 1, 'omitnan');
%     % plot(1:num_thresholds, mean_entropy, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
%     title('Shannon Entropy by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Shannon Entropy', 'FontSize', 12);
%     grid on;
% 
%     % Create correlation plots for Shannon entropy vs Nugent scores
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Find the threshold with the highest absolute Spearman correlation
%     [~, max_corr_idx] = max(abs(entropy_nugent_corr.SpearmanRho));
%     best_threshold = thresholds(max_corr_idx);
% 
%     % Scatter plot for the threshold with highest correlation
%     scatter(shannon_entropy(:, max_corr_idx), patient_mean_nugent, 100, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     [sorted_entropy, idx] = sort(shannon_entropy(:, max_corr_idx));
%     sorted_nugent = patient_mean_nugent(idx);
%     coef = polyfit(sorted_entropy, sorted_nugent, 1);
%     xfit = linspace(min(sorted_entropy), max(sorted_entropy), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Add correlation statistics
%     [rho, pval] = corr(shannon_entropy(:, max_corr_idx), patient_mean_nugent, 'Type', 'Spearman');
%     text(min(sorted_entropy) + 0.1*(max(sorted_entropy)-min(sorted_entropy)), ...
%          max(patient_mean_nugent) - 0.2*(max(patient_mean_nugent)-min(patient_mean_nugent)), ...
%          sprintf('Threshold = %.1f\nSpearman rho = %.3f\np-value = %.4f', best_threshold, rho, pval), ...
%          'FontSize', 12, 'BackgroundColor', [1 1 1 0.7]);
% 
%     title('Shannon Entropy vs Nugent Score (Best Correlation)', 'FontSize', 14);
%     xlabel('Shannon Entropy (bits)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Create correlation plot for Shannon entropy vs link percentages
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Find the threshold with the highest absolute Spearman correlation
%     [~, max_corr_idx] = max(abs(entropy_link_corr.SpearmanRho));
%     best_threshold = thresholds(max_corr_idx);
% 
%     % Scatter plot for the threshold with highest correlation
%     scatter(shannon_entropy(:, max_corr_idx), link_percentages(:, max_corr_idx), 100, 'filled', 'b');
%     hold on;
% 
%     % Add trend line
%     [sorted_entropy, idx] = sort(shannon_entropy(:, max_corr_idx));
%     sorted_links = link_percentages(idx, max_corr_idx);
%     coef = polyfit(sorted_entropy, sorted_links, 1);
%     xfit = linspace(min(sorted_entropy), max(sorted_entropy), 100);
%     yfit = polyval(coef, xfit);
%     plot(xfit, yfit, 'r-', 'LineWidth', 2);
% 
%     % Add correlation statistics
%     [rho, pval] = corr(shannon_entropy(:, max_corr_idx), link_percentages(:, max_corr_idx), 'Type', 'Spearman');
%     text(min(sorted_entropy) + 0.1*(max(sorted_entropy)-min(sorted_entropy)), ...
%          max(link_percentages(:, max_corr_idx)) - 0.2*(max(link_percentages(:, max_corr_idx))-min(link_percentages(:, max_corr_idx))), ...
%          sprintf('Threshold = %.1f\nSpearman rho = %.3f\np-value = %.4f', best_threshold, rho, pval), ...
%          'FontSize', 12, 'BackgroundColor', [1 1 1 0.7]);
% 
%     title('Shannon Entropy vs Link Percentage (Best Correlation)', 'FontSize', 14);
%     xlabel('Shannon Entropy (bits)', 'FontSize', 12);
%     ylabel('Link Percentage (%)', 'FontSize', 12);
%     grid on;
% 
%     % Print summary statistics
%     fprintf('\nSummary Statistics:\n');
%     fprintf('Total Patients with Valid Data: %d\n', sum(valid_indices));
% 
%     % Find the threshold with the highest correlation between entropy and Nugent
%     [max_corr, max_idx] = max(abs(entropy_nugent_corr.SpearmanRho));
%     fprintf('\nBest Threshold for Entropy-Nugent Correlation: %.1f (Spearman rho = %.3f, p = %.4f)\n', ...
%         thresholds(max_idx), entropy_nugent_corr.SpearmanRho(max_idx), entropy_nugent_corr.SpearmanP(max_idx));
% 
%     % Find the threshold with the highest correlation between entropy and link percentage
%     [max_corr, max_idx] = max(abs(entropy_link_corr.SpearmanRho));
%     fprintf('Best Threshold for Entropy-Link Correlation: %.1f (Spearman rho = %.3f, p = %.4f)\n', ...
%         thresholds(max_idx), entropy_link_corr.SpearmanRho(max_idx), entropy_link_corr.SpearmanP(max_idx));
% end
% 
% analyze_entropy_by_threshold(vaginal_data);

%% 8 heatmap link and no-link probabilities sorted alphabetically 


% function analyze_species_pair_probabilities(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.3;      % Threshold for score
%     thres_R = 0.9;      % Threshold for region coverage
%     time_interval = 1;
% 
%     % Get species names from column headers
%     species_names = vaginal_data.Properties.VariableNames(11:end);
%     num_species_total = length(species_names);
% 
%     fprintf('Analyzing %d patients and %d species...\n', num_patients, num_species_total);
% 
%     % Create arrays to store pair data
%     % Each row will represent one pair occurrence across all patients
%     all_pairs = [];  % Will store [species1_idx, species2_idx, patient_idx, has_link]
% 
%     % For each patient
%     valid_patient_count = 0;
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         % Get indices of species present in this patient (>1% abundance)
%         present_species_idx = find(keep_species);
%         num_species = length(present_species_idx);
% 
%         if num_species < 2
%             continue;
%         end
% 
%         valid_patient_count = valid_patient_count + 1;
% 
%         % For each pair of present species
%         for i = 1:num_species
%             for j = (i+1):num_species
%                 % Get the global indices of these species
%                 global_i = present_species_idx(i);
%                 global_j = present_species_idx(j);
% 
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Determine if this pair has a link (1) or not (0)
%                 if has_link_ij || has_link_ji
%                     has_link = 1;
%                 else
%                     has_link = 0;
%                 end
% 
%                 % Store this pair occurrence
%                 all_pairs = [all_pairs; global_i, global_j, pat, has_link];
%             end
%         end
%     end
% 
%     fprintf('Found %d valid patients with at least 2 species at >1%% abundance\n', valid_patient_count);
%     fprintf('Recorded %d pair occurrences across all patients\n', size(all_pairs, 1));
% 
%     % Group by species pairs and calculate link probabilities
%     unique_pairs = unique(all_pairs(:, 1:2), 'rows');
%     num_unique_pairs = size(unique_pairs, 1);
% 
%     % Initialize arrays for results
%     pair_names = cell(num_unique_pairs, 1);
%     times_present = zeros(num_unique_pairs, 1);
%     times_linked = zeros(num_unique_pairs, 1);
%     link_probabilities = zeros(num_unique_pairs, 1);
%     no_link_probabilities = zeros(num_unique_pairs, 1);
% 
%     % For each unique pair
%     for p = 1:num_unique_pairs
%         sp1 = unique_pairs(p, 1);
%         sp2 = unique_pairs(p, 2);
% 
%         % Find all occurrences of this pair
%         pair_occurrences = all_pairs(:, 1) == sp1 & all_pairs(:, 2) == sp2;
% 
%         % Count occurrences and links
%         times_present(p) = sum(pair_occurrences);
%         times_linked(p) = sum(all_pairs(pair_occurrences, 4));
% 
%         % Calculate link probability and no-link probability
%         link_probabilities(p) = times_linked(p) / times_present(p);
%         no_link_probabilities(p) = 1 - link_probabilities(p);
% 
%         % Create pair name
%         pair_names{p} = sprintf('%s - %s', species_names{sp1}, species_names{sp2});
%     end
% 
%     % Create a table
%     pair_table = table(pair_names, times_present, times_linked, link_probabilities, no_link_probabilities, ...
%                       'VariableNames', {'SpeciesPair', 'TimesPresent', 'TimesLinked', 'LinkProbability', 'NoLinkProbability'});
% 
%     % Sort by link probability
%     pair_table = sortrows(pair_table, 'LinkProbability', 'descend');
% 
%     % Display the table
%     disp('Species Pairs and Their Link Probabilities:');
%     disp(pair_table);
% 
%     % Create probability matrices for easier visualization
%     % First, identify all unique species that appear in pairs
%     all_species_indices = unique([unique_pairs(:,1); unique_pairs(:,2)]);
%     num_species_in_pairs = length(all_species_indices);
% 
%     % Get species names for matrix labels and clean them up
%     matrix_species_names = cell(num_species_in_pairs, 1);
%     for i = 1:num_species_in_pairs
%         species_name = species_names{all_species_indices(i)};
% 
%         % Clean up species name for display
%         clean_name = species_name;
% 
%         % Replace underscores with spaces
%         clean_name = strrep(clean_name, '_', ' ');
% 
%         % Fix common genus abbreviations
%         if startsWith(clean_name, 'L ')
%             clean_name = ['Lactobacillus', clean_name(2:end)];
%         end
% 
%         matrix_species_names{i} = clean_name;
%     end
% 
%     % Sort species alphabetically
%     [sorted_species_names, sort_idx] = sort(matrix_species_names);
%     sorted_species_indices = all_species_indices(sort_idx);
% 
%     % Create the link probability matrix with alphabetically sorted species
%     link_prob_matrix = zeros(num_species_in_pairs, num_species_in_pairs);
%     no_link_prob_matrix = zeros(num_species_in_pairs, num_species_in_pairs);
%     presence_count_matrix = zeros(num_species_in_pairs, num_species_in_pairs);
% 
%     % Map global species indices to matrix indices for faster lookup
%     species_index_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
%     for i = 1:num_species_in_pairs
%         species_index_map(sorted_species_indices(i)) = i;
%     end
% 
%     % Fill the matrices
%     for p = 1:num_unique_pairs
%         sp1 = unique_pairs(p, 1);
%         sp2 = unique_pairs(p, 2);
% 
%         % Check if both species are in our map (they should be)
%         if isKey(species_index_map, sp1) && isKey(species_index_map, sp2)
%             % Get matrix indices
%             mat_i = species_index_map(sp1);
%             mat_j = species_index_map(sp2);
% 
%             % Fill both directions (symmetric matrix)
%             link_prob_matrix(mat_i, mat_j) = link_probabilities(p);
%             link_prob_matrix(mat_j, mat_i) = link_probabilities(p);
% 
%             no_link_prob_matrix(mat_i, mat_j) = no_link_probabilities(p);
%             no_link_prob_matrix(mat_j, mat_i) = no_link_probabilities(p);
% 
%             presence_count_matrix(mat_i, mat_j) = times_present(p);
%             presence_count_matrix(mat_j, mat_i) = times_present(p);
%         end
%     end
% 
%     % Create heatmap of link probabilities with alphabetically sorted species
%     figure('Position', [100, 100, 1200, 1000]);
%     h = heatmap(sorted_species_names, sorted_species_names, link_prob_matrix);
%     h.Title = 'Link Probability Matrix';
%     h.XLabel = '';
%     h.YLabel = '';
%     h.ColorbarVisible = 'on';
% 
%     % Create heatmap of no-link probabilities with alphabetically sorted species
%     figure('Position', [100, 100, 1200, 1000]);
%     h = heatmap(sorted_species_names, sorted_species_names, no_link_prob_matrix);
%     h.Title = 'No-Link Probability Matrix';
%     h.XLabel = '';
%     h.YLabel = '';
%     h.ColorbarVisible = 'on';
% 
%     % Create heatmap of presence counts with alphabetically sorted species
%     figure('Position', [100, 100, 1200, 1000]);
%     h = heatmap(sorted_species_names, sorted_species_names, presence_count_matrix);
%     h.Title = 'Pair Presence Count Matrix';
%     h.XLabel = '';
%     h.YLabel = '';
%     h.ColorbarVisible = 'on';
% 
%     % Return the matrices and species names for further analysis
%     link_prob_result.link_prob_matrix = link_prob_matrix;
%     link_prob_result.no_link_prob_matrix = no_link_prob_matrix;
%     link_prob_result.presence_count_matrix = presence_count_matrix;
%     link_prob_result.species_names = sorted_species_names;
%     link_prob_result.pair_table = pair_table;
% 
%     % Save the results to a MAT file for later use
%     save('species_pair_probabilities.mat', 'link_prob_result');
%     disp('Results saved to species_pair_probabilities.mat');
% end
% 
% analyze_species_pair_probabilities(vaginal_data);
 



%% 9 single species level link percent analysis 
% function analyze_species_link_percentages(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_S = 0.3;      % Threshold for score
%     thres_R = 0.9;      % Threshold for region coverage
%     time_interval = 1;
% 
%     % Get species names from column headers
%     species_names = vaginal_data.Properties.VariableNames(11:end);
%     num_species_total = length(species_names);
% 
%     fprintf('Analyzing %d patients and %d species...\n', num_patients, num_species_total);
% 
%     % Create arrays to track individual species statistics
%     % For each species, track how many possible connections it could have and how many it actually forms
%     species_stats = struct();
% 
%     % For each patient
%     valid_patient_count = 0;
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         % Get indices of species present in this patient (>1% abundance)
%         present_species_idx = find(keep_species);
%         num_species = length(present_species_idx);
% 
%         if num_species < 2
%             continue;
%         end
% 
%         valid_patient_count = valid_patient_count + 1;
% 
%         % For each species in this patient
%         for i = 1:num_species
%             global_i = present_species_idx(i);
%             species_name = species_names{global_i};
% 
%             % Initialize this species in our tracking structure if needed
%             if ~isfield(species_stats, species_name)
%                 species_stats.(species_name) = struct('possible_links', 0, 'actual_links', 0, 'patient_link_percentages', []);
%             end
% 
%             % This species could potentially link with all other species in this patient
%             possible_links_in_patient = num_species - 1;
%             actual_links_in_patient = 0;
% 
%             % Check links with all other species
%             for j = 1:num_species
%                 if i == j
%                     continue; % Skip self
%                 end
% 
%                 global_j = present_species_idx(j);
% 
%                 % Check i->j direction
%                 [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Check j->i direction
%                 [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                     (1:size(rel_abun_1pct,1))', time_interval);
%                 has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                 % Count this link if it exists in either direction
%                 if has_link_ij || has_link_ji
%                     actual_links_in_patient = actual_links_in_patient + 1;
%                 end
%             end
% 
%             % Update species statistics
%             species_stats.(species_name).possible_links = species_stats.(species_name).possible_links + possible_links_in_patient;
%             species_stats.(species_name).actual_links = species_stats.(species_name).actual_links + actual_links_in_patient;
% 
%             % Calculate and store link percentage for this species in this patient
%             if possible_links_in_patient > 0
%                 link_percentage = (actual_links_in_patient / possible_links_in_patient) * 100;
%                 species_stats.(species_name).patient_link_percentages = [species_stats.(species_name).patient_link_percentages, link_percentage];
%             end
%         end
%     end
% 
%     fprintf('Found %d valid patients with at least 2 species at >1%% abundance\n', valid_patient_count);
% 
%     % Create a table to display species-level statistics
%     all_species = fieldnames(species_stats);
%     num_species = length(all_species);
% 
%     species_names_table = cell(num_species, 1);
%     possible_links_total = zeros(num_species, 1);
%     actual_links_total = zeros(num_species, 1);
%     overall_link_percentages = zeros(num_species, 1);
%     avg_patient_link_percentages = zeros(num_species, 1);
%     patient_appearances = zeros(num_species, 1);
% 
%     for s = 1:num_species
%         species_name = all_species{s};
%         species_data = species_stats.(species_name);
% 
%         % Clean up species name for display
%         clean_name = species_name;
% 
%         % Replace underscores with spaces
%         clean_name = strrep(clean_name, '_', ' ');
% 
%         % Fix common genus abbreviations
%         if startsWith(clean_name, 'L ')
%             clean_name = ['Lactobacillus', clean_name(2:end)];
%         end
% 
%         % Store the cleaned name
%         species_names_table{s} = clean_name;
% 
%         possible_links_total(s) = species_data.possible_links;
%         actual_links_total(s) = species_data.actual_links;
% 
%         % Calculate overall link percentage (total links / total possible)
%         if possible_links_total(s) > 0
%             overall_link_percentages(s) = (actual_links_total(s) / possible_links_total(s)) * 100;
%         end
% 
%         % Calculate average of per-patient link percentages
%         if ~isempty(species_data.patient_link_percentages)
%             avg_patient_link_percentages(s) = mean(species_data.patient_link_percentages);
%             patient_appearances(s) = length(species_data.patient_link_percentages);
%         end
%     end
% 
%     % Create and display the table
%     species_table = table(species_names_table, patient_appearances, possible_links_total, actual_links_total, ...
%                          overall_link_percentages, avg_patient_link_percentages, ...
%                          'VariableNames', {'Species', 'PatientAppearances', 'PossibleLinks', 'ActualLinks', ...
%                                           'OverallLinkPercentage', 'AvgPatientLinkPercentage'});
% 
%     % Filter out species with 0% link percentage
%     species_table = species_table(species_table.OverallLinkPercentage > 0, :);
% 
%     % Sort by overall link percentage (descending)
%     species_table = sortrows(species_table, 'OverallLinkPercentage', 'descend');
% 
%     disp('Species-Level Link Statistics (excluding 0% links):');
%     disp(species_table);
% 
%     % Create horizontal bar chart of overall link percentages for non-zero species
%     figure('Position', [100, 100, 1000, 1200]);  % Taller figure for horizontal bars
% 
%     % Get number of species to plot
%     num_species_to_plot = height(species_table);
% 
%     % Create horizontal bar chart
%     barh(flip(species_table.OverallLinkPercentage));  % flip to have highest at top
%     title('Overall Link Percentage by Species', 'FontSize', 14);
%     ylabel('Species', 'FontSize', 12);
%     xlabel('Percentage of Possible Links Formed (%)', 'FontSize', 12);
% 
%     % Set y-ticks and labels (flipped order to have highest at top)
%     yticks(1:num_species_to_plot);
%     yticklabels(flip(species_table.Species));
% 
%     % Set x-axis from 0 to 100
%     xlim([0 100]);
% 
%     % Remove grid
%     grid off;
% 
%     % Save the results to a MAT file for later use
%     save('species_link_percentages.mat', 'species_table', 'species_stats');
%     disp('Results saved to species_link_percentages.mat');
% end
% 
% analyze_species_link_percentages(vaginal_data);






% %% 3rd task: pairs percent correlated with nugent score
% function analyze_link_percent_sensitivity(vaginal_data)
%     % Get unique subject IDs
%     subject_ids = unique(vaginal_data.SubjectID);
%     num_patients = length(subject_ids);
% 
%     % Parameters for RDS
%     thres_noise = 0;
%     thres_R = 0.9;  % Threshold for region coverage
%     time_interval = 1;
% 
%     % Define thresholds for sensitivity analysis
%     thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%     num_thresholds = length(thresholds);
% 
%     % Initialize arrays for storing results
%     all_link_percentages = zeros(num_patients, num_thresholds);
%     total_possible_pairs = zeros(num_patients, 1);
%     patient_mean_nugent = zeros(1, num_patients);
%     patient_median_nugent = zeros(1, num_patients);
%     valid_data = false(1, num_patients);  % Track which patients have valid data
% 
%     % Initialize correlation results table
%     corr_results = table('Size', [num_thresholds, 5], ...
%                          'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
%                          'VariableNames', {'Threshold', 'SpearmanMean', 'SpearmanMeanP', 'PearsonMean', 'PearsonMeanP'});
% 
%     % For each patient
%     for pat = 1:num_patients
%         subject_id = subject_ids(pat);
%         rows = vaginal_data.SubjectID == subject_id;
%         abun = table2array(vaginal_data(rows, 11:end));
% 
%         % Store Nugent scores
%         patient_mean_nugent(pat) = mean(vaginal_data.NugentScore(rows), 'omitnan');
%         patient_median_nugent(pat) = median(vaginal_data.NugentScore(rows), 'omitnan');
% 
%         if size(abun, 1) < 3
%             continue;
%         end
% 
%         % Process abundances
%         rel_abun = abun ./ sum(abun, 2);
%         mean_rel = mean(rel_abun, 1);
%         keep_species = mean_rel > 0.01;
%         rel_abun_1pct = rel_abun(:, keep_species);
% 
%         num_species = size(rel_abun_1pct, 2);
%         if num_species < 2
%             continue;
%         end
% 
%         % Mark this patient as having valid data
%         valid_data(pat) = true;
% 
%         % Calculate total possible pairs
%         total_possible = (num_species * (num_species - 1))/2;
%         total_possible_pairs(pat) = total_possible;
% 
%         % For each threshold
%         for t = 1:num_thresholds
%             thres_S = thresholds(t);
%             link_count = 0;
% 
%             % Check for links between species pairs
%             for i = 1:num_species
%                 for j = (i+1):num_species
%                     % Check i->j direction
%                     [score_list_ij, ~, ~] = RDS_dim1(rel_abun_1pct(:,i), rel_abun_1pct(:,j), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ij = check_significance(score_list_ij, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % Check j->i direction
%                     [score_list_ji, ~, ~] = RDS_dim1(rel_abun_1pct(:,j), rel_abun_1pct(:,i), ...
%                         (1:size(rel_abun_1pct,1))', time_interval);
%                     has_link_ji = check_significance(score_list_ji, thres_noise, thres_S, thres_R, size(rel_abun_1pct,1));
% 
%                     % If either direction has a link, count it
%                     if has_link_ij || has_link_ji
%                         link_count = link_count + 1;
%                     end
%                 end
%             end
% 
%             % Calculate and store percentage
%             if total_possible > 0
%                 all_link_percentages(pat, t) = (link_count / total_possible) * 100;
%             else
%                 all_link_percentages(pat, t) = 0;
%             end
%         end
%     end
% 
%     % Remove any patients with no data
%     valid_indices = valid_data & ~isnan(patient_mean_nugent) & ~isnan(patient_median_nugent);
%     all_link_percentages = all_link_percentages(valid_indices, :);
%     patient_mean_nugent = patient_mean_nugent(valid_indices);
%     patient_median_nugent = patient_median_nugent(valid_indices);
%     subject_ids = subject_ids(valid_indices);
% 
%     % Calculate correlations for each threshold
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
% 
%         % Spearman correlation with mean Nugent
%         [rho, pval] = corr(link_percentages, patient_mean_nugent', 'Type', 'Spearman');
% 
%         % Pearson correlation with mean Nugent
%         [r, p] = corrcoef(link_percentages, patient_mean_nugent');
%         pearson_corr = r(1,2);
%         pearson_p = p(1,2);
% 
%         % Store in results table
%         corr_results.Threshold(t) = thresholds(t);
%         corr_results.SpearmanMean(t) = rho;
%         corr_results.SpearmanMeanP(t) = pval;
%         corr_results.PearsonMean(t) = pearson_corr;
%         corr_results.PearsonMeanP(t) = pearson_p;
%     end
% 
%     % Print summary statistics for default threshold (0.2)
%     fprintf('\nSummary Statistics (Threshold = 0.2):\n');
%     fprintf('Total Patients: %d\n', sum(valid_indices));
%     fprintf('Mean Link Percentage: %.2f%%\n', mean(all_link_percentages(:, 2), 'omitnan'));
%     fprintf('Mean Nugent Score: %.2f\n', mean(patient_mean_nugent, 'omitnan'));
% 
%     % Display correlation results table
%     disp('Correlation Results by Threshold:');
%     disp(corr_results);
% 
%     % Create a figure for the correlation table
%     figure('Position', [100, 100, 800, 400]);
% 
%     % Create a uitable
%     t = uitable('Data', table2array(corr_results), ...
%                 'ColumnName', corr_results.Properties.VariableNames, ...
%                 'RowName', arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false), ...
%                 'Units', 'normalized', ...
%                 'Position', [0.05, 0.05, 0.9, 0.9]);
% 
%     % Set title
%     annotation('textbox', [0.5, 0.95, 0, 0], 'String', 'Correlation Results by Threshold (Link Percentage)', ...
%                'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
% 
%     % Create 4x4 correlation plots with multiple regression lines
%     figure('Position', [100, 100, 1200, 900]);
% 
%     % Define colors for different thresholds
%     colors = jet(num_thresholds);
% 
%     % 1. Mean Nugent vs Link Percentage (Spearman)
%     subplot(2, 2, 1);
%     hold on;
% 
%     % Plot scatter for default threshold (0.2)
%     default_idx = 2; % 0.2 is the second threshold
%     scatter(all_link_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         [sorted_links, idx] = sort(link_percentages);
%         sorted_nugent = patient_mean_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Mean Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 2. Mean Nugent vs Link Percentage (Pearson)
%     subplot(2, 2, 2);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_mean_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         coef = polyfit(link_percentages, patient_mean_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Mean Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Mean Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 3. Median Nugent vs Link Percentage (Spearman)
%     subplot(2, 2, 3);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add trend lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         [sorted_links, idx] = sort(link_percentages);
%         sorted_nugent = patient_median_nugent(idx);
%         coef = polyfit(sorted_links, sorted_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Median Nugent (Spearman)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % 4. Median Nugent vs Link Percentage (Pearson)
%     subplot(2, 2, 4);
%     hold on;
% 
%     % Plot scatter for default threshold
%     scatter(all_link_percentages(:, default_idx), patient_median_nugent, 70, 'filled', 'MarkerFaceColor', colors(default_idx,:));
% 
%     % Add regression lines for all thresholds
%     for t = 1:num_thresholds
%         link_percentages = all_link_percentages(:, t);
%         coef = polyfit(link_percentages, patient_median_nugent, 1);
%         xfit = linspace(min(link_percentages), max(link_percentages), 100);
%         yfit = polyval(coef, xfit);
%         plot(xfit, yfit, 'Color', colors(t,:), 'LineWidth', 2);
%     end
% 
%     title('Link Percentage vs Median Nugent (Pearson)', 'FontSize', 14);
%     xlabel('Link Percentage (%)', 'FontSize', 12);
%     ylabel('Median Nugent Score', 'FontSize', 12);
%     grid on;
% 
%     % Add a common colorbar for all subplots
%     colormap(jet);
%     cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
%     cb.Ticks = linspace(0, 1, num_thresholds);
%     cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), thresholds, 'UniformOutput', false);
%     ylabel(cb, 'Threshold Value', 'FontSize', 12);
% 
%     % Create a figure showing average link percentages by threshold
%     figure('Position', [100, 100, 800, 600]);
% 
%     % Calculate mean link percentages for each threshold
%     mean_link_percentages = mean(all_link_percentages, 1);
% 
%     % Create bar chart
%     bar(thresholds, mean_link_percentages);
%     title('Mean Link Percentage by Threshold', 'FontSize', 14);
%     xlabel('Threshold Value', 'FontSize', 12);
%     ylabel('Mean Link Percentage (%)', 'FontSize', 12);
%     grid on;
% 
%     % Add text labels with percentages
%     for t = 1:num_thresholds
%         text(thresholds(t), mean_link_percentages(t) + 2, sprintf('%.1f%%', mean_link_percentages(t)), ...
%              'HorizontalAlignment', 'center', 'FontSize', 10);
%     end
% 
%     % Create a figure showing link percentages by Nugent category
%     figure('Position', [100, 100, 1200, 600]);
% 
%     % Define Nugent score categories
%     low_nugent = [0, 3];    % Healthy
%     mid_nugent = [4, 6];    % Intermediate
%     high_nugent = [7, 10];  % BV
% 
%     % Categorize patients
%     low_indices = patient_mean_nugent >= low_nugent(1) & patient_mean_nugent <= low_nugent(2);
%     mid_indices = patient_mean_nugent >= mid_nugent(1) & patient_mean_nugent <= mid_nugent(2);
%     high_indices = patient_mean_nugent >= high_nugent(1) & patient_mean_nugent <= high_nugent(2);
% 
%     % Calculate means for each category and threshold
%     nugent_categories = {'Low (0-3)', 'Mid (4-6)', 'High (7-10)'};
%     mean_by_category = zeros(3, num_thresholds);
% 
%     for t = 1:num_thresholds
%         mean_by_category(1, t) = mean(all_link_percentages(low_indices, t), 'omitnan');
%         mean_by_category(2, t) = mean(all_link_percentages(mid_indices, t), 'omitnan');
%         mean_by_category(3, t) = mean(all_link_percentages(high_indices, t), 'omitnan');
%     end
% 
%     % Create grouped bar chart
%     bar(mean_by_category);
%     title('Mean Link Percentage by Nugent Category and Threshold', 'FontSize', 14);
%     xlabel('Nugent Category', 'FontSize', 12);
%     ylabel('Mean Link Percentage (%)', 'FontSize', 12);
%     set(gca, 'XTick', 1:3, 'XTickLabel', nugent_categories);
%     legend(arrayfun(@(x) sprintf('Threshold = %.1f', x), thresholds, 'UniformOutput', false), 'Location', 'best');
%     grid on;
% 
%     % Add patient counts
%     text(1, max(mean_by_category(1,:)) + 5, sprintf('n=%d', sum(low_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
%     text(2, max(mean_by_category(2,:)) + 5, sprintf('n=%d', sum(mid_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
%     text(3, max(mean_by_category(3,:)) + 5, sprintf('n=%d', sum(high_indices)), 'HorizontalAlignment', 'center', 'FontSize', 12);
% end


%% 10 classify with the nugent 
% pending 


