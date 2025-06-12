function plot_link_correlations(results, parameters)
   
    figure('Position', [100, 100, parameters.plot_width*2, parameters.plot_height*2]);
    
    % Mean Nugent vs Link % using SPEARMAN 
    subplot(2, 2, 1);
    create_correlation_plot(results.link_percent, results.mean_nugent, 'Spearman', 'Mean Nugent Score');
    title('Link % vs Mean Nugent (Spearman)', 'FontSize', 14);
    
    % Mean Nugent vs Link % using PEARSON 
    subplot(2, 2, 2);
    create_correlation_plot(results.link_percent, results.mean_nugent, 'Pearson', 'Mean Nugent Score');
    title('Link % vs Mean Nugent (Pearson)', 'FontSize', 14);
    
    % Median Nugent vs Link % using SPEARMAN
    subplot(2, 2, 3);
    create_correlation_plot(results.link_percent, results.median_nugent,'Spearman', 'Median Nugent Score');
    title('Link % vs Median Nugent (Spearman)', 'FontSize', 14);
    
    % Median Nugent vs Link % using PEARSON 
    subplot(2, 2, 4);
    create_correlation_plot(results.link_percent, results.median_nugent,'Pearson', 'Median Nugent Score');
    title('Link % vs Median Nugent (Pearson)', 'FontSize', 14);

    % Shared plot creation function
    function create_correlation_plot(input, output, corr_type, ylabel_text)
        scatter(input, output, 70, 'filled', 'b');
        hold on;
        
        % linear fit logic 
        [sorted_input, idx] = sort(input);
        sorted_output = output(idx);
        coef = polyfit(sorted_input, sorted_output, 1); %linear coefficient 
        input_fit = linspace(min(input), max(input), 100); 
        output_fit = polyval(coef, input_fit); %plot linear trend 
        plot(input_fit, output_fit, 'r-', 'LineWidth', 2); %plot input and output 
        
        % correlations computation 
        switch lower(corr_type)
            case 'spearman'
                [rho, pval] = corr(input, output, 'Type', 'Spearman');
                correlation_report = sprintf('Spearman \\rho = %.2f\np = %.4f', rho, pval);
            case 'pearson'
                [r, p] = corrcoef(input, output);
                rho = r(1,2);  % Now using rho/pval for Pearson too
                pval = p(1,2);
                correlation_report = sprintf('Pearson r = %.2f\np = %.4f', rho, pval);
        end 
       
        %correlation report text location  
        text(min(input) + 0.1*(max(input)-min(input)), max(ouput) - 0.1*(max(output)-min(output)), ...
            correlation_report, 'FontSize', 12, 'BackgroundColor', [1 1 1 0.8]);

        %labels 
        grid off;
        xlabel('Link Percentage (%)', 'FontSize', 12);
        ylabel(ylabel_text, 'FontSize', 12);
        set(gca, 'FontSize', 12); 
    end
end
