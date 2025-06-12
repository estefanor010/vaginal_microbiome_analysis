function [parameters] = parameters()
    parameters.plot_width = 700;
    parameters.plot_height = 700;

    parameters.thres_noise = 0;
    parameters.thres_S = 0.3;
    parameters.thres_R = 0.9;
    parameters.time_interval = 1;
    parameters.thresholds = 0.1:0.1:0.9;
    
    parameters.abundance_threshold = 0.01;
    parameters.filter_method = 'max';

    
    
end