function has_link = calculate_links(X, Y, parameters)
    %Input
        % X: temporal variable
        % Y: temporal variable 
        % parameters: 
    %Output 
        % has_link (boolean): yes/no link 
    
    
    
    [score_list_ij, ~, ~] = RDS_dim1(X, Y, (1:size(X,1))', parameters.time_interval);
    has_link_ij = check_significance(score_list_ij, parameters.thres_noise, parameters.thres_S, parameters.thres_R, size(X,1));
    
    [score_list_ji, ~, ~] = RDS_dim1(Y, X, (1:size(X,1))', parameters.time_interval);
    has_link_ji = check_significance(score_list_ji, parameters.thres_noise, parameters.thres_S, parameters.thres_R, size(X,1));
    
    has_link = has_link_ij || has_link_ji;
end
