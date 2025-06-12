function is_significant = check_significance(score_list, thres_noise, thres_S, thres_R, T)
    is_significant = false;
    
    for type = 1:2
        score = reshape(score_list(:,:,type), [T, T]);
        
        %noise threshold 
        local = find(abs(score) > thres_noise);

        if ~isempty(local)
            s = sum(score(local)) / sum(abs(score(local)));
            r = length(local) / (T * T / 2);

            %satisfies threshold and region 
            if abs(s) >= thres_S && r >= thres_R
                is_significant = true;
                return;
            end
        end
    end
end
