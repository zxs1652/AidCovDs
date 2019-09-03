function out_Spd = fea2SPD(feature_Matrix)
%     tmp_FeatureMatrix = zscore(feature_Matrix);
    tmp_FeatureMatrix = feature_Matrix;
    temp_Cov = cov(tmp_FeatureMatrix);
    [~,S] = eig(temp_Cov);
    [temp_min,~] = min(diag(S));
    if temp_min <= 10^(-9)
        tmp_FeatureMatrix = feature_Matrix + rand(size(feature_Matrix));
        temp_Cov = cov(tmp_FeatureMatrix);
        [~,S] = eig(temp_Cov);
        [temp_min,~] = min(diag(S));
        while temp_min <= 10^(-9)
            temp_Cov = temp_Cov + 0.001*trace(temp_Cov)*eye(size(temp_Cov));
            [~,S] = eig(temp_Cov);
            [temp_min,~] = min(diag(S));
        end  
    elseif norm(temp_Cov,'fro') == 0 && trace(temp_Cov) == 0
        temp_Cov = size(temp_Cov);
    end 
    out_Spd = temp_Cov+0.001*trace(temp_Cov)*eye(size(temp_Cov));
end