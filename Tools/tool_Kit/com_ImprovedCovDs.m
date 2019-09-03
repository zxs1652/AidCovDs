%% Genetating the Riemannian CovDs
% block_Matrix,option,type_Block
% Author: Kai-Xuan Chen
% Date(1): 2018.06.08
% Date(2): 2018.09.28

function final_ImprovedCovDs = com_ImprovedCovDs(block_Matrix,option,type_Block)
    
    [ImprovedCovDs, ~] = com_ImprovedCovDs_Arccos(block_Matrix,option,type_Block);  
    switch option.kernel_Type
        case 'LogEk'
            [dim,~,~] = size(block_Matrix);
            ImprovedCovDs = ImprovedCovDs/((dim + 1)^2 - 1);
        case 'LogEk.Arc'
            ImprovedCovDs = ImprovedCovDs;
    end

%     [~,S] = eig(CSPD_Matrix);
%     [temp_min,~] = min(diag(S));
%     Tr = trace(Improve_CSPD);
    if option.if_JudgeiCovDsEig
        [~,S] = eig(ImprovedCovDs);
        [temp_min,~] = min(diag(S));
%         while temp_min <= option.min_Eig_iCovDs
        while temp_min <= option.min_Eig_iCovDs
            ImprovedCovDs = ImprovedCovDs + (option.min_Eig_iCovDs)*trace(ImprovedCovDs)*eye(size(ImprovedCovDs));
%             ImprovedCovDs = ImprovedCovDs + tem_Perturbation*trace(ImprovedCovDs)*eye(size(ImprovedCovDs));
            [~,S] = eig(ImprovedCovDs);
            [temp_min,~] = min(diag(S));
        end
    else
        ImprovedCovDs = ImprovedCovDs + (option.min_Eig_iCovDs)*trace(ImprovedCovDs)*eye(size(ImprovedCovDs));   
%         ImprovedCovDs = ImprovedCovDs + tem_Perturbation*trace(ImprovedCovDs)*eye(size(ImprovedCovDs));      
    end
    final_ImprovedCovDs = ImprovedCovDs;
end