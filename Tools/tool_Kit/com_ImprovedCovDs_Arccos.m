function [iCovDs, salient_Matrix] = com_ImprovedCovDs_Arccos(block_Matrix,option,type_Block)
    [dim,~,num_Block] = size(block_Matrix);
    dim_LogDomain = (dim + 1)^2;
    iCovDs = zeros(num_Block,num_Block);
    salient_Matrix = zeros(num_Block,num_Block);
    for block_th = 1:num_Block
        current_Block = block_Matrix(:,:,block_th);
    %                 current_Block = current_Block/255;  
    %                 current_Block = current_Block;  
        %%     embed  Gauss into SPD        
        current_Cov = embed_Gauss2SPD_Sim(current_Block,option);     
        log_current_Cov = logm(current_Cov);
%         log_current_Cov = log_current_Cov - mean2(log_current_Cov).*ones(size(log_current_Cov));
        mean_C = mean(log_current_Cov,2);
        mean_R = mean(log_current_Cov,1);
        log_current_Cov = log_current_Cov + (1/(size(log_current_Cov,1)*size(log_current_Cov,1)))*mean2(log_current_Cov).*ones(size(log_current_Cov))...
            - 1/(size(log_current_Cov,1))*(repmat(mean_C,[1,size(log_current_Cov,2)]) + repmat(mean_R,[size(log_current_Cov,1),1]));
        log_Matrix(:,:,block_th) = log_current_Cov;
    end    
    switch option.kernel_Type
        case 'LogEk'   % first release    
            vectors = reshape(log_Matrix, size(log_Matrix, 1) * size(log_Matrix, 2), size(log_Matrix, 3))';
            iCovDs= vectors*vectors';  
            vectors = bsxfun(@rdivide, vectors, sqrt(sum(vectors.^2, 2)));
            salient_Matrix = vectors*vectors';         
        case 'LogEk.Arc'  % second release 
            vectors = reshape(log_Matrix, size(log_Matrix, 1) * size(log_Matrix, 2), size(log_Matrix, 3));
            vectors_Fro = sqrt(sum(vectors.*vectors,1)); 
        %     vectors = vectors./repmat(vectors_Fro,[size(vectors,1),1]);
        %     vectors_Fro = vectors_Fro./vectors_Fro;
        %     vectors_Norm = vectors./repmat(vectors_Fro,[size(vectors,1),1]);
            inner_Matrix = vectors_Fro'*vectors_Fro;
        %     inner_Matrix = vectors'*vectors;
            jtheta_Matrix0 = zeros(num_Block,num_Block); jtheta_Matrix1 = zeros(num_Block,num_Block);
            jtheta_Matrix2 = zeros(num_Block,num_Block); jtheta_Matrix3 = zeros(num_Block,num_Block);
            for i_th = 1:num_Block
                vector_I = vectors(:,i_th);
                norm_I = vectors_Fro(i_th);
                for j_th = i_th:num_Block
                    vector_J = vectors(:,j_th);
                    norm_J = vectors_Fro(j_th);
                    current_Theta = acos((vector_I'*vector_J)/(norm_I*norm_J));
                    sin_Theta = real(sin(current_Theta));
                    cos_Theta = real(cos(current_Theta));
        % J0(¦È) = ¦Ð - ¦È,
                    if length(option.vector_W)>=1 && option.vector_W(1,1)~=0
                        jtheta_Matrix0(i_th,j_th) = pi - current_Theta;
                        jtheta_Matrix0(j_th,i_th) = jtheta_Matrix0(i_th,j_th);
                    end
        % J1(¦È) = (¦Ð - ¦È) cos ¦È + sin¦È,
                    if length(option.vector_W)>=2 && option.vector_W(1,2)~=0
                        jtheta_Matrix1(i_th,j_th) = (pi - current_Theta)*cos_Theta - sin_Theta;
                        jtheta_Matrix1(j_th,i_th) = jtheta_Matrix1(i_th,j_th);
                    end
        % J2(¦È) = (¦Ð - ¦È)(1 + 2 cos^2 ¦È) + 3sin¦È cos ¦È,
                    if length(option.vector_W)>=3 && option.vector_W(1,3)~=0
                        jtheta_Matrix2(i_th,j_th) = (pi - current_Theta)*(1 + 2*cos_Theta^2) + 3*sin_Theta*cos_Theta;
                        jtheta_Matrix2(j_th,i_th) = jtheta_Matrix2(i_th,j_th);
                    end
        % J3(¦È) = (¦Ð - ¦È)(9 sin^2 ¦È cos ¦È + 15cos^3 ¦È) + 4sin^3 ¦È + 15sin¦È cos^2 ¦È.
                    if length(option.vector_W)>=4 && option.vector_W(1,4)~=0
                        jtheta_Matrix3(i_th,j_th) = (pi - current_Theta)*(9*sin_Theta^2*cos_Theta+15*cos_Theta^3) + 4*sin_Theta^3 + 15*sin_Theta*cos_Theta^2;
                        jtheta_Matrix3(j_th,i_th) = jtheta_Matrix3(i_th,j_th);
                    end
                end        
            end  
            jtheta_Matrix0 = real(jtheta_Matrix0);jtheta_Matrix1 = real(jtheta_Matrix1);
            jtheta_Matrix2 = real(jtheta_Matrix2);jtheta_Matrix3 = real(jtheta_Matrix3);
            k_Matrix0 = (1/pi)*inner_Matrix.^(0).*(jtheta_Matrix0); 
            k_Matrix1 = (1/pi)*inner_Matrix.^(1).*(jtheta_Matrix1); 
            k_Matrix2 = (1/pi)*inner_Matrix.^(2).*(jtheta_Matrix2); 
            k_Matrix3 = (1/pi)*inner_Matrix.^(3).*(jtheta_Matrix3); 

            iCovDs = option.vector_W(1,1)*k_Matrix0 + option.vector_W(1,2)*k_Matrix1 + option.vector_W(1,3)*k_Matrix2 + option.vector_W(1,4)*k_Matrix3; 
    end
% %% third release 
%     vectors = reshape(log_Matrix, size(log_Matrix, 1) * size(log_Matrix, 2), size(log_Matrix, 3));
%     vectors_Norm = sqrt(sum(vectors.*vectors,1)); 
%     R0 = 0; R1 = 1; R2 = 2; R3 = 3;
%     k_Matrix0 = zeros(num_Block,num_Block); k_Matrix1 = zeros(num_Block,num_Block);
%     k_Matrix2 = zeros(num_Block,num_Block); k_Matrix3 = zeros(num_Block,num_Block);
%     for i_th = 1:num_Block
%         vector_I = vectors(:,i_th);
%         norm_I = vectors_Norm(i_th);
%         for j_th = i_th:num_Block
%             vector_J = vectors(:,j_th);
%             norm_J = vectors_Norm(j_th);
%             theta_Euclidean = acos((vector_I'*vector_J)/(norm_I*norm_J));
%             sin_Theta = sin(theta_Euclidean);
%             cos_Theta = cos(theta_Euclidean);
% % J0(¦È) = ¦Ð - ¦È,
% % r = 0 ,
%             j_R0_L0 = pi - theta_Euclidean;
%             k_R0_L1_IJ = 1/pi*norm_I^R0*norm_J^R0*j_R0_L0;
%             k_R0_L1_II = 1/pi*norm_I^R0*norm_I^R0*j_R0_L0;
%             k_R0_L1_JJ = 1/pi*norm_J^R0*norm_J^R0*j_R0_L0;
%             [k_R0_L2_IJ,k_R0_L2_II,k_R0_L2_JJ] = com_arccos_infinite(k_R0_L1_IJ,k_R0_L1_II,k_R0_L1_JJ,R0);
%             [k_R0_L3_IJ,k_R0_L3_II,k_R0_L3_JJ] = com_arccos_infinite(k_R0_L2_IJ,k_R0_L2_II,k_R0_L2_JJ,R0);
%             [k_R0_L4_IJ,k_R0_L4_II,k_R0_L4_JJ] = com_arccos_infinite(k_R0_L3_IJ,k_R0_L3_II,k_R0_L3_JJ,R0);           
%             k_Matrix0(i_th,j_th) = k_R0_L4_IJ;
%             k_Matrix0(j_th,i_th) = k_Matrix0(i_th,j_th);
% % J1(¦È) = (¦Ð - ¦È) cos ¦È + sin¦È,
%             j_R1_L0 = (pi - theta_Euclidean)*cos_Theta - sin_Theta;
%             k_R1_L1_IJ = 1/pi*norm_I^R1*norm_J^R1*j_R1_L0;
%             k_R1_L1_II = 1/pi*norm_I^R1*norm_I^R1*j_R1_L0;
%             k_R1_L1_JJ = 1/pi*norm_J^R1*norm_J^R1*j_R1_L0;
%             [k_R1_L2_IJ,k_R1_L2_II,k_R1_L2_JJ] = com_arccos_infinite(k_R1_L1_IJ,k_R1_L1_II,k_R1_L1_JJ,R1);
%             [k_R1_L3_IJ,k_R1_L3_II,k_R1_L3_JJ] = com_arccos_infinite(k_R1_L2_IJ,k_R1_L2_II,k_R1_L2_JJ,R1);
%             [k_R1_L4_IJ,k_R1_L4_II,k_R1_L4_JJ] = com_arccos_infinite(k_R1_L3_IJ,k_R1_L3_II,k_R1_L3_JJ,R1);         
%             k_Matrix1(i_th,j_th) = k_R1_L4_IJ;
%             k_Matrix1(j_th,i_th) = k_Matrix1(i_th,j_th);
% % J2(¦È) = (¦Ð - ¦È)(1 + 2 cos^2 ¦È) + 3sin¦È cos ¦È,
%             j_R2_L0 = (pi - theta_Euclidean)*(1 + 2*cos_Theta^2) + 3*sin_Theta*cos_Theta;
%             k_R2_L1_IJ = 1/pi*norm_I^R2*norm_J^R2*j_R2_L0;
%             k_R2_L1_II = 1/pi*norm_I^R2*norm_I^R2*j_R2_L0;
%             k_R2_L1_JJ = 1/pi*norm_J^R2*norm_J^R2*j_R2_L0;
%             [k_R2_L2_IJ,k_R2_L2_II,k_R2_L2_JJ] = com_arccos_infinite(k_R2_L1_IJ,k_R2_L1_II,k_R2_L1_JJ,R2);
%             [k_R2_L3_IJ,k_R2_L3_II,k_R2_L3_JJ] = com_arccos_infinite(k_R2_L2_IJ,k_R2_L2_II,k_R2_L2_JJ,R2);
%             [k_R2_L4_IJ,k_R2_L4_II,k_R2_L4_JJ] = com_arccos_infinite(k_R2_L3_IJ,k_R2_L3_II,k_R2_L3_JJ,R2); 
%             k_Matrix2(i_th,j_th) = k_R2_L4_IJ;
%             k_Matrix2(j_th,i_th) = k_Matrix2(i_th,j_th);
% % J3(¦È) = (¦Ð - ¦È)(9 sin^2 ¦È cos ¦È + 15cos^3 ¦È) + 4sin^3 ¦È + 15sin¦È cos^2 ¦È.
%             j_R3_L0 = (pi - theta_Euclidean)*(9*sin_Theta^2*cos_Theta+15*cos_Theta^3) + 4*sin_Theta^3 + 15*sin_Theta*cos_Theta^2;
%             k_R3_L1_IJ = 1/pi*norm_I^R3*norm_J^R3*j_R3_L0;
%             k_R3_L1_II = 1/pi*norm_I^R3*norm_I^R3*j_R3_L0;
%             k_R3_L1_JJ = 1/pi*norm_J^R3*norm_J^R3*j_R3_L0;
%             [k_R3_L2_IJ,k_R3_L2_II,k_R3_L2_JJ] = com_arccos_infinite(k_R3_L1_IJ,k_R3_L1_II,k_R3_L1_JJ,R3);
%             [k_R3_L3_IJ,k_R3_L3_II,k_R3_L3_JJ] = com_arccos_infinite(k_R3_L2_IJ,k_R3_L2_II,k_R3_L2_JJ,R3);
%             [k_R3_L4_IJ,k_R3_L4_II,k_R3_L4_JJ] = com_arccos_infinite(k_R3_L3_IJ,k_R3_L3_II,k_R3_L3_JJ,R3); 
%             k_Matrix3(i_th,j_th) = k_R3_L4_IJ;
%             k_Matrix3(j_th,i_th) = k_Matrix3(i_th,j_th);
%         end        
%     end  
%     k_Matrix0 = real(k_Matrix0);k_Matrix1 = real(k_Matrix1);
%     k_Matrix2 = real(k_Matrix2);k_Matrix3 = real(k_Matrix3);
%     switch option.order
%         case 0
%             iCovDs_Arccos = (k_Matrix0);
%         case 1
%             iCovDs_Arccos = (k_Matrix0+k_Matrix1);
%         case 2
%             iCovDs_Arccos = (k_Matrix0+k_Matrix1+k_Matrix2);
%         case 3
%             iCovDs_Arccos = (k_Matrix0+k_Matrix1+k_Matrix2+k_Matrix3);        
%     end
% 
end