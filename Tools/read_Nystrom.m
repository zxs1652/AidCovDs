%%   Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
%   version 2.0 -- April/2019 
%   version 1.0 -- December/2017 
% 
%   input : 
%         option : the struct of parameters  
%   output : 
%         Cov_RandomD: the cell of SPD matrices
%         Log_RandomD: the cell of the logarithm of SPD matrices      



function [Cov_RandomD,Log_RandomD,Inv_RandomD] = read_Nystrom(option)
%%  Training for Nystrom
    for basis_th = 1:option.num_Basis*option.times_Basis
        index_Class = randperm(option.num_Class);                                           % random index for class
        current_Class = [option.root_Path,option.pre_Class,num2str(index_Class(1)),'\'];    % select one class
        index_Set = randperm(option.num_Sample);                                            % random index for set
        current_Set = [current_Class,option.pre_Set,num2str(index_Set(1)),'\'];             % select one set
        current_Images = dir(strcat(current_Set,'*',option.type_Image));
        index_Image =  randperm(length(current_Images));                                    % random index for images
        current_Image = imread([current_Set,current_Images(index_Image(1)).name]);          % select one image     
        switch option.type_F
            case ('Gabor')
                if size(current_Image,3)==3
                    current_Image = rgb2gray(current_Image);                % RGB to gray
                end
                gabor_Image = ZF_GaborFilter(double(current_Image),option.scale,option.orientation,option.width,option.width);  % extract Gabor feature                                                            
                gabor_Image = reshape(gabor_Image,[size(gabor_Image,1)*size(gabor_Image,2),option.scale*option.orientation]);
                gabor_Image = double(gabor_Image);                      
                Cov = fea2SPD(gabor_Image);                                 % Generate traditional CovDs
                Log = logm(Cov);                                            % Generate the logarithm matrix of CovDs
            case ('Sift')
                if size(current_Image,3)==3
                    current_Image = rgb2gray(current_Image);                % RGB to gray
                end
                sift = getDenseSIFT(current_Image,'step', option.step, 'scales', option.scales);
                sift_Image =  sift.descr;                                   % extract Sift feature
                sift_Image = double(sift_Image);
                Cov = fea2SPD(sift_Image');                                 % Generate traditional CovDs
                Log = logm(Cov);                                            % Generate the logarithm matrix of CovDs
            case ('Local')
                local_Image = compute_LocalFea(option,current_Image);       % extract Local feature
                local_Image = reshape(local_Image,[size(local_Image,1)*size(local_Image,2),size(local_Image,3)]);
                Cov = fea2SPD(local_Image);                                 % Generate traditional CovDs
                Log = logm(Cov);                                            % Generate the logarithm matrix of CovDs
        end
        Cov_Basis(:,:,basis_th) = Cov;                                      % Train CovDs for Nystorm
        Log_Basis(:,:,basis_th) = Log;                                      % Train logarithm matrices for Nystorm
    end
    gallery_Matrix = reshape(Log_Basis,[size(Log_Basis,1)*size(Log_Basis,2),option.num_Basis*option.times_Basis])';   
    gallery_Matrix = bsxfun(@rdivide, gallery_Matrix, sqrt(sum(gallery_Matrix.^2, 2))); % Constructing matrix of gallery samples
    gallery_Kernel = gallery_Matrix*gallery_Matrix';                        % kernel matrix of sampled images for Nystorm
    [V,E] = compute_svd(gallery_Kernel);
    if option.num_Basis > rank(gallery_Kernel)
        option.num_Basis = rank(gallery_Kernel);
    end
    p_Nystrom = diag(1./sqrt(diag(E(1:option.num_Basis,1:option.num_Basis))))*V(:,1:option.num_Basis)'; % projection matrix for Nystorm
    
%%  Generate approximate infinite-dimensional covariance descriptors (AidCovDs)
    for class_th = 1:option.num_Class
        current_Class = [option.root_Path,option.pre_Class,num2str(class_th),'\']; % path of current class
        fprintf('--read %d th class --------\n',class_th);
        for set_th = 1:option.num_Sample
            fprintf('-------read %d th set --------\n',set_th);
            current_Set = [current_Class,option.pre_Set,num2str(set_th),'\'];       % path of current set
            current_Images = dir(strcat(current_Set,'*',option.type_Image));        % images of current set
            % initialize matrix of probe samples
            switch option.type_F
                case ('Gabor')
                    probe_Matrix = zeros(option.scale*option.orientation,option.scale*option.orientation,length(current_Images));
                case ('Sift')
                    probe_Matrix = zeros(128,128,length(current_Images));
                case ('Local')
                    switch option.type_C
                        case ('RGB')
                            probe_Matrix = zeros(11,11,length(current_Images));
                        case ('Gray')
                            probe_Matrix = zeros(5,5,length(current_Images));
                    end
            end
            for image_th = 1:length(current_Images);
                current_Image = imread([current_Set,current_Images(image_th).name]);% matrix of current image
                switch option.type_F
                    case ('Gabor')
                        if size(current_Image,3)==3
                            current_Image = rgb2gray(current_Image);        % RGB to gray
                        end
                        gabor_Image = ZF_GaborFilter(double(current_Image),option.scale,option.orientation,option.width,option.width); % extract Gabor feature  
                        gabor_Image = reshape(gabor_Image,[size(gabor_Image,1)*size(gabor_Image,2),option.scale*option.orientation]);
                        gabor_Image = double(gabor_Image);                      
                        Cov = fea2SPD(gabor_Image);                         % Generate traditional CovDs
                        Log = logm(Cov);                                    % Generate the logarithm matrix of CovDs
                    case ('Sift')
                        if size(current_Image,3)==3
                            current_Image = rgb2gray(current_Image);        % RGB to gray
                        end
                        sift = getDenseSIFT(current_Image,'step', option.step, 'scales', option.scales);
                        sift_Image =  sift.descr;                                   % extract Sift feature
                        sift_Image = double(sift_Image);
                        Cov = fea2SPD(sift_Image');                         % Generate traditional CovDs
                        Log = logm(Cov);                                    % Generate the logarithm matrix of CovDs
                    case ('Local')
                        local_Image = compute_LocalFea(option,current_Image);           % extract local feature
                        local_Image = reshape(local_Image,[size(local_Image,1)*size(local_Image,2),size(local_Image,3)]);
                        Cov = fea2SPD(local_Image);                         % Generate traditional CovDs
                        Log = logm(Cov);                                    % Generate the logarithm matrix of CovDs
                end
                Log = Log/(sqrt(sum(sum(Log.*Log))));
                probe_Matrix(:,:,image_th) = Log;                           % Constructing matrix of probe samples            
            end 
            probe_Matrix = reshape(probe_Matrix,[size(probe_Matrix,1)*size(probe_Matrix,2),size(probe_Matrix,3)]);
            probe_Gallery_Kernel = gallery_Matrix*probe_Matrix;             % kernel matrix between gallery and probe samples 
            reduce_Probe_Gallery = p_Nystrom*probe_Gallery_Kernel;          % approximate feature matrix
            reduce_Probe_Gallery = zscore(reduce_Probe_Gallery')';
            Cov_Set = cov(reduce_Probe_Gallery');                           % AidCovDs of current set
            Cov_Set = Cov_Set+0.001*trace(Cov_Set)*eye(size(Cov_Set));      % add perturbation
            Log_Set = logm(Cov_Set);
            Inv_Set = Cov_Set^(-1);
            Cov_RandomD{class_th,set_th} = Cov_Set;
            Log_RandomD{class_th,set_th} = Log_Set; 
            Inv_RandomD{class_th,set_th} = Inv_Set; 
        end        
    end
end