function [Cov_RandomD,Log_RandomD] = read_NystromKmeans(option)
% ,Cov_RieD,Log_RieD
%     Cov_Basis = zeros(option.scale*option.orientation,option.scale*option.orientation,option.num_Basis*option.times_Basis);
%     Log_Basis = zeros(option.scale*option.orientation,option.scale*option.orientation,option.num_Basis*option.times_Basis);
    num_Tr = option.num_Basis*option.times_Basis;
    num_Eclass = fix(num_Tr/option.num_Class);
    num_Random = num_Tr-(num_Eclass*option.num_Class);
    for basis_th = 1:num_Random
        index_Class = randperm(option.num_Class);
        current_Class = [option.root_Path,option.pre_Class,num2str(index_Class(1)),'\'];
        index_Set = randperm(option.num_Sample);
        current_Set = [current_Class,option.pre_Set,num2str(index_Set(1)),'\'];
        current_Images = dir(strcat(current_Set,'*',option.type_Image));
        index_Image =  randperm(length(current_Images));
        current_Image = imread([current_Set,current_Images(index_Image(1)).name]);
        if size(current_Image,3)==3
            current_Image = rgb2gray(current_Image);
        else
            current_Image = current_Image;
        end
        switch option.type_F
            case ('Gabor')
                gabor_Image = ZF_GaborFilter(double(current_Image),option.scale,option.orientation,option.width,option.width);
                gabor_Image = reshape(gabor_Image,[size(gabor_Image,1)*size(gabor_Image,2),option.scale*option.orientation]);
                gabor_Image = zscore(gabor_Image);
                Cov = cov(gabor_Image);
                Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                Log = logm(Cov);
            case ('Sift')
                sift_Image =  DSIFT(current_Image,option.size,option.step);
                sift_Image = zscore(double(sift_Image)');
                Cov = cov(sift_Image);
                Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                Log = logm(Cov);
        end
        Log_Basis(:,:,basis_th) = Log;
    end
    for class_th = 1:option.num_Class
        current_Class = [option.root_Path,option.pre_Class,num2str(class_th),'\'];
        fprintf('--read %d th class for k-means--------\n',class_th);
        image_InClass = 0;
        for set_th = 1:option.num_Sample
            fprintf('-------read %d th set for k-means--------\n',set_th);
            current_Set = [current_Class,option.pre_Set,num2str(set_th),'\'];
            current_Images = dir(strcat(current_Set,'*',option.type_Image));
            
%             switch option.type_F
%                 case ('Gabor')
%                     image_Matrix = zeros(option.scale*option.orientation,option.scale*option.orientation,length(current_Images));
%                 case ('Sift')
%                     image_Matrix = zeros(128,128,length(current_Images));
%             end
            for image_th = 1:length(current_Images);
                current_Image = imread([current_Set,current_Images(image_th).name]);
                image_InClass = image_InClass + 1;
                if size(current_Image,3)==3
                    current_Image = rgb2gray(current_Image);
                else
                    current_Image = current_Image;
                end
                switch option.type_F
                    case ('Gabor')
                        gabor_Image = ZF_GaborFilter(double(current_Image),option.scale,option.orientation,option.width,option.width);
                        gabor_Image = reshape(gabor_Image,[size(gabor_Image,1)*size(gabor_Image,2),option.scale*option.orientation]);
                        gabor_Image = zscore(gabor_Image);
                        Cov = cov(gabor_Image);
                        Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                        Log = logm(Cov);
                    case ('Sift')
                        sift_Image =  DSIFT(current_Image,option.size,option.step);
                        sift_Image = zscore(double(sift_Image)');
                        Cov = cov(sift_Image);
                        Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                        Log = logm(Cov);
                end
                Log = Log/(sqrt(sum(sum(Log.*Log))));
                image_Matrix(:,:,image_InClass) = Log;                
            end                     
        end
        mean_Logs = mean(image_Matrix,3);
        num_Current_Class = size(image_Matrix,3);
        for log_th = 1:num_Current_Class
           dis_Current_Class(1,log_th) = norm(image_Matrix(:,:,log_th)-mean_Logs,'fro');
        end
        [~,index] = sort(dis_Current_Class);
        index = index(1,1:num_Eclass);
        Log_Basis(:,:,(num_Random+1+(class_th-1)*num_Eclass:num_Random+class_th*num_Eclass)) = image_Matrix(:,:,index);
    end
    
    gallery_Matrix = reshape(Log_Basis,[size(Log_Basis,1)*size(Log_Basis,2),option.num_Basis*option.times_Basis])';    
    gallery_Matrix = bsxfun(@rdivide, gallery_Matrix, sqrt(sum(gallery_Matrix.^2, 2)));%和第24行的效果一样
    gallery_Kernel = gallery_Matrix*gallery_Matrix';
    [V,E,~] = svd(gallery_Kernel);
    p_Nystrom = diag(1./sqrt(diag(E(1:option.num_Basis,1:option.num_Basis))))*V(:,1:option.num_Basis)';
    num_Images = zeros(option.num_Class,option.num_Sample);
    for class_th = 1:option.num_Class
        current_Class = [option.root_Path,option.pre_Class,num2str(class_th),'\'];
        fprintf('--read %d th class --------\n',class_th);
        for set_th = 1:option.num_Sample
            fprintf('-------read %d th set --------\n',set_th);
            current_Set = [current_Class,option.pre_Set,num2str(set_th),'\'];
            current_Images = dir(strcat(current_Set,'*',option.type_Image));
            num_Images(class_th,set_th) = length(current_Images);
            switch option.type_F
                case ('Gabor')
                    probe_Matrix = zeros(option.scale*option.orientation,option.scale*option.orientation,length(current_Images));
                case ('Sift')
                    probe_Matrix = zeros(128,128,length(current_Images));
            end
            for image_th = 1:length(current_Images);
                current_Image = imread([current_Set,current_Images(image_th).name]);
                if size(current_Image,3)==3
                    current_Image = rgb2gray(current_Image);
                else
                    current_Image = current_Image;
                end
                switch option.type_F
                    case ('Gabor')
                        gabor_Image = ZF_GaborFilter(double(current_Image),option.scale,option.orientation,option.width,option.width);
                        gabor_Image = reshape(gabor_Image,[size(gabor_Image,1)*size(gabor_Image,2),option.scale*option.orientation]);
                        gabor_Image = zscore(gabor_Image);
                        Cov = cov(gabor_Image);
                        Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                        Log = logm(Cov);
                    case ('Sift')
                        sift_Image =  DSIFT(current_Image,option.size,option.step);
                        sift_Image = zscore(double(sift_Image)');
                        Cov = cov(sift_Image);
                        Cov = Cov+0.001*trace(Cov)*eye(size(Cov));
                        Log = logm(Cov);
                end
                Log = Log/(sqrt(sum(sum(Log.*Log))));
                probe_Matrix(:,:,image_th) = Log;                
            end 
            probe_Matrix = reshape(probe_Matrix,[size(probe_Matrix,1)*size(probe_Matrix,2),size(probe_Matrix,3)]);
            probe_Gallery_Kernel = gallery_Matrix*probe_Matrix;
            reduce_Probe_Gallery = p_Nystrom*probe_Gallery_Kernel;
            reduce_Probe_Gallery = zscore(reduce_Probe_Gallery')';
            Cov_Set = cov(reduce_Probe_Gallery');
            Cov_Set = Cov_Set+0.001*trace(Cov_Set)*eye(size(Cov_Set));
            Log_Set = logm(Cov_Set);
            Cov_RandomD{class_th,set_th} = Cov_Set;
            Log_RandomD{class_th,set_th} = Log_Set;            
        end        
    end
end