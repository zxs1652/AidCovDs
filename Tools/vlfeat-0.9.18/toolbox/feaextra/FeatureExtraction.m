function  features = FeatureExtraction(image_dir,data_dir,feat_prms,exprs) 


sub_image_dir = fullfile(image_dir, exprs.sub_dir{1});
fnames = dir(fullfile(sub_image_dir, '*.png')); 
imageFileList = cell( length(fnames),1);
for f = 1: length(fnames)
    imageFileList{f} = fnames(f).name;
end

imageFName = fullfile(sub_image_dir, imageFileList{1});
I = imread(imageFName);
[RawFeas RawPos] = GenerateRawFea(I,feat_prms);
single_desc_dim = (size(RawFeas,1)*(size(RawFeas,1)+1)/2);
feat_dim = sum((feat_prms.spm).^2) * single_desc_dim;
features = zeros(feat_dim, numel(imageFileList)*exprs.all_dir_per_category, 'single');
shift = feat_prms.constrant.*eye(size(RawFeas,1));
true_mat = true(size(RawFeas,1),size(RawFeas,1));
in_triu = triu(true_mat);
fprintf('Extracting %s features--Processing %s', feat_prms.feature_type,image_dir);
for ss = 1:length(exprs.sub_dir)
     sub_image_dir = fullfile(image_dir, exprs.sub_dir{ss});
     fnames = dir(fullfile(sub_image_dir, '*.png')); 
     imageFileList = cell( length(fnames),1);
     for f = 1: length(fnames)
         imageFileList{f} = fnames(f).name;
     end
     sub_data_dir = fullfile(data_dir,exprs.sub_dir{ss});
    for f =1:size(imageFileList,1)
        
        if ~mod(f, 5)
            fprintf('.');
        end

        curr_filename = imageFileList{f};
        suffix_pos = findstr(curr_filename, '.png');
        outFName = fullfile(sub_data_dir, [curr_filename(1:suffix_pos-1) '.mat']);
        if feat_prms.canSkip && size(dir(outFName), 1)
            load(outFName);
            features(:,f) = Fea_each_im;
            continue;
        end

        imageFName = fullfile(sub_image_dir, curr_filename);
        I = imread(imageFName);

        [hgt wid ~] = size(I);
        if max(hgt,wid) > feat_prms.maxImageSize    
            I = imresize(I, feat_prms.maxImageSize/max(hgt,wid), 'bicubic'); 
        end

    %     % resize image
    %     [im_h, im_w, ~] = size(I);
    %     
    %     if max(im_h, im_w) > feat_prms.maxImageSize,
    %         I = imresize(I, feat_prms.maxImageSize/max(im_h, im_w), 'bicubic');
    %         [im_h, im_w] = size(I);
    %     end;
    %     if min(im_h, im_w) < feat_prms.minImageSize,
    %         I = imresize(I, feat_prms.minImageSize/min(im_h, im_w), 'bicubic');
    %     end;

        [RawFeas RawPos] = GenerateRawFea(I,feat_prms);

        Fea_each_im = IndividualModelIm(I,RawFeas,RawPos,feat_prms,single_desc_dim,in_triu,shift);

        %% Normalization 
        Fea_each_im = squash_features(Fea_each_im,feat_prms.norm);

        Fea_each_im = single(Fea_each_im);
        features(:,(ss-1)*size(imageFileList,1)+f) = Fea_each_im;
        sp_make_dir(outFName);
        save(outFName, 'Fea_each_im');   
    end    
end
fprintf('Completed.\n');

%features = [];