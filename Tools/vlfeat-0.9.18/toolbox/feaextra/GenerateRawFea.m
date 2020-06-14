function  [RawD RawF] = GenerateRawFea(I,feat_prms) 

if strcmp(feat_prms.raw_feature,'sift')==1
    
    enrichments = {};
    [RawD,RawF] = compute_shape_invariant_feats(I, 'really_dense_sift', {enrichments}, 'gray', [], feat_prms.step, feat_prms.patch_size);
    
elseif strcmp(feat_prms.raw_feature,'esift')==1
    
    enrichments = {'rgb', 'hsv', 'lab', 'xy_fullimg', 'scale_fullimg'};
    [RawD,RawF] = compute_shape_invariant_feats(I,'really_dense_sift', {enrichments}, 'gray', [], feat_prms.step, feat_prms.patch_size);
    
elseif strcmp(feat_prms.raw_feature,'ehog')==1
    
    enrichments = {'rgb', 'hsv', 'lab', 'xy_fullimg', 'scale_fullimg'};
    [RawD,RawF] = compute_shape_invariant_feats(I,'really_dense_hog', {enrichments}, 'gray', [], feat_prms.step, feat_prms.patch_size);
    
elseif strcmp(feat_prms.raw_feature,'elbp')==1
    
    enrichments = {'rgb', 'hsv', 'lab', 'xy_fullimg', 'scale_fullimg'};
    [RawD,RawF] = compute_shape_invariant_feats(I,'really_dense_lbp', {enrichments}, 'gray', [], feat_prms.step, feat_prms.patch_size);
    
elseif strcmp(feat_prms.raw_feature,'L2ECM')==1
    
    [RawD, RawF] = computelocalfeat(I, feat_prms.patch_size(2)^2,feat_prms.patch_size,'CM');
    
elseif strcmp(feat_prms.raw_feature,'L2EMG')==1
    
    [RawD, RawF] = computelocalfeat(I, feat_prms.patch_size(2)^2,feat_prms.patch_size,'MG');
    
end