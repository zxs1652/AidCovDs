% function [beta] = sc_approx_pooling(feaSet, B, pyramid, gamma, knn)
function [beta] = sc_approx_pooling(features, imageSize, codebook, pyramid, gamma, knn)
%================================================
% 
% Usage:
% Compute the linear spatial pyramid feature using sparse coding. 
%
% Inputss:
%   feaSet        -structure defining the feature set of an image   
%                       .feaArr     local feature array extracted from the
%                                   image, column-wise
%                       .x          x locations of each local feature, 2nd
%                                   dimension of the matrix
%                       .y          y locations of each local feature, 1st
%                               dimension of the matrix
%                       .width      width of the image
%                       .height     height of the image
%   codebook      -sparse dictionary, column-wise
%   pyramid       -defines structure of pyramid 
%   gamma         -sparsity regularization parameter
%   knn           -k nearest neighbors selected for sparse coding
% 
% Output:
%   beta          -multiscale max pooling feature
%
% Written by Jianchao Yang @ NEC Research Lab America (Cupertino)
% Mentor: Kai Yu
% July 2008
%
% Revised May. 2010
%===============================================

dSize = size(codebook, 2);
nSmp = size(features.descr, 2);
img_width = imageSize(2);
img_height = imageSize(1);
idxBin = zeros(nSmp, 1);

sc_codes = zeros(dSize, nSmp);

% compute the local feature for each local feature
D = features.descr'*codebook;
IDX = zeros(nSmp, knn);
for ii = 1:nSmp,
	d = D(ii, :);
	[dummy, idx] = sort(d, 'descend');
	IDX(ii, :) = idx(1:knn);
end

for ii = 1:nSmp,
    y = features.descr(:, ii);
    idx = IDX(ii, :);
    BB = codebook(:, idx);
    sc_codes(idx, ii) = feature_sign(BB, y, 2*gamma);
end

sc_codes = abs(sc_codes);

% spatial levels
pLevels = length(pyramid);
% spatial bins on each level
pBins = pyramid.^2;
% total spatial bins
tBins = sum(pBins);

beta = zeros(dSize, tBins);
bId = 0;

for iter1 = 1:pLevels,
    
    nBins = pBins(iter1);
    
    wUnit = img_width / pyramid(iter1);
    hUnit = img_height / pyramid(iter1);
    
    % find to which spatial bin each local descriptor belongs
    xBin = ceil(features.frame(1,:) / wUnit);
    yBin = ceil(features.frame(2,:) / hUnit);
    idxBin = (yBin - 1)*pyramid(iter1) + xBin;
    
    for iter2 = 1:nBins,     
        bId = bId + 1;
        sidxBin = find(idxBin == iter2);
        if isempty(sidxBin),
            continue;
        end      
        beta(:, bId) = max(sc_codes(:, sidxBin), [], 2);
    end
end

if bId ~= tBins,
    error('Index number error!');
end

beta = beta(:);
beta = beta./sqrt(sum(beta.^2));
