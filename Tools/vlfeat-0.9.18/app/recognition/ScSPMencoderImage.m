function descrs = ScSPMencoderImage(encoder, im, varargin)
% ENCODEIMAGE   Apply an encoder to an image
%   DESCRS = ENCODEIMAGE(ENCODER, IM) applies the ENCODER
%   to image IM, returning a corresponding code vector PSI.
%
%   IM can be an image, the path to an image, or a cell array of
%   the same, to operate on multiple images.
%
%   ENCODEIMAGE(ENCODER, IM, CACHE) utilizes the specified CACHE
%   directory to store encodings for the given images. The cache
%   is used only if the images are specified as file names.
%
%   See also: TRAINENCODER().

% Author: Andrea Vedaldi

% Copyright (C) 2013 Andrea Vedaldi
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

opts.cacheDir = [] ;
opts.cacheChunkSize = 512 ;
opts = vl_argparse(opts,varargin) ;

if ~iscell(im), im = {im} ; end

% break the computation into cached chunks
startTime = tic ;
descrs = cell(1, numel(im)) ;
numChunks = ceil(numel(im) / opts.cacheChunkSize) ;

for c = 1:numChunks
  n  = min(opts.cacheChunkSize, numel(im) - (c-1)*opts.cacheChunkSize) ;
  chunkPath = fullfile(opts.cacheDir, sprintf('chunk-%03d.mat',c)) ;
  if ~isempty(opts.cacheDir) && exist(chunkPath)
    fprintf('%s: loading descriptors from %s\n', mfilename, chunkPath) ;
    load(chunkPath, 'data') ;
  else
    range = (c-1)*opts.cacheChunkSize + (1:n) ;
    fprintf('%s: processing a chunk of %d images (%3d of %3d, %5.1fs to go)\n', ...
      mfilename, numel(range), ...
      c, numChunks, toc(startTime) / (c - 1) * (numChunks - c + 1)) ;
    data = processChunk(encoder, im(range)) ;
    if ~isempty(opts.cacheDir)
      save(chunkPath, 'data') ;
    end
  end
  descrs{c} = data ;
  clear data ;
end
descrs = cat(1,descrs{:}) ;

% --------------------------------------------------------------------
function psi = processChunk(encoder, im)
% --------------------------------------------------------------------
psi = cell(1,numel(im)) ;
% means = bsxfun(@times, encoder.means, 1./max(1e-12,sqrt(sum(encoder.means.^2)))) ;
% covariances = bsxfun(@times, encoder.covariances, 1./max(1e-12, sqrt(sum(encoder.covariances.^2)))) ;
% codebook = [means;covariances];
codebook = encoder.words;
% codebook = bsxfun(@times, codebook, 1./max(1e-12, sqrt(sum(codebook.^2)))) ;
if numel(im) > 1 & matlabpool('size') > 1
  parfor i = 1:numel(im)
    psi{i} = encodeOne(encoder, im{i},codebook) ;
  end
else
  % avoiding parfor makes debugging easier
  for i = 1:numel(im)
    psi{i} = encodeOne(encoder, im{i},codebook) ;
  end
end
psi = cat(1, psi{:}) ;

% --------------------------------------------------------------------
function psi = encodeOne(encoder, im, codebook)
% --------------------------------------------------------------------

im = encoder.readImageFn(im) ;

features = encoder.extractorFn(im) ;
features.descr = encoder.projection * bsxfun(@minus, features.descr, encoder.projectionCenter) ;
if encoder.renormalize
    features.descr = bsxfun(@times, features.descr, 1./max(1e-12, sqrt(sum(features.descr.^2)))) ;
end

% descrs = features.descr;
% descrs = bsxfun(@times, descrs, 1./max(sqrt(sum(descrs.^2)), 1e-12)) ;
imageSize = size(im) ;
% features.frame(1,:) = features.frame(1,:) / imageSize(2);
% features.frame(2,:) = features.frame(2,:) / imageSize(1);
% psi = {} ;
% X = [];
% arr=[];x=[];y=[];

% for i = 1:size(encoder.subdivisions,2)
%   minx = encoder.subdivisions(1,i) * imageSize(2) ;
%   miny = encoder.subdivisions(2,i) * imageSize(1) ;
%   maxx = encoder.subdivisions(3,i) * imageSize(2) ;
%   maxy = encoder.subdivisions(4,i) * imageSize(1) ;
% 
%   ok = ...
%     minx <= features.frame(1,:) & features.frame(1,:) < maxx  & ...
%     miny <= features.frame(2,:) & features.frame(2,:) < maxy ;
% 
%   descrs = encoder.projection * bsxfun(@minus, ...
%                                        features.descr(:,ok), ...
%                                        encoder.projectionCenter) ;
%   if sum(ok) > 3
% %       descrs = features.descr(:,ok);
%       if encoder.renormalize
%         descrs = bsxfun(@times, descrs, 1./max(1e-12, sqrt(sum(descrs.^2)))) ;
%       end
% 
%     %   w = size(im,2) ;
%     %   h = size(im,1) ;
%     %   frames = features.frame(1:2,ok) ;
%     %   frames = bsxfun(@times, bsxfun(@minus, frames, [w;h]/2), 1./[w;h]) ;
% 
%     %   descrs = extendDescriptorsWithGeometry(encoder.geometricExtension, frames, descrs) ;
% 
% %       v = var(descrs')' ;
% %       [means, covariances] = vl_gmm(descrs, 1, 'verbose', 'Initialization', 'kmeans', 'CovarianceBound', double(max(v)*0.0001), 'NumRepetitions', 1) ;
%       means = mean(descrs');
%       covariances = std(descrs');
%       means = means / max(1e-12,sqrt(sum(means.^2)));
%       covariances = covariances / max(1e-12,sqrt(sum(covariances.^2)));
%       X = [means,covariances];   
%       X = X / max(sqrt(sum(X.^2)), 1e-12) ;
%       arr = [arr;X];
%       x = [x,(encoder.subdivisions(1,i) + encoder.subdivisions(3,i)) / 2];
%       y = [y,(encoder.subdivisions(2,i) + encoder.subdivisions(4,i)) / 2];
%       
% %   z = means;
% 
% %   switch encoder.type
% %     case 'bovw'
% %       [words,distances] = vl_kdtreequery(encoder.kdtree, encoder.words, ...
% %                                          descrs, ...
% %                                          'MaxComparisons', 100) ;
% %       z = vl_binsum(zeros(encoder.numWords,1), 1, double(words)) ;
% %       z = sqrt(z) ;
% % 
% %     case 'fv'
% %       z = vl_fisher(descrs, ...
% %                     encoder.means, ...
% %                     encoder.covariances, ...
% %                     encoder.priors, ...
% %                     'Improved') ;
% %     case 'vlad'
% %       [words,distances] = vl_kdtreequery(encoder.kdtree, encoder.words, ...
% %                                          descrs, ...
% %                                          'MaxComparisons', 15) ;
% %       assign = zeros(encoder.numWords, numel(words), 'single') ;
% %       assign(sub2ind(size(assign), double(words), 1:numel(words))) = 1 ;
% %       z = vl_vlad(descrs, ...
% %                   encoder.words, ...
% %                   assign, ...
% %                   'SquareRoot', ...
% %                   'NormalizeComponents') ;
% %   end
%   
%   end
% end

%% LLC
% knn = 5;
% beta = 1e-4;
% SPM_layouts = {'1x1','2x2','3x1'};
% fea = LLC_coding_appr(codebook',arr,knn,beta);
% % psi = LLC_pooling(fea,features.frame(1,:),features.frame(2,:),SPM_layouts);
% psi = LLC_pooling(fea,x,y,SPM_layouts);

%% ScSPM 
pyramid = [1,2,4];
gamma = 0.15;
knn = 200;
psi = sc_approx_pooling(features, imageSize, codebook, pyramid, gamma, knn);
psi = psi';


% function psi = LLC_pooling(fea,fea_x,fea_y,SPM_layouts)
% psi=[];
% spm_subdivisions=[];
% for i = 1:numel(SPM_layouts)
%   t = sscanf(SPM_layouts{i},'%dx%d') ;
%   m = t(1) ;
%   n = t(2) ;
%   [x,y] = meshgrid(linspace(0,1,n+1), linspace(0,1,m+1)) ;
%   x1 = x(1:end-1,1:end-1) ;
%   y1 = y(1:end-1,1:end-1) ;
%   x2 = x(2:end,2:end) ;
%   y2 = y(2:end,2:end) ;
%   spm_subdivisions = cat(2, spm_subdivisions, [x1(:)' ; y1(:)' ; x2(:)' ; y2(:)'] ) ;
% end
% 
% for i = 1:size(spm_subdivisions,2)
%   minx = spm_subdivisions(1,i);
%   miny = spm_subdivisions(2,i);
%   maxx = spm_subdivisions(3,i);
%   maxy = spm_subdivisions(4,i);
% 
%   ok = minx <= fea_x & fea_x < maxx  & miny <= fea_y & fea_y < maxy ;
%   temp = max(fea(ok,:),[],1);
% %   temp =  temp / max(1e-12,sqrt(sum(temp.^2)));   %不归一化更好一点
%   psi = [psi,temp];
% end

