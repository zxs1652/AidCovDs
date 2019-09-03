function  Fea_each_im = IndividualModelIm(I,D,F,feat_prms,single_desc_dim,in_triu,shift)

y = F(1,:);
x = F(2,:);
y = y - min(y) + 1 ;
x = x - min(x) + 1;

pgrid = feat_prms.spm.^2;
Fea_each_im = zeros(sum(pgrid)*single_desc_dim,1);
% pool local sift descriptors over a spatial pyramid

weights = (1./feat_prms.spm); % divided by the number of grids at the corresponding level
weights = weights/sum(weights);

counter = 1;
for s = 1:length(feat_prms.spm)
       width = (size(I,2)-1)/feat_prms.spm(s);
       height = (size(I,1)-1)/feat_prms.spm(s);
       xgrid = ceil(x/height);
       ygrid = ceil(y/width);
       allgrid = (ygrid -1 )*feat_prms.spm(s) + xgrid;
       for t = 1:pgrid(s)
           range = counter:counter+single_desc_dim-1;
           ind = find(allgrid == t);
           if(numel(ind)>1)
               theD = D(:,ind);
               if strcmp(feat_prms.feature_type,'mean')==1
                    regionD = mean(theD);
               elseif strcmp(feat_prms.feature_type,'cov_mean')==1
                   % second-order pooling
                   regionD = real(logm(((1/size(theD,2)).*(theD *theD')) + shift)); 
                   regionD = regionD(in_triu);              
               elseif strcmp(feat_prms.feature_type,'cov')==1
                   regionD = real(logm((cov(theD')) + shift)); 
                   regionD = regionD(in_triu);  
               elseif strcmp(feat_prms.feature_type,'gaussian')==1
                   true_mat = true(size(D,1)+1,size(D,1)+1);
                   in_triu = triu(true_mat);
                   regionD = ones(size(theD,1));
                   regionD(1:size(theD,1),1:size(theD,1)) = (cov(theD') + shift); 
                   regionD(1:size(theD,1),1+size(theD,1)) =  mean(theD);
                   regionD(1+size(theD,1),1:size(theD,1)) =  mean(theD)';
                   const = max(det(P),esp(4))^(-1/(size(theD,1)+1));
                   regionD = const.*regionD;
                   regionD = real(logm(regionD) + shift); 
                   regionD = regionD(in_triu);  
               end
               Fea_each_im(range) = weights(s) * regionD(:);
           end
           counter = counter + single_desc_dim;
       end
end
end