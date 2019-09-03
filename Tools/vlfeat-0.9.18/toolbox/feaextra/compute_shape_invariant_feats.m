function [D,F, feat_ranges, variable_grids] = compute_shape_invariant_feats(I, main_feat, enrichments, color_type, pca_basis, STEP, base_scales, F, variable_grids)
    DefaultVal('*variable_grids', 'false');
    DefaultVal('*F', '[]');    
    D = [];
    feat_ranges = {};
    
    fmax_i = size(I,1);
    fmin_i = 1;
    fmax_j = size(I,2);
    fmin_j = 1;
    width = (fmax_j - fmin_j);
    height = (fmax_i - fmin_i);
    
    if(strcmp(main_feat, 'really_dense_sift'))
        fuh = cell(1,numel(base_scales));
        duh = fuh;
        theI = I;
        total_pad = 0;
        for k=1:numel(base_scales)
            if(k~=1)
                pad = ceil(1.5*(base_scales(k)-base_scales(k-1)));
                total_pad = total_pad + pad;
                theI = padarray(theI, repmat(pad,1,2), 0, 'both');
            end                
            [fuh{k}, duh{k}] = vl_phow(im2single(theI), 'FloatDescriptors', true, 'Sizes', base_scales(k), 'STEP', STEP, 'Magnif', inf, 'Color', color_type, 'ContrastThreshold', 0,'Verbose',false);                
                
            fuh{k}(1:2,:) = fuh{k}(1:2,:) - total_pad;
        end
        var.F = cell2mat(fuh);
        var.D = cell2mat(duh);        
    
        %%modified by zenghui 2014/8/19
        var.D = var.D/255.0;
%         var.D = sqrt(var.D) ;
%     	var.D = bsxfun(@times, var.D, 1./max(1e-5,sqrt(sum(var.D.^2,1)))) ;
  
    elseif(strcmp(main_feat, 'really_dense_lbp'))
        if size(I,3)>1
            Isingle = single(rgb2gray(I));
        else
            Isingle = single(I);
        end
   
        var.D = [];
        var.F = [];
        variable_grids = true;
        for j=1:numel(base_scales)
            Ilbp = vl_lbp(Isingle, base_scales(j));
            [a,b] = find(Ilbp(:,:,1)~=inf); 
            newD = reshape(Ilbp, numel(a), size(Ilbp,3))';
            var.D = [var.D newD];
            a = a*base_scales(j);
            b = b*base_scales(j);
            newF = [b'; a'; base_scales(j)*ones(1,numel(a))];
            var.F = [var.F newF];
        end
      elseif(strcmp(main_feat, 'really_dense_hog'))
        if size(I,3)>1
            Isingle = single(rgb2gray(I));
        else
            Isingle = single(I);
        end
        
        var.D = [];
        var.F = [];
        variable_grids = true;
        for j=1:numel(base_scales)
            Ihog = vl_hog(Isingle, base_scales(j));
            [a,b] = find(Ihog(:,:,1)~=inf); 
            newD = reshape(Ihog, numel(a), size(Ihog,3))';
            var.D = [var.D newD];
            a = a*base_scales(j);
            a = min(a,size(Isingle,1)); 
            b = b*base_scales(j);
            b = min(b,size(Isingle,2)); 
            newF = [b'; a'; base_scales(j)*ones(1,numel(a))];
            var.F = [var.F newF];
        end
    else
        error('no such raw local descriptor');
    end
    
    D = var.D;
    F = var.F;
        
    all_enrich = [enrichments{:}];
    un_enrich = unique_no_sort(all_enrich);
    ranges = cell(numel(un_enrich)+1,1);
    ranges{1} = 1:size(D,1);
    
    for i=1:numel(un_enrich)
        % features that can be computed independently of mask
        if(strcmp(un_enrich{i}, 'texton_filterbank'))
            I_gray = im2single(rgb2gray(I));

            FB = FbMake(2,1);

            I_filtered = FbApply2d(I_gray, FB, 'same')*2; % multiply to get better range

            pos = F(1:2,un_f_range);
            if(~variable_grids)
                var.D = repmat(single(mean_nb_I(I_filtered, pos(2,:), pos(1,:), 0)), 1 ,numel(base_scales));
            end
        elseif(strcmp(un_enrich{i}, 'siftres'))
            var.D = [];
            for i=1:numel(base_scales)                    
                duh = reshape(D(1:128, F(end,:)==base_scales(i)), 128, max(F(1,:)/STEP), max(F(2,:)/STEP));
                duh = permute(duh, [3 2 1]);
                duh = vl_imsmooth(duh, 2);
                [Fx, Fy] = gradient(duh);
                new_duh = abs(Fy)+abs(Fx);
                new_duh = permute(new_duh, [3 2 1]);
                var.D = [var.D reshape(new_duh,128,size(duh,1)*size(duh,2))];
            end
        elseif(strcmp(un_enrich{i}, 'rgb_hist'))                
            for i=1:numel(base_scales)         
                error('not ready');
            end                
        elseif(strcmp(un_enrich{i}, 'xy_fullimg'))
            x = (F(1,:)-fmin_j) / max(1,width);
            y = (F(2,:)-fmin_i) / max(1,height);
            x2 = (F(1,:)-fmin_j) / max(1,height);
            y2 = (F(2,:)-fmin_i) / max(1,width);
            var.D = [x; y; x2; y2];
        elseif(strcmp(un_enrich{i}, 'scale_fullimg'))
            scale_1 = (F(end,:)*24)/height;
            scale_2 = (F(end,:)*24)/width;
            var.D = [scale_1; scale_2];
        else                
            % the following assume there some other features already there that populated F
            if(~variable_grids)
                un_f_range = 1:(size(F,2)/numel(base_scales));
            else
                un_f_range = 1:size(F,2);
            end
            pos = F(1:2,un_f_range);

            if(strcmp(un_enrich{i}, 'rgb'))                    
                thisI = I;
                if(size(thisI,3) == 1)
                    thisI = repmat(thisI, [1 1 3]);
                end                                                    

                m = single(mean_nb_I(thisI, pos(2,:), pos(1,:), 0))/255.0;
            elseif(strcmp(un_enrich{i}, 'hsv'))                                       
                thisI = I;
                if(size(thisI,3) == 1)
                    thisI = repmat(thisI, [1 1 3]);
                end

                thisI = rgb2hsv(thisI);                             
                m = single(mean_nb_I(thisI, pos(2,:), pos(1,:), 0));                    
            elseif(strcmp(un_enrich{i}, 'lab'))                                     
                thisI = I;
                if(size(thisI,3) == 1)
                    thisI = repmat(thisI, [1 1 3]);
                end

                thisI = rgb2lab(thisI);
                thisI(:,:,1) = thisI(:,:,1)/100;
                thisI(:,:,2) = (thisI(:,:,2)+110)/220;
                thisI(:,:,3) = (thisI(:,:,3)+110)/220;

                m = mean_nb_I(thisI, pos(2,:), pos(1,:), 0);
            elseif(strcmp(un_enrich{i}, 'ycbcr'))
                thisI = I;
                if(size(thisI,3) == 1)
                    thisI = repmat(thisI, [1 1 3]);
                end

                thisI = single(rgb2ycbcr(thisI))/255.0;
                m = single(mean_nb_I(thisI, pos(2,:), pos(1,:), 0));
            elseif(strcmp(un_enrich{i}, 'xyz'))
                thisI = I;
                if(size(thisI,3) == 1)
                    thisI = repmat(thisI, [1 1 3]);
                end

                C = makecform('srgb2xyz');
                thisI = applycform(thisI,C);
                thisI = single(thisI)/65536.0;

                m = single(mean_nb_I(thisI, pos(2,:), pos(1,:), 0));
            elseif(strcmp(un_enrich{i}, 'entropy'))
                pos = F(1:2,:);

                thisI = entropyfilt(rgb2gray(I));

                m = mean_nb_I(thisI, pos(2,:), pos(1,:), 0);
                m = single(m/6.0);  
            elseif(strcmp(un_enrich{i}, 'harris'))
                Igray = im2single(rgb2gray(I));
                si = [1.0 2.0 3.0];
                Iharr = zeros(size(Igray,1), size(Igray,2), numel(si));
                for k=1:3
                    Iharr(:,:,k) = vl_harris(Igray, si(k));
                end

                pos = F(1:2,:);
                m = mean_nb_I(Iharr, pos(2,:), pos(1,:), 0);
            elseif(strcmp(un_enrich{i}, 'grad'))
                pos = F(1:2,:);                    

                [grad_X, grad_Y] = gradient(rgb2gray(im2double(I)));                    
                YX = grad_X./grad_Y;
                A = ((atan(YX)+(pi/2))*180)/pi;
                A = A/180;
                B = ((atan2(grad_Y,grad_X)+pi)*180)/pi;
                B = B/360;
                grad_X = grad_X * 2; % get a range similar to the other features
                grad_Y = grad_Y* 2;

                I_grad = grad_X;
                I_grad = cat(3, I_grad, grad_Y);

                I_grad = reshape(I_grad, size(thisI,1)*size(thisI,2), size(I_grad,3))';
                ids = sub2ind(size(thisI(:,:,1)), pos(2,:), pos(1,:));
                m = I_grad(:, ids);
            elseif(strcmp(un_enrich{i}, 'colornames'))
                load('w2c.mat');
                N_CN = 11;
                I_double = double(I);
                colornames = zeros(size(I,1), size(I,2), N_CN);
                for cn=1:N_CN
                    colornames(:,:,cn) = im2c(I_double, w2c, cn);
                end
                m = mean_nb_I(colornames, pos(2,:), pos(1,:), 0);                
            else
                % not available
                m = [];
            end

            if(~variable_grids)
                var.D = repmat(single(m), 1, numel(base_scales));
            else
                var.D = single(m);
            end

            var.F = [];
        end        
               
        ranges{i+1} = max(cell2mat(ranges'))+1:max(cell2mat(ranges'))+size(var.D,1);
        D = [D; var.D];
        F = [F; var.F];        
    end       
    
    if(~isempty(pca_basis))
        error('fix this wrt the remaining lines of code');
        D = pca_basis*D;
    end    
    
    % there is one top-level cell per pooling domain (figure/ground/etc)
    for i=1:numel(enrichments)
        [~, ids] = intersect([main_feat un_enrich], [main_feat enrichments{i}]);
        feat_ranges{i} = sort(cell2mat(ranges(ids)'), 'ascend');
    end
end