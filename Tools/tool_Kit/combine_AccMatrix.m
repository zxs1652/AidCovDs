function out_MST = combine_AccMatrix(option,num_Folds)
    temp_AccMatrix = [];temp_TiMatrix = [];
    for fold_th = 1:num_Folds
        temp_Option = set_Option('YTC',fold_th);
        res_DisPath = temp_Option.res_Output;
        if ~exist(res_DisPath,'file')
            return;
        end
        load(res_DisPath);
        temp_AccMatrix = [temp_AccMatrix,accracy_Matrix];
        temp_TiMatrix = [temp_TiMatrix,time_Matrix];
        clear('res_DisPath','accracy_Matrix','time_Matrix');
        
    end
    accracy_Matrix = temp_AccMatrix;
    time_Matrix = temp_TiMatrix;
    [MST.mean_Accuracy,mean_Accuracy] = deal(mean(accracy_Matrix,2)); % mean accracy
    [MST.std_Accuracy,std_Accuracy] = deal(std(accracy_Matrix,0,2));  % standard deviation
    [MST.mean_Time,mean_Time] = deal(mean(time_Matrix,2));            % time matrix
    if ~exist(option.res_Path,'dir')
        mkdir(option.res_Path);
    end
    res_AllFold_Path = [option.res_Path,'\',option.name_Dataset,'_numGallery',num2str(option.num_Gallery),'_resized',num2str(option.resized_Row),'M',num2str(option.resized_Col),'_blockSize',...
        num2str(option.block_Row),'M',num2str(option.block_Col),'_stepSize',num2str(option.step_Row),'M',num2str(option.step_Col),'_beta',num2str(option.beta_Gauss),'.mat'];
    out_MST = MST;
    save(res_AllFold_Path,'MST','mean_Accuracy','std_Accuracy','accracy_Matrix','time_Matrix'); 
end