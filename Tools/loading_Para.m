function [option,Log_RandomD,output_DisMatrix] = loading_Para(option)
    if ~exist(option.res_Path,'dir')
        mkdir (option.res_Path);
    end
    option.num_Gallery = option.num_Train;
    option.num_Probe = option.num_Test;
    option.label_Gallery = option.label_Train;
    option.label_Probe = option.label_Test;  
    load (option.mat_File);
    if ~exist(option.dis_Matrix_Path,'file')
        workers = 4;
        pool = parpool(workers);    
        parfor worker_th = 1:workers
            temp_Dis_Matrix = compute_Dis_Par(option,worker_th,workers); % generate distance matrix, if it does not exist.   
        end
        delete(pool); 
        dis_Matrix = combine_Dis_Matrix(option,workers);
    else
        load (option.dis_Matrix_Path);
    end
    output_DisMatrix = dis_Matrix;
end