%% Author: Kai-Xuan Chen 
% Date: 2018.08.27
% combine the distance matrix of parpool

function out_Matrix = combine_Dis_Matrix(option,workers)
    
    num_Total_Samples = option.num_Sample*option.num_Class;
    dis_Matrix_A = zeros(num_Total_Samples,num_Total_Samples);
    dis_Matrix_S = zeros(num_Total_Samples,num_Total_Samples);
    dis_Matrix_J = zeros(num_Total_Samples,num_Total_Samples);
    dis_Matrix_L = zeros(num_Total_Samples,num_Total_Samples);
    for worker_th = 1:workers
        current_Worker_Dis_Matrix = [option.mat_Path,'\','disMatrix','_','par',num2str(worker_th),'_',option.name_Dataset,'.mat'];
        load(current_Worker_Dis_Matrix);
        dis_Matrix_A = dis_Matrix_A + dis_Matrix.A; 
        dis_Matrix_S = dis_Matrix_S + dis_Matrix.S;
        dis_Matrix_J = dis_Matrix_J + dis_Matrix.J;
        dis_Matrix_L = dis_Matrix_L + dis_Matrix.L;
        clear('dis_Matrix','current_Worker');
    end
    dis_Matrix.A = dis_Matrix_A;
    dis_Matrix.S = dis_Matrix_S;
    dis_Matrix.J = dis_Matrix_J;
    dis_Matrix.L = dis_Matrix_L;
    dis_Matrix_Output = option.dis_Matrix_Path;
    save(dis_Matrix_Output,'dis_Matrix');
    out_Matrix = dis_Matrix;
end