%%  Riemannian kernel based nystrom method for approximate infinite-dimensional covariance descriptors (AidCovDs) with application to image set classification
%   Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
%   version 2.0 -- April/2019 
%   version 1.0 -- December/2017 
%
%   Please cite the following paper (more theoretical and technical details) if your are using this code:
% 
%   Kai-Xuan Chen, Xiao-Jun Wu, Rui Wang, Josef Kittler. Riemannian kernel based nystrom method for approximate infinite- 
%   dimensional covariance descriptors with application to image set classification. In 2018 24th International Conference
%   on Pattern Recognition (ICPR), pages 651¨C656. IEEE, 2018.

clear;
clc;
workers = 4;
option = set_Option('ETH');             % set Hy-Parameters
if ~exist(option.mat_File,'file')
    [Cov_RandomD,Log_RandomD,Inv_RandomD] = read_Nystrom(option);     % Generate approximate infinite-dimensional covariance descriptors (AidCovDs) 
    %   output : 
    %         Cov_RandomD: the cell of SPD matrices
    %         Log_RandomD: the cell of the logarithm of SPD matrices
    %         Inv_RandomD: the cell of the inverse matrices of SPD matrices
    %   input:
    %         option : the struct of parameters 

    save(option.mat_File,'Cov_RandomD','Log_RandomD','Inv_RandomD');
end

if ~exist(option.dis_Matrix_Path,'file')
    pool = parpool(workers);    
%     parfor worker_th = 1:workers
    for worker_th = 1:workers
        temp_Dis_Matrix = compute_Dis_Par(option,worker_th,workers); % Generate distance matrix, if it does not exist.   
    end
    delete(pool);   
end
dis_Matrix = combine_Dis_Matrix(option,workers); % combine distance matrix