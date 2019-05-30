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
option = set_Option('ETH');         % set Hy-Parameters
load (option.mat_File);             % load data for classifiers
load In_Matrix_ETH.mat;             % random index matrix
row_total = size(In_Matrix,1);      % In_Matrix
for ite_th = 1:option.num_Ite
    fprintf('------------ite_th :  %d------------------\n',ite_th);
    ind_Begin = (ite_th-1)*option.num_Class+1;    % The beginning of an arbitrary index 
    if ind_Begin >row_total
        ind_Begin = mod(ind_Begin,row_total);     % Re-assign if the 'ind_Begin' is greater than the row number of index matrix
    end
    ind_End = ind_Begin + option.num_Class -1;    % compute the end of index via 'ind_Begin'
    if ind_End > row_total                        % If the end of index is greater than the number of rows, re-assign the 'ind_Begin' and 'ind_End'
        ind_Begin = mod(ind_End,row_total);
        ind_End = ind_Begin + option.num_Class -1;
    end
    rand_Matrix = In_Matrix(ind_Begin:ind_End,:);  % The index matrix of the current iteration 

    tic;        % CDL (COV-LDA)
    accracy_Matrix(1,ite_th) = Classify_CDL_LDA(Log_RandomD,rand_Matrix,option);
    time_Matrix(1,ite_th) = toc;

    tic;        % CDL (COV-PLS)
    accracy_Matrix(2,ite_th) = Classify_CDL_PLS(Log_RandomD,rand_Matrix,option);
    time_Matrix(2,ite_th) = toc;

end
[MST.mean_Accuracy,mean_Accuracy] = deal(mean(accracy_Matrix,2)); % mean accracy
[MST.std_Accuracy,std_Accuracy] = deal(std(accracy_Matrix,0,2));  % standard deviation
[MST.mean_Time,mean_Time] = deal(mean(time_Matrix,2));            % time matrix
save(option.res_File,'MST','mean_Accuracy','std_Accuracy','accracy_Matrix','time_Matrix');