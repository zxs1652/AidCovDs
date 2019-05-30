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
option = set_Option('ETH');             % set Hy-Parameters
[Cov_RandomD,Log_RandomD] = read_Nystrom(option);     % Generate approximate infinite-dimensional covariance descriptors (AidCovDs) 
%   output : 
%         Cov_RandomD: the cell of SPD matrices
%         Log_RandomD: the cell of the logarithm of SPD matrices
%   input:
%         option : the struct of parameters 

save(option.mat_File,'Cov_RandomD','Log_RandomD');
