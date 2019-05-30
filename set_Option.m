%%  set option for different datasets
%   input : 
%         param  : the name of dataset
%   output:
%         option : the struct of parameters  
%%   Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
%   version 2.0 -- April/2019 
%   version 1.0 -- December/2017 

function option = set_Option(name)

    addpath(genpath('Tools'));  
    addpath(genpath('Mat_File'));  
    switch name
        case 'ETH'
            option.mat_Path = '.\Data\mat_ETH';                     % relative path of the resulting data
            option.mat_File = [option.mat_Path,'\','ETH.mat'];      % '.mat' files of the resulting data
            option.res_Path = '.\Data\res_ETH';                     % relative path of accuracies
            option.res_File = [option.res_Path,'\','ETH.mat'];      % '.mat' files of accuracies
            option.root_Path = '.\ETH-80\';% rootpath of dataset in your computer
            option.num_Class = 8;                                   % number of categories
            option.num_Sample = 10;                                 % number of sample in each class
            option.num_Train = 5;                                   % number of gallery sample
            option.num_Test = option.num_Sample - option.num_Train; % number of probe sample
            option.label_Train = reshape(ones(option.num_Train,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Train]);
                                                                    % label of gallery sample
            option.label_Test = reshape(ones(option.num_Test,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Test]);
                                                                    % label of probe sample
            option.pre_Class = '';                                  % prefix string of each class
            option.pre_Set = '';                                    % prefix string of each image set  
            option.type_Image = '.png';                             % suffix of images
            option.type_C = 'RGB';                                  % image type
            option.num_Basis = 100;                                 % basis number for nystrom
            option.times_Basis = 1.2;                               % train number for nystrom
            option.num_Ite = 100;                                   % number of iterations  
            option.scale = 5;                                       % parameter for Gabor
            option.orientation = 8;                                 % parameter for Gabor
            option.width = 39;                                      % parameter for Gabor
            option.type_F = 'Sift';                                 % feature type 
            option.size = 4;                                        % parameter for Sift
            option.step = 4;                                        % parameter for Sift
            
    end
    if ~exist(option.mat_Path,'dir')     
        mkdir (option.mat_Path)
    end
    if ~exist(option.res_Path,'dir')     
        mkdir (option.res_Path)
    end
end