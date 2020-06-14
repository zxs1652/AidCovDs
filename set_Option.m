%%  set option for different datasets
%   Written by Kai-Xuan Chen (e-mail: kaixuan_chen_jsh@163.com)
%   version 2.0 -- April/2019 
%   version 1.0 -- December/2017 
% 
%   input : 
%         param  : the name of dataset
%   output:
%         option : the struct of parameters  


function option = set_Option(name)
    root = pwd;
    addpath(fullfile(root, 'Tools'));
    addpath(fullfile(root, 'Tools\Gabor'));
    addpath(fullfile(root, 'Tools\vlfeat-0.9.18\toolbox'));
    addpath(fullfile(root, 'Tools\vlfeat-0.9.18\app\recognition'));
    addpath(genpath(fullfile(root, 'Tools\tool_Kit\KSVM')));
    addpath(genpath(fullfile(root, 'Tools\tool_Kit\CDL')));
    switch name
        case 'CG'
            option.mat_Path = '.\Data\mat_CG';                      % relative path of the resulting data
            option.mat_File = [option.mat_Path,'\','CG.mat'];       % '.mat' files of the resulting data
            option.res_Path = '.\Data\res_CG';                      % relative path of accuracies
            option.res_File = [option.res_Path,'\','CG.mat'];       % '.mat' files of accuracies
            option.root_Path = 'E:\WORKSPACE\DATASET\CGesture\';    % rootpath of dataset in your computer
            option.num_Class = 9;                                   % number of categories
            option.num_Sample = 100;                                % number of sample in each class
            option.num_Train = 20;                                  % number of gallery sample
            option.num_Test = option.num_Sample - option.num_Train; % number of probe sample
            option.label_Train = reshape(ones(option.num_Train,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Train]);
                                                                    % label of gallery sample
            option.label_Test = reshape(ones(option.num_Test,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Test]);
                                                                    % label of probe sample
            option.pre_Class = '';                                  % prefic string of each class
            option.pre_Set = 's';                                   % prefix string of each image set  
            option.type_Image = '.jpg';                             % suffix of images
            option.type_C = 'RGB';                                  % image type
            option.num_Basis = 100;                                  % basis number for nystrom
            option.times_Basis = 1.5;                               % train number for nystrom
            option.num_Ite = 120;                                   % number of iterations                           
            option.scale = 5;                                       % parameter for Gabor
            option.orientation = 8;                                 % parameter for Gabor
            option.width = 39;                                      % parameter for Gabor
            option.type_F = 'Sift';                                % feature type 
            option.step = 4;                                        % parameter for Sift
            option.scales = [2,1.4142,1,0.7071,0.5];                % parameter for Sift
            
        case 'ETH'
            option.mat_Path = '.\Data\mat_ETH';
            option.mat_File = [option.mat_Path,'\','ETH.mat'];
            option.res_Path = '.\Data\res_ETH';
            option.res_File = [option.res_Path,'\','ETH.mat'];
            option.root_Path = '.\ETH-80\';
            option.num_Class = 8;
            option.num_Sample = 10;
            option.num_Train = 5;
            option.num_Test = option.num_Sample - option.num_Train;
            option.label_Train = reshape(ones(option.num_Train,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Train]);
            option.label_Test = reshape(ones(option.num_Test,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Test]);
            option.pre_Class = '';
            option.pre_Set = '';
            option.type_Image = '.png';
            option.type_C = 'RGB';
            option.num_Basis = 200;
            option.times_Basis = 1.2;
            option.num_Ite = 100;
            option.scale = 5;
            option.orientation = 8;
            option.width = 31;
            option.type_F = 'Sift';
            option.step = 4;                                        % parameter for Sift
            option.scales = [2,1.4142,1,0.7071,0.5];                % parameter for Sift
            
        case 'Virus'
            option.mat_Path = '.\Data\mat_Virus';
            option.mat_File = [option.mat_Path,'\','Virus.mat'];
            option.res_Path = '.\Data\res_Virus';
            option.res_File = [option.res_Path,'\','Virus.mat'];
            option.root_Path = 'E:\WORKSPACE\DATASET\Virus\Virus\';
            option.num_Class = 15;
            option.num_Sample = 5;
            option.num_Train = 2;
            option.num_Test = option.num_Sample - option.num_Train;
            option.label_Train = reshape(ones(option.num_Train,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Train]);
            option.label_Test = reshape(ones(option.num_Test,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Test]);
            option.pre_Class = '';
            option.pre_Set = '';
            option.type_Image = '.png';
            option.type_C = 'Gray';
            option.num_Basis = 50;
            option.times_Basis = 1.7;
            option.num_Ite = 10;
            option.scale = 5;
            option.orientation = 8;
            option.width = 13;
            option.type_F = 'Gabor';
            option.step = 4;                                        % parameter for Sift
            option.scales = [2,1.4142,1,0.7071,0.5];                % parameter for Sift
            
        case 'MDSD'
            option.mat_Path = '.\Data\mat_MDSD';
            option.mat_File = [option.mat_Path,'\','MDSD.mat'];
            option.res_Path = '.\Data\res_MDSD';
            option.res_File = [option.res_Path,'\','MDSD.mat'];
            option.root_Path = 'E:\WORKSPACE\DATASET\MDSD\';
            option.num_Class = 13;
            option.num_Sample = 10;
            option.num_Train = 7;
            option.num_Test = option.num_Sample - option.num_Train;
            option.label_Train = reshape(ones(option.num_Train,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Train]);
            option.label_Test = reshape(ones(option.num_Test,option.num_Class)*diag([1:option.num_Class]),[1,option.num_Class*option.num_Test]);
            option.pre_Class = 'c';
            option.pre_Set = 's';
            option.type_Image = '.jpg';
            option.type_C = 'RGB';
            option.num_Basis = 200;
            option.times_Basis = 1.2;
            option.num_Ite = 100;
            option.scale = 5;
            option.orientation = 8;
            option.width = 39;
            option.type_F = 'Local';
            option.step = 4;                                        % parameter for Sift
            option.scales = [2,1.4142,1,0.7071,0.5];                % parameter for Sift
    end
    if ~exist(option.mat_Path,'dir')     
        mkdir (option.mat_Path)
    end
    if ~exist(option.res_Path,'dir')     
        mkdir (option.res_Path)
    end
    option.name_Dataset = name;
    option.dis_Matrix_Path = [option.mat_Path,'\','disMatrix','_',name,'.mat'];
end