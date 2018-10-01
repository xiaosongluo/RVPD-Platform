%--------------------------------------------------------------------------
% Parameters Setting
%--------------------------------------------------------------------------
% Workshop Parameters
datasetBaseFolder = '..\Datasets';
sourceBaseFolder = '..\Methods';
resultBaseFolder = '..\Results';

% Select the Data Set
itemSet = dir(datasetBaseFolder);
itemName = {itemSet.name};
selection = listdlg('PromptString','select a dataset:', 'SelectionMode','single','ListString',itemName);
dataset = itemSet(selection).name;

% Select the Method
itemSet = dir(sourceBaseFolder);
itemName = {itemSet.name};
selection = listdlg('PromptString','select a method:', 'SelectionMode','single','ListString',itemName);
method = itemSet(selection).name;

% Task Parameters
datasetFolder = fullfile(datasetBaseFolder, dataset);
sourceFolder = fullfile(sourceBaseFolder, method);
resultFolder = fullfile(resultBaseFolder, dataset, method);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
end

% Select the Image File
[filename, pathname] = uigetfile(datasetFolder);
ind = strfind(filename, '.');
index = filename(1 : ind-1);
id = str2double(index);

%--------------------------------------------------------------------------
% Main
%--------------------------------------------------------------------------
% Add Path to Worksapce
addpath(sourceFolder);

[total, scale, gxa, gya] = configure(dataset);

imageFile = fullfile(datasetFolder, [index, '.jpg']);
reportFile = fullfile(resultFolder, [index, '.mat']);
if exist(imageFile,'file')
    srcImage = imresize(imread(imageFile), scale);
    tic;
    [vpx, vpy] = vp(srcImage, resultFolder, index);
    cost = toc;
    gx=gxa(id);
    gy=gya(id);
    error = geterror(vpx, vpy, gx, gy, srcImage, resultFolder, index);
    save(reportFile, 'vpx', 'vpy', 'gx', 'gy', 'cost', 'error');
    disp(['Image ', index, ' in ', dataset, ' : done in ', num2str(cost), ' seconds with error of ', num2str(error)]);
end

% Remove Path to Worksapce
rmpath(sourceFolder);