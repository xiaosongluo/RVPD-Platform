function [total, scale, gx, gy] = configure(dataset)
% Workshop Parameters
datasetBaseFolder = '..\Datasets';
datasetFolder = fullfile(datasetBaseFolder, dataset);

switch dataset 
    case 'Moghadam'
        % Get Total Test Number
        total = 501;
        
        % Get Scale
        scale =[240, 320];
        
        % Get Gx & Gy
        fid = fopen(fullfile(datasetFolder, 'ground_truth.txt'));
        data = fscanf(fid, '%d');
        fclose(fid);
        gx = data(2:2:total*2);
        gy = data(1:2:total*2);
    case 'Fang'
        % Get Total Test Number
        total = 103;
        
        % Get Scale
        scale =[240, 320];
        
        % Get Gx & Gy
        gx = zeros(total);
        gy = zeros(total);
    case 'Others'
        % Get Total Test Number
        total = 99;
        
        % Get Scale
        scale =[320, 640];
        
        % Get Gx & Gy
        gx = zeros(total);
        gy = zeros(total);
end 