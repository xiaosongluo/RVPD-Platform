for id = 1 : 1 :103
    index = int2str(id);
    
    datasetFolder = '..\Datasets\Fang';
    resultFolder = '..\Results\Temp';
    dataFolder='..\Results\Fang';
    
    imageFile = fullfile(datasetFolder, [index, '.jpg']);
    doim = imread(imageFile);
    
    %Proposed
    dataFile = fullfile(dataFolder,'Proposed', [index, '.mat']);
    load(dataFile);
    doim = mask(vpx, vpy, doim, 255, 0, 0);
    %Kong
    dataFile = fullfile(dataFolder,'Kong', [index, '.mat']);
    load(dataFile);
    doim = mask(vpx*2, vpy*2, doim, 0, 255, 0);
    %Moghamdam
    dataFile = fullfile(dataFolder,'Moghadam', [index, '.mat']);
    load(dataFile);
    doim = mask(vpx*2, vpy*2, doim, 0, 0, 255);
    
    imwrite(uint8(doim), fullfile(resultFolder, [index, '.jpg']), 'jpg');
end