function [vpY, vpX] = vp(originImg, resultFolder, index)

% Get Origin Image
imwrite(originImg, fullfile(resultFolder, [index, 'OriginImg.jpg']), 'jpg');

% Get Gray Image
grayImg = rgb2gray(originImg);
imwrite(grayImg, fullfile(resultFolder, [index, 'GrayImg.jpg']), 'jpg');

% Construct Kernels
[convolution, halfKerSize] = vp_convolution(grayImg);

% Get Dominant Orientation & Confidence
orientation = vp_dominant(convolution, halfKerSize);
imwrite(uint8(orientation), fullfile(resultFolder, [index, 'Orientation.jpg']), 'jpg');

% Show Orientation
orientbar = vp_orientbar(grayImg, orientation);
imwrite(uint8(orientbar), fullfile(resultFolder, [index, 'OrientionBar.jpg']), 'jpg');

% Voting for Vanishing Point
[votingMap, vpX, vpY] = vp_vote(orientation, halfKerSize);
imwrite(uint8(votingMap), fullfile(resultFolder, [index, 'VotingMap.jpg']), 'jpg');