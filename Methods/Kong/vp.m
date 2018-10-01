function [vpY, vpX] = vp(originImg, resultFolder, index)

norient = 36;

% Get Origin Image
imwrite(originImg, fullfile(resultFolder, [index, 'OriginImg.jpg']), 'jpg');

% Get Gray Image
grayImg = rgb2gray(originImg);
imwrite(grayImg, fullfile(resultFolder, [index, 'GrayImg.jpg']), 'jpg');

% Construct Kernels
[convolution, halfKerSize] = vp_convolution(grayImg, norient);

% Get Dominant Orientation & Confidence
[orientation, confidence] = vp_dominant(convolution, halfKerSize);
imwrite(uint8(orientation), fullfile(resultFolder, [index, 'Orientation.jpg']), 'jpg');
imwrite(uint8(confidence), fullfile(resultFolder, [index, 'Confidence.jpg']), 'jpg');

% Get Confidence Overlap
overlap = (confidence > 255*0.3);

% Show Orientation
orientbar = vp_orientbar(grayImg, orientation);
imwrite(uint8(orientbar), fullfile(resultFolder, [index, 'OrientionBar.jpg']), 'jpg');

% Get Voter Area
voter = vp_voter(grayImg, overlap);

% Show Voter Area
voterarea = vp_voterarea(grayImg, voter);
imwrite(uint8(voterarea), fullfile(resultFolder, [index, 'VoterArea.jpg']), 'jpg');

% Voting for Vanishing Point
[votingMap, vpX, vpY] = vp_vote(orientation, voter, halfKerSize);
imwrite(uint8(votingMap), fullfile(resultFolder, [index, 'VotingMap.jpg']), 'jpg');