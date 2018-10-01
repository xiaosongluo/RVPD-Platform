function [vpY, vpX] = vp(originImg, resultFolder, index)

norient = 36;
r = 12;
threshold = 0.05;

% Get Origin Image
imwrite(originImg, fullfile(resultFolder, [index, 'OriginImg.jpg']), 'jpg');

% Get Gray Image
grayImg = rgb2gray(originImg);
mgrayImg = medfilt2(grayImg,[5 5]);
imwrite(grayImg, fullfile(resultFolder, [index, 'GrayImg.jpg']), 'jpg');

% Get Edge Image
edgeImg = edge(grayImg,'canny');
imwrite(edgeImg, fullfile(resultFolder, [index, 'EdgeImg.jpg']), 'jpg');

% Get Orientation
orient = mwld_o(mgrayImg, norient, r);
imwrite(uint8(orient), fullfile(resultFolder, [index, 'Orient.jpg']), 'jpg');

% Get Differential Excitation
diffex = mwld_de(mgrayImg, r);
imwrite(uint8(diffex*255), fullfile(resultFolder, [index, 'DiffEx.jpg']), 'jpg');

% Get Confidence Overlap
overlap = (diffex >= threshold);

% Show Orientation
orientbar = vp_orientbar(grayImg, orient);
imwrite(uint8(orientbar), fullfile(resultFolder, [index, 'OrientBar.jpg']), 'jpg');

% Show Voter Area
voter = vp_voter(edgeImg, orient, overlap);
voterarea = vp_voterarea(grayImg, voter);
imwrite(uint8(voterarea), fullfile(resultFolder, [index, 'VoterArea.jpg']), 'jpg');

% Voting for Vanishing Point
[votingMap, vpX, vpY] = vp_vote(orient, voter);
imwrite(uint8(votingMap), jet(256), fullfile(resultFolder, [index, 'VotingMap.jpg']), 'jpg');