function distance = geterror(vX, vY, gX, gY, originImg, resultFolder, index)

doim = originImg;

%Ground Truth
doim = mask(gX, gY, doim, 0, 255, 0);

%Result Point
doim = mask(vX, vY, doim, 255, 0, 255);

imwrite(uint8(doim), fullfile(resultFolder, [index, 'Compare.jpg']), 'jpg');

% Get Answer
[imgH, imgW] = size(originImg);
diag = sqrt(imgH^2 + imgW^2);
dist = sqrt((vX - gX)^2 + (vY - gY)^2);
distance =  dist / diag;