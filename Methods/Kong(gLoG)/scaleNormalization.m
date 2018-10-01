function [normalizedImg] = scaleNormalization(img, lowBound, uppBound)

%%% img is a grey-level one
imgH = size(img,1);
imgW = size(img,2);

maxV = max(max(img));
minV = min(min(img));

normalizedImg = lowBound + (img - minV)*(uppBound-lowBound)/(maxV-minV);
