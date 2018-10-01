function voterarea = vp_voterarea(grayImg, voter)

[imgH, imgW] = size(grayImg);

doim = zeros(imgH, imgW, 3);
doim(:,:,1) = grayImg;
doim(:,:,2) = grayImg;
doim(:,:,3) = grayImg;
doim(:,:,1) = doim(:,:,1).*(1 - voter) + voter * 255;
doim(:,:,2) = doim(:,:,2).*(1 - voter);
doim(:,:,3) = doim(:,:,3).*(1 - voter);

voterarea = uint8(doim);