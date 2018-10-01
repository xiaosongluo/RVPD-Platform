function voter = vp_voter(edgeImg, overlap)

[imgH, imgW] =size(edgeImg);

voter = ones(imgH,imgW).*overlap;