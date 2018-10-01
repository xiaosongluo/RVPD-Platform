function [vpY, vpX] = vp(originImg, resultFolder, index)

norient = 36;

colorImg = originImg;
    
if size(colorImg,3)>1
    img = rgb2gray(colorImg);
else
    img = colorImg;
end

[vpY, vpX] = vpeDOIm_gLoG(img,colorImg, norient, resultFolder, str2double(index));