function [vpY, vpX] = vp(originImg, resultFolder, index)

norient = 36;

colorImg = originImg;
    
if size(colorImg,3)>1
    img = rgb2gray(colorImg);
else
    img = colorImg;
end

outputPath_final = [resultFolder,'\'];

[vpY, vpX] = vpeDOIm_gabor_road_detection(img,colorImg, norient, outputPath_final, str2double(index));