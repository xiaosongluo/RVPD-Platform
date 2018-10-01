

inputPath = 'C:\Users\hkong\myWork2011\TMI_blob\TestImages\Road Images\'; %%% input image path
outputPath_final = 'C:\Users\hkong\myWork2011\TMI_blob\TestImages\vp2\'; %%% result path


for i=429:580
    i
    
    if i<10
        imgName = [inputPath, 'image','000', int2str(i), '.jpg'];
    elseif i<100
        imgName = [inputPath, 'image','00', int2str(i), '.jpg'];
    elseif i<1000
        imgName = [inputPath, 'image','0', int2str(i), '.jpg'];
    else
        imgName = [inputPath, 'image',int2str(i), '.jpg'];
    end
    
    
    if exist(imgName,'file')>0
        
        
        
    colorImg = imread(imgName);
    if size(colorImg,3)>1
       img = rgb2gray(colorImg);
    else
        img = colorImg;
    end

    
    img = imresize(img, [180, 240]); 
    colorImg = imresize(colorImg, [180, 240]); %%[180, 240]
    
    imgH = 180;
    imgW = 240;
    
    imgName = [outputPath_final, int2str(i), 'vp1VotingMap_ini.jpg'];
    tmpVotingMap = double(imread(imgName));


www=5;
sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end


max_votingMap = max(max(sum_votingMap));
[r,c] = find(sum_votingMap==max_votingMap);


if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
else
    r = r;
    c = c;
end



doim = colorImg;

coord3 = [r,c];
if (coord3(1)<=7)
        coord3(1) = 8;
    end
    if (coord3(1)>=imgH-8)
        coord3(1)=imgH-8;
    end
    
    if (coord3(2)<=7)
        coord3(2) = 8;
    end
    
    if (coord3(2)>=imgW-8)
        coord3(2) = imgW-8;
    end
    y = coord3(1);
    x = coord3(2);
    doim(y-3:y+3,x-7:x+7,1) = 255;
    doim(y-7:y+7,x-3:x+3,1) = 255;
    doim(y-3:y+3,x-7:x+7,2) = 0;
    doim(y-7:y+7,x-3:x+3,2) = 0;
    doim(y-3:y+3,x-7:x+7,3) = 255;
    doim(y-7:y+7,x-3:x+3,3) = 255;
    
imwrite(uint8(doim), [outputPath_final,int2str(i),'vp1Map_ini.jpg'], 'jpg');










imgName = [outputPath_final, int2str(i), 'vp1VotingMap_ini1.jpg'];
    tmpVotingMap = double(imread(imgName));


www=5;
sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end


max_votingMap = max(max(sum_votingMap));
[r,c] = find(sum_votingMap==max_votingMap);


if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
else
    r = r;
    c = c;
end



doim = colorImg;

coord3 = [r,c];
if (coord3(1)<=7)
        coord3(1) = 8;
    end
    if (coord3(1)>=imgH-8)
        coord3(1)=imgH-8;
    end
    
    if (coord3(2)<=7)
        coord3(2) = 8;
    end
    
    if (coord3(2)>=imgW-8)
        coord3(2) = imgW-8;
    end
    y = coord3(1);
    x = coord3(2);
    doim(y-3:y+3,x-7:x+7,1) = 255;
    doim(y-7:y+7,x-3:x+3,1) = 255;
    doim(y-3:y+3,x-7:x+7,2) = 0;
    doim(y-7:y+7,x-3:x+3,2) = 0;
    doim(y-3:y+3,x-7:x+7,3) = 0;
    doim(y-7:y+7,x-3:x+3,3) = 0;
    
imwrite(uint8(doim), [outputPath_final,int2str(i),'vp1Map_ini1.jpg'], 'jpg');






imgName = [outputPath_final, int2str(i), 'vp1VotingMap_ini_dual.jpg'];
    tmpVotingMap = double(imread(imgName));


www=5;
sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end


max_votingMap = max(max(sum_votingMap));
[r,c] = find(sum_votingMap==max_votingMap);


if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
else
    r = r;
    c = c;
end



doim = colorImg;

coord3 = [r,c];
if (coord3(1)<=7)
        coord3(1) = 8;
    end
    if (coord3(1)>=imgH-8)
        coord3(1)=imgH-8;
    end
    
    if (coord3(2)<=7)
        coord3(2) = 8;
    end
    
    if (coord3(2)>=imgW-8)
        coord3(2) = imgW-8;
    end
    y = coord3(1);
    x = coord3(2);
    doim(y-3:y+3,x-7:x+7,1) = 0;
    doim(y-7:y+7,x-3:x+3,1) = 0;
    doim(y-3:y+3,x-7:x+7,2) = 0;
    doim(y-7:y+7,x-3:x+3,2) = 0;
    doim(y-3:y+3,x-7:x+7,3) = 255;
    doim(y-7:y+7,x-3:x+3,3) = 255;
    
imwrite(uint8(doim), [outputPath_final,int2str(i),'vp1Map_ini_dual.jpg'], 'jpg');



    end
end
