function [vpY, vpX, vpY1, vpX1, vpX_dual, vpY_dual] = vpeDOIm_gLoG(grayImg,colorImg, norient, outputPath_final, numOfValidFiles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grayImgCopy = grayImg;
[imCopyH,imCopyW] = size(grayImgCopy);

largestDistence = sqrt(imCopyH*imCopyH + imCopyW*imCopyW);

colorImgCopy = colorImg;


outputPath = outputPath_final;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grayImg = grayImgCopy;
edgeImg = edge(grayImg,'canny');
imwrite(edgeImg,[outputPath, int2str(numOfValidFiles),'vpEdgeImg0.jpg'], 'jpg');






grayImgCopy = grayImg;
grayImg = medfilt2(grayImg,[5 5]);
[imCopyH,imCopyW] = size(grayImgCopy);
[imgH, imgW] = size(grayImgCopy);






%%%%% remove some vertical short (e.g., 30-pixel long) edges 
outlierBinary = ones(imCopyH,imCopyW);
candidateVotingMap = edgeImg; %%%zeros(imCopyH,imCopyW);
tmpVerBarLen = 30; %% odd number
halftmpVerBarLen = floor(tmpVerBarLen/2);
for tt=1+5:imCopyW-5
    for rr= halftmpVerBarLen+1:imCopyH-halftmpVerBarLen
        tmpSum = sum(edgeImg(rr-halftmpVerBarLen:rr+halftmpVerBarLen,tt));
        if tmpSum>=20
            candidateVotingMap(rr-halftmpVerBarLen:rr+halftmpVerBarLen,tt-5:tt+5)=0; %%labeled as 1
            outlierBinary(rr-halftmpVerBarLen:rr+halftmpVerBarLen,tt-5:tt+5)=0;
        end
    end
end
edgeImg = candidateVotingMap;




edgeImg(:,1:8)=0;
   edgeImg(:,imgW-7:imgW)=0;
   edgeImg(1:round(imCopyH*0.1),:)=0;

%%%%%%%%%% save the outliers image
outliers = zeros(imCopyH,imCopyW,3);
outliers(:,:,1) = grayImgCopy;
outliers(:,:,2) = grayImgCopy;
outliers(:,:,3) = grayImgCopy;
outliers(:,:,1) = outliers(:,:,1).*(1-candidateVotingMap) + candidateVotingMap*255;
% imwrite(uint8(outliers), [outputPath,int2str(numOfValidFiles),'vpOutliersImg.jpg'],'jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% construct kernels
angleRange = 180;
angleInterval = angleRange/norient;
tmpW = min([imCopyH,imCopyW]); %%min([imCopyH,imCopyW]);
lamda = 2^(round(log2(tmpW)-5));
kerSize = round(10*lamda/pi);

if mod(kerSize,2)==0
    kerSize = kerSize + 1;
    halfKerSize = floor(kerSize/2);
else
    halfKerSize = floor(kerSize/2);
end



cosTheta = zeros(angleRange, 1);
sinTheta = zeros(angleRange, 1);
for theta = 1:angleRange
    cosTheta(theta) = cos((theta-1)*pi/180);
    sinTheta(theta) = sin((theta-1)*pi/180);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% image convolution with gLoG filter banks %%%%%%%%%%%%%
grayImg1 = 255-double(grayImg);
imwrite(uint8(grayImg1),[outputPath, int2str(numOfValidFiles),'vpNegGrayImg.jpg'], 'jpg');

thetaStep = 3.14159/norient;
smallestSigma = 2;
largestSigma = 8;
complexResponse = zeros(imCopyH,imCopyW,norient);
complexResponse1 = zeros(imCopyH,imCopyW,norient);
ori_indx = 0;
for theta=0: thetaStep : pi-thetaStep
    ori_indx = ori_indx+1;
    [responseMap] = convolve_gLoG(double(grayImg), smallestSigma, largestSigma, theta);
    complexResponse(:,:,ori_indx) = responseMap;
    [responseMap] = convolve_gLoG(double(grayImg1), smallestSigma, largestSigma, theta);
    complexResponse1(:,:,ori_indx) = responseMap;
end




%%%%%%%%%%%%%%% the dominant orientation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
confidenceMap = zeros(imCopyH,imCopyW);
confidenceMap1 = zeros(imCopyH,imCopyW);

orientationMap = zeros(imCopyH,imCopyW);
orientationMap1 = zeros(imCopyH,imCopyW);
for i=1+halfKerSize:imCopyH-halfKerSize
    for j=1+halfKerSize:imCopyW-halfKerSize
%         maxV = max(complexResponse(i,j,:));
%         minV = min(complexResponse(i,j,:));
%         complexResponseVector = (complexResponse(i,j,:)-minV)*100/(maxV-minV);
%         maxLoc = find(complexResponseVector==100);
%         if length(maxLoc)>0
%             maxLoc = round(mean(maxLoc));
%             if (maxLoc>2)&&(maxLoc<=34)
%                 complexResponseVector(maxLoc-2:maxLoc+2)=0;
%             elseif (maxLoc<=2)
%                 complexResponseVector(maxLoc:maxLoc+2)=0;
%             elseif (maxLoc>=35)
%                 complexResponseVector(maxLoc-2:maxLoc)=0;
%             end
% 
%             [a,b] = sort(complexResponseVector,'descend');
%             tmpContrast = (100 - mean(a(5:15)))^2;
% 
%             
%             if maxV>1
%                 confidenceMap(i,j) = tmpContrast;
%             end
%         end
        
        maxma = max(complexResponse(i,j,:));
        indx = find(complexResponse(i,j,:)==maxma);
        orientationMap(i,j) = mean(indx); %%indx(1);
        
        maxma = max(complexResponse1(i,j,:));
        indx = find(complexResponse1(i,j,:)==maxma);
        orientationMap1(i,j) = mean(indx); %%indx(1);
    end
end
orientationMap = (orientationMap-1)*angleInterval;
orientationMap1 = (orientationMap1-1)*angleInterval;



maxV = max(max(confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize)));
minV = min(min(confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize)));
aaTmp = confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize);
aaTmp = (aaTmp-minV)*255/(maxV-minV);
confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize) = aaTmp;
imwrite(uint8(confidenceMap),[outputPath, int2str(numOfValidFiles),'vpConfidence.jpg'], 'jpg');


confidenceMapBinary = (confidenceMap>80);

doim = zeros(imCopyH, imCopyW, 3);
doim(:,:,1) = grayImgCopy;
doim(:,:,2) = grayImgCopy;
doim(:,:,3) = grayImgCopy;
doim(:,:,1) = doim(:,:,1).*(1-confidenceMapBinary) + confidenceMapBinary*255;
doim(:,:,2) = doim(:,:,2).*(1-confidenceMapBinary);
doim(:,:,3) = doim(:,:,3).*(1-confidenceMapBinary);
imwrite(uint8(doim),[outputPath, int2str(numOfValidFiles),'vpConfidenceOverlap.jpg'], 'jpg');


%%%%%%%%%%%%%%% show the orientation bar %%%%%%%%%%%%%%%%%%%%%%%
doim = grayImgCopy;
doim = zeros(imCopyH, imCopyW, 3);
doim(:,:,1) = grayImgCopy;
doim(:,:,2) = grayImgCopy;
doim(:,:,3) = grayImgCopy;

doim1 = doim;

orientationMapDisplay = 180-orientationMap;
orientationMapDisplay1 = 180-orientationMap1;
% for i=10:8:imCopyH-10
%     for j=10:8:imCopyW-10
%         ori = orientationMapDisplay(i,j); %%(orientationMap(i,j)-1)*angleInterval;
%         if (ori<=95)&&(ori>=85)
%             yy = i;
%             xx = j;
%             doim(yy:yy+7,xx-1:xx+1,1) = 255;
%             doim(yy:yy+7,xx-1:xx+1,2) = 0;
%             doim(yy:yy+7,xx-1:xx+1,3) = 0;
%         elseif ori==0
%             doim(i-1:i+1,j:j+7,1) = 255;
%             doim(i-1:i+1,j:j+7,2) = 0;
%             doim(i-1:i+1,j:j+7,3) = 0;
%         else
%             kk = tan(ori*pi/180);
%             for xx=j:j+7
%                 yy = round(kk*(xx-j) + i);
%                 if (yy>=i-7)&&(yy<=i+7) %%(yy>=1)&&(yy<=imCopyH)&&(xx>=2)&&(xx<=imCopyW-1) %%
%                 doim(yy,xx-1:xx+1,1) = 255;
%                 doim(yy,xx-1:xx+1,2) = 0;
%                 doim(yy,xx-1:xx+1,3) = 0;
%                 end
%             end
%         end
%         
%         ori1 = orientationMapDisplay1(i,j); %%(orientationMap(i,j)-1)*angleInterval;
%         if (ori1<=95)&&(ori1>=85)
%             yy = i;
%             xx = j;
%             doim1(yy:yy+7,xx-1:xx+1,1) = 255;
%             doim1(yy:yy+7,xx-1:xx+1,2) = 0;
%             doim1(yy:yy+7,xx-1:xx+1,3) = 0;
%         elseif ori1==0
%             doim1(i-1:i+1,j:j+7,1) = 255;
%             doim1(i-1:i+1,j:j+7,2) = 0;
%             doim1(i-1:i+1,j:j+7,3) = 0;
%         else
%             kk = tan(ori1*pi/180);
%             for xx=j:j+7
%                 yy = round(kk*(xx-j) + i);
%                 if (yy>=i-7)&&(yy<=i+7) %%(yy>=1)&&(yy<=imCopyH)&&(xx>=2)&&(xx<=imCopyW-1) %%
%                 doim1(yy,xx-1:xx+1,1) = 255;
%                 doim1(yy,xx-1:xx+1,2) = 0;
%                 doim1(yy,xx-1:xx+1,3) = 0;
%                 end
%             end
%         end
%             
%     end
% end

% imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vpOrientationBarImg.jpg'], 'jpg');
% imwrite(uint8(doim1), [outputPath,int2str(numOfValidFiles),'vpOrientationBarImg1.jpg'], 'jpg');

%%%%%%%%%%%%%%% show the orientation bar on resized image %%%%%%%%%%%%%%%%%%%%%%%
doim = zeros(imCopyH, imCopyW, 3);
doim(:,:,1) = grayImgCopy;
doim(:,:,2) = grayImgCopy;
doim(:,:,3) = grayImgCopy;
doim_resized = imresize(doim,2,'bilinear');
doim1_resized = doim_resized;
orientationMapDisplay_resized = imresize(orientationMapDisplay,2,'bilinear'); 
orientationMapDisplay1_resized = imresize(orientationMapDisplay1,2,'bilinear');

barLen = 4;
barWid = 1;
for i=10:8:imCopyH*2-10
    for j=10:8:imCopyW*2-10
        ori = orientationMapDisplay_resized(i,j); 
        if (ori==90)
            yy = i;
            xx = j;
            doim_resized(yy:yy+barLen,xx-barWid:xx,1) = 255;
            doim_resized(yy:yy+barLen,xx-barWid:xx,2) = 0;
            doim_resized(yy:yy+barLen,xx-barWid:xx,3) = 0;
        elseif (ori==180)||(ori==0)
            doim_resized(i-barWid:i,j:j+barLen,1) = 255;
            doim_resized(i-barWid:i,j:j+barLen,2) = 0;
            doim_resized(i-barWid:i,j:j+barLen,3) = 0;
        else
            if (ori<=45)||(ori>=135)
                kk = tan(ori*pi/180);
                for xx=j:j+barLen
                    yy = round(kk*(xx-j) + i);
                    if (yy>=i-barLen)&&(yy<=i+barLen) 
                    doim_resized(yy,xx-barWid:xx,1) = 255;
                    doim_resized(yy,xx-barWid:xx,2) = 0;
                    doim_resized(yy,xx-barWid:xx,3) = 0;
                    end
                end
            elseif (ori>45)&&(ori<135)
                kk = tan(ori*pi/180);
                for yy=i:i+barLen
                    xx = round((yy-i)/kk + j);
                    if (xx>=j-barLen)&&(xx<=j+barLen) 
                    doim_resized(yy-barWid:yy,xx,1) = 255;
                    doim_resized(yy-barWid:yy,xx,2) = 0;
                    doim_resized(yy-barWid:yy,xx,3) = 0;
                    end
                end
            end
        end
        
        ori1 = orientationMapDisplay1_resized(i,j); 
        if (ori1==90)
            yy = i;
            xx = j;
            doim1_resized(yy:yy+barLen,xx-barWid:xx,1) = 255;
            doim1_resized(yy:yy+barLen,xx-barWid:xx,2) = 0;
            doim1_resized(yy:yy+barLen,xx-barWid:xx,3) = 0;
        elseif (ori1==180)||(ori1==0)
            doim1_resized(i-barWid:i,j:j+barLen,1) = 255;
            doim1_resized(i-barWid:i,j:j+barLen,2) = 0;
            doim1_resized(i-barWid:i,j:j+barLen,3) = 0;
        else
            if (ori1<=45)||(ori1>=135)
                kk = tan(ori1*pi/180);
                for xx=j:j+barLen
                    yy = round(kk*(xx-j) + i);
                    if (yy>=i-barLen)&&(yy<=i+barLen) 
                    doim1_resized(yy,xx-barWid:xx,1) = 255;
                    doim1_resized(yy,xx-barWid:xx,2) = 0;
                    doim1_resized(yy,xx-barWid:xx,3) = 0;
                    end
                end
            elseif (ori1>45)&&(ori1<135)
                kk = tan(ori1*pi/180);
                for yy=i:i+barLen
                    xx = round((yy-i)/kk + j);
                    if (xx>=j-barLen)&&(xx<=j+barLen) 
                    doim1_resized(yy-barWid:yy,xx,1) = 255;
                    doim1_resized(yy-barWid:yy,xx,2) = 0;
                    doim1_resized(yy-barWid:yy,xx,3) = 0;
                    end
                end
            end
        end
            
    end
end

imwrite(uint8(doim_resized), [outputPath,int2str(numOfValidFiles),'vpOrientationBarImg.jpg'], 'jpg');
imwrite(uint8(doim1_resized), [outputPath,int2str(numOfValidFiles),'vpOrientationBarImg1.jpg'], 'jpg');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to deal with shadows, remove those edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pixels whose orientations are oriental or
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nearly oriental
orientationMapCopyNew = orientationMapDisplay;
nonHorizontal_edges = ones(imCopyH,imCopyW);
% nonHorizontal_edges((orientationMapCopyNew<=180)&(orientationMapCopyNew>=165))=0;
nonHorizontal_edges((orientationMapCopyNew==180))=0;



%%%%%%%%%%%%%%% voting for vainishing point estimation %%%%%%%%%%%
edgeImg = ones(imCopyH,imCopyW).*outlierBinary.*nonHorizontal_edges;
[rowInd, colInd] = find(double(edgeImg)==1);
numOfEdgePixels = sum(sum(double(edgeImg)));

imwrite(uint8(orientationMap), [outputPath,int2str(numOfValidFiles),'vpOrientationMap.jpg'], 'jpg');
imwrite(uint8(orientationMap1), [outputPath,int2str(numOfValidFiles),'vpOrientationMap1.jpg'], 'jpg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
votingMap = zeros(imCopyH, imCopyW);
votingMap1 = zeros(imCopyH, imCopyW);
interval = 1;
uppPercent = 0.9;
borderWidth = halfKerSize;


half_largestDistence = largestDistence*0.5;

halfImgH = round(imCopyH*0.4);  %%% change to smaller values if you want to get better results for those very low vanishing points

for i=1:interval:round(imCopyH*0.85) 
   for j=borderWidth+1:interval:imCopyW-borderWidth
           for ind=1:numOfEdgePixels
              
               if (rowInd(ind)>i)&&(rowInd(ind)<i+halfImgH)
               tmpI = rowInd(ind)-i;
               tmpJ = colInd(ind)-j;
               tempDist = sqrt(tmpJ*tmpJ+tmpI*tmpI);
                   
                   alpha = acos(tmpJ/tempDist)*180/pi;
                   theta = orientationMapDisplay(rowInd(ind), colInd(ind)); %%orientationMap(ii, jj);
                   angleDiffer = abs(alpha - theta);
                   
                   theta1 = orientationMapDisplay1(rowInd(ind), colInd(ind)); %%orientationMap(ii, jj);
                   angleDiffer1 = abs(alpha - theta1);
                   

                   if angleDiffer<=angleInterval
                          votingMap(i,j) = votingMap(i,j) + exp(-tempDist*angleDiffer/half_largestDistence); 
                   end
                   
                   if angleDiffer1<=angleInterval
                          votingMap1(i,j) = votingMap1(i,j) + exp(-tempDist*angleDiffer1/half_largestDistence); 
                   end
               end
           end
   end
end

votingMap_dual = votingMap + votingMap1;

max_votingMap = max(max(votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
min_votingMap = min(min(votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
tmpVotingMap = (votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)-min_votingMap)*255/(max_votingMap-min_votingMap);
votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap;

max_votingMap = max(max(votingMap1(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
min_votingMap = min(min(votingMap1(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
tmpVotingMap = (votingMap1(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)-min_votingMap)*255/(max_votingMap-min_votingMap);
votingMap1(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap;

max_votingMap = max(max(votingMap_dual(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
min_votingMap = min(min(votingMap_dual(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
tmpVotingMap = (votingMap_dual(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)-min_votingMap)*255/(max_votingMap-min_votingMap);
votingMap_dual(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpVotingMap = votingMap;
tmpVotingMap(round(0.85*imCopyH):imCopyH,:)=0;
% tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth).*(1-candidateVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth));
votingMap = tmpVotingMap;

www=5;
sum_votingMap = zeros(imCopyH,imgW);
for iii=1+www:imCopyH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end
tmpVotingMap = sum_votingMap;

max_votingMap = max(max(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)));
[r,c] = find(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)==max_votingMap);

if length(r)>1
    r = round(mean(r));
    c = round(mean(c))+borderWidth;
else
    r = r;
    c = c+borderWidth;
end



votingMap = double(uint8(votingMap));

sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(votingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
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

vpX=c;
vpY=r;


doim = colorImgCopy;
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
    
votingMap(round(0.8*imCopyH):imCopyH,:)=0;
imwrite(uint8(votingMap), [outputPath,int2str(numOfValidFiles),'vp1VotingMap_ini.jpg'], 'jpg');
    
imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vp1Map_ini.jpg'], 'jpg');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
tmpVotingMap = votingMap1;
tmpVotingMap(round(0.85*imCopyH):imCopyH,:)=0;
% tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth).*(1-candidateVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth));
votingMap1 = tmpVotingMap;

www=5;
sum_votingMap = zeros(imCopyH,imgW);
for iii=1+www:imCopyH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end
tmpVotingMap = sum_votingMap;

max_votingMap = max(max(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)));
[r,c] = find(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)==max_votingMap);

if length(r)>1
    r = round(mean(r));
    c = round(mean(c))+borderWidth;
else
    r = r;
    c = c+borderWidth;
end


votingMap1 = double(uint8(votingMap1));

sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(votingMap1(iii-www:iii+www,jjj-www:jjj+www)));
        
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


vpX1=c;
vpY1=r;


doim = colorImgCopy;
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
    
    votingMap1(round(0.8*imCopyH):imCopyH,:)=0;
imwrite(uint8(votingMap1), [outputPath,int2str(numOfValidFiles),'vp1VotingMap_ini1.jpg'], 'jpg');

    
imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vp1Map_ini1.jpg'], 'jpg');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
tmpVotingMap = votingMap_dual;
tmpVotingMap(round(0.85*imCopyH):imCopyH,:)=0;
% tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = tmpVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth).*(1-candidateVotingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth));
votingMap_dual = tmpVotingMap;


www=5;
sum_votingMap = zeros(imCopyH,imgW);
for iii=1+www:imCopyH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(tmpVotingMap(iii-www:iii+www,jjj-www:jjj+www)));
        
    end
end
tmpVotingMap = sum_votingMap;

max_votingMap = max(max(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)));
[r,c] = find(tmpVotingMap(1:1:round(imCopyH*uppPercent),borderWidth+1:1:imgW-borderWidth)==max_votingMap);

if length(r)>1
    r = round(mean(r));
    c = round(mean(c))+borderWidth;
else
    r = r;
    c = c+borderWidth;
end



votingMap_dual = double(uint8(votingMap_dual));

sum_votingMap = zeros(imgH,imgW);
for iii=1+www:imgH-www
    for jjj=1+www:imgW-www
        sum_votingMap(iii,jjj)=sum(sum(votingMap_dual(iii-www:iii+www,jjj-www:jjj+www)));
        
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

vpX_dual=c;
vpY_dual=r;


doim = colorImgCopy;
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
    
    votingMap_dual(round(0.8*imCopyH):imCopyH,:)=0;
imwrite(uint8(votingMap_dual), [outputPath,int2str(numOfValidFiles),'vp1VotingMap_ini_dual.jpg'], 'jpg');
imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vp1Map_ini_dual.jpg'], 'jpg');


