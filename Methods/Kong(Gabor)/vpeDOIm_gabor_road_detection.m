function [vpY, vpX] = vpeDOIm_gabor_road_detection(grayImg,colorImg, norient, outputPath_final, numOfValidFiles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grayImgCopy = grayImg;
[imCopyH,imCopyW] = size(grayImgCopy);

largestDistence = sqrt(imCopyH*imCopyH + imCopyW*imCopyW);

colorImgCopy = colorImg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            candidateVotingMap(rr-halftmpVerBarLen:rr+halftmpVerBarLen,tt-5:tt+5)=0; 
            outlierBinary(rr-halftmpVerBarLen:rr+halftmpVerBarLen,tt-5:tt+5)=0;
        end
    end
end


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

oddKernel = zeros(kerSize, kerSize, norient);
evenKernel = zeros(kerSize, kerSize, norient);
delta = kerSize/9;
tmpDelta = -1/(delta*delta*8);
c = 2*pi/lamda; %%%%
cc = 2.2; %%%%


cosTheta = zeros(angleRange, 1);
sinTheta = zeros(angleRange, 1);
for theta = 1:angleRange
    cosTheta(theta) = cos((theta-1)*pi/180);
    sinTheta(theta) = sin((theta-1)*pi/180);
end




%%%%%%%%%%%%%%%%%%%%%% kernels generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for theta = 90+1:angleInterval:angleRange-angleInterval+1+90
    tmpTheta = (theta - 1)*pi/180;
    for y= -halfKerSize:halfKerSize
        ySinTheta = y*sin(tmpTheta);
        yCosTheta = y*cos(tmpTheta);
        for x=-halfKerSize:halfKerSize
            xCosTheta = x*cos(tmpTheta);
            xSinTheta = x*sin(tmpTheta);
            a = xCosTheta+ySinTheta;
            b = -xSinTheta+yCosTheta;
            oddKernel(y+halfKerSize+1,x+halfKerSize+1,(theta - 1-90)/angleInterval+1) = exp(tmpDelta*(4*a*a+b*b))*(sin(c*a)-exp(-cc*cc/2));
            evenKernel(y+halfKerSize+1,x+halfKerSize+1,(theta - 1-90)/angleInterval+1) = exp(tmpDelta*(4*a*a+b*b))*(cos(c*a)-exp(-cc*cc/2));
        end
    end
end

%%%%%%%%%%%%%%%%% normalize kernels %%%%%%%%%%%%%%%%%%%%%
normalizedOddKernel = zeros(kerSize, kerSize, norient);
normalizedEvenKernel = zeros(kerSize, kerSize, norient);
for i=1:norient
    tmpKernel = oddKernel(:,:,i)-mean(mean(oddKernel(:,:,i)));
    tmpKernel = tmpKernel/(norm(tmpKernel));
    normalizedOddKernel(:,:,i) = tmpKernel;
    
    tmpKernel = evenKernel(:,:,i)-mean(mean(evenKernel(:,:,i)));
    tmpKernel = tmpKernel/(norm(tmpKernel));
    normalizedEvenKernel(:,:,i) = tmpKernel;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% image convolution with gabor filter banks %%%%%%%%%%%%%

filteredImgsOdd = zeros(imCopyH,imCopyW,norient*1);
filteredImgsEven = zeros(imCopyH,imCopyW,norient*1);
complexResponse = zeros(imCopyH,imCopyW,norient*1);
for i=1:norient
    filteredImgsOdd(:,:,i) = conv2(double(grayImg), normalizedOddKernel(:,:,i), 'same');
    filteredImgsEven(:,:,i) = conv2(double(grayImg), normalizedEvenKernel(:,:,i), 'same');
    complexResponse(:,:,i) = filteredImgsOdd(:,:,i).*filteredImgsOdd(:,:,i) + filteredImgsEven(:,:,i).*filteredImgsEven(:,:,i);
    complexResponse(:,:,i) = complexResponse(:,:,i)/(kerSize*kerSize);
end


%%%%%%%%%%%%%%% the dominant orientation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
confidenceMap = zeros(imCopyH,imCopyW);

orientationMap = zeros(imCopyH,imCopyW);
for i=1+halfKerSize:imCopyH-halfKerSize
    for j=1+halfKerSize:imCopyW-halfKerSize
        maxV = max(complexResponse(i,j,:));
        minV = min(complexResponse(i,j,:));
        complexResponseVector = (complexResponse(i,j,:)-minV)*100/(maxV-minV);
        maxLoc = find(complexResponseVector==100);
        if length(maxLoc)>0
            maxLoc = round(mean(maxLoc));
            if (maxLoc>2)&&(maxLoc<=34)
                complexResponseVector(maxLoc-2:maxLoc+2)=0;
            elseif (maxLoc<=2)
                complexResponseVector(maxLoc:maxLoc+2)=0;
            elseif (maxLoc>=35)
                complexResponseVector(maxLoc-2:maxLoc)=0;
            end

            [a,b] = sort(complexResponseVector,'descend');
            tmpContrast = (100 - mean(a(5:15)))^2;

            
            if maxV>1
                confidenceMap(i,j) = tmpContrast;
            end
        end
        
        maxma = max(complexResponse(i,j,:));
        indx = find(complexResponse(i,j,:)==maxma);
        orientationMap(i,j) = mean(indx); %%indx(1);
    end
end
orientationMap = (orientationMap-1)*angleInterval;



maxV = max(max(confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize)));
minV = min(min(confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize)));
aaTmp = confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize);
aaTmp = (aaTmp-minV)*255/(maxV-minV);
confidenceMap(1+halfKerSize:imCopyH-halfKerSize,1+halfKerSize:imCopyW-halfKerSize) = aaTmp;
imwrite(uint8(confidenceMap),[outputPath, int2str(numOfValidFiles),'vpConfidence.jpg'], 'jpg');


confidenceMapBinary = (confidenceMap>30);

doim = zeros(imCopyH, imCopyW, 3);
doim(:,:,1) = grayImgCopy;
doim(:,:,2) = grayImgCopy;
doim(:,:,3) = grayImgCopy;
doim(:,:,1) = doim(:,:,1).*(1-confidenceMapBinary) + confidenceMapBinary*255;
doim(:,:,2) = doim(:,:,2).*(1-confidenceMapBinary);
doim(:,:,3) = doim(:,:,3).*(1-confidenceMapBinary);
imwrite(uint8(doim),[outputPath, int2str(numOfValidFiles),'vpConfidenceOverlap.jpg'], 'jpg');


%%%%%%%%%%%%%%% show the orientation bar on resized image %%%%%%%%%%%%%%%%%%%%%%%
orientationMapDisplay = orientationMap;
doim = zeros(imCopyH, imCopyW, 3);
doim(:,:,1) = grayImgCopy;
doim(:,:,2) = grayImgCopy;
doim(:,:,3) = grayImgCopy;
doim_resized = imresize(doim,2,'bilinear');
orientationMapDisplay_resized = imresize(orientationMapDisplay,2,'bilinear'); 

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
        
        
            
    end
end

imwrite(uint8(doim_resized), [outputPath,int2str(numOfValidFiles),'vpOrientationBarImg.jpg'], 'jpg');




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




votingMap = zeros(imCopyH, imCopyW);
interval = 1;
uppPercent = 0.9;
borderWidth = halfKerSize;



half_largestDistence = largestDistence*0.5;

halfImgH = round(imCopyH*0.4);  %%% change to smaller values if you want to get better results for those extremely low vanishing points

for i=1:interval:round(imCopyH*0.9) 
   for j=borderWidth+1:interval:imCopyW-borderWidth
           for ind=1:numOfEdgePixels
              
               if (rowInd(ind)>i)&&(rowInd(ind)<i+halfImgH)
               tmpI = rowInd(ind)-i;
               tmpJ = colInd(ind)-j;
               tempDist = sqrt(tmpJ*tmpJ+tmpI*tmpI);
                   
                   alpha = acos(tmpJ/tempDist)*180/pi;
                   theta = orientationMapDisplay(rowInd(ind), colInd(ind)); %%orientationMap(ii, jj);
                   angleDiffer = abs(alpha - theta);
                   

                   if angleDiffer<=angleInterval
                          votingMap(i,j) = votingMap(i,j) + exp(-tempDist*angleDiffer/half_largestDistence); 
                   end
               end
           end
   end
end



max_votingMap = max(max(votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
min_votingMap = min(min(votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
votingMap1 = (votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth)-min_votingMap)*255/(max_votingMap-min_votingMap);
votingMap(1:interval:round(imCopyH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = votingMap1;





tmpVotingMap = votingMap;
tmpVotingMap(round(0.9*imCopyH):imCopyH,:)=0;
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
    
    imwrite(uint8(votingMap), [outputPath,int2str(numOfValidFiles),'vp1VotingMap_ini.jpg'], 'jpg');
    
imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vp1Map_ini.jpg'], 'jpg');




























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%     thr = adaptiveThreholding1(votingMap(1:round(uppPercent*imCopyH),:),1);
%     %%imwrite(uint8(segImg), [outputPath,int2str(numOfValidFiles),'vp1VotingMap_ini_seg.jpg'], 'jpg');
% 
%     segImg = (votingMap>thr+30);
%     
% 
%             %%%% get the largest connected area
%         [clusterLabels, numOfAllClusters] = bwlabel(segImg);
%            if numOfAllClusters>=1  %%% always larger than or equal to 1
%                allClustersData = regionprops(clusterLabels,'basic');
%                clusterArea = 0;
%                clusterIndex = 0;
%                alargestSize = 0;
%                for xx=1:numOfAllClusters
%                    tmpArea = allClustersData(xx).Area;
%                    if tmpArea>alargestSize
%                        clusterIndex = xx;
%                        alargestSize = tmpArea;
%                    end
%                end
%            end
%        segImg = (clusterLabels==clusterIndex);
% 
%     
%      
%     spikeLength = length(find(sum(segImg)>0));
%     spikeVector = find(sum(segImg)>0);
%     aTempVec = zeros(1,spikeLength);
%     zeroLength = 0;
%     for iii = 1:spikeLength
%         aTempMinVal = min(find(segImg(:,spikeVector(iii))==1));
%         if aTempMinVal<=3
%             zeroLength = zeroLength + 1;
%         end
%     end
% 
% 
%     sum_segImg = sum(segImg);
%     tta = size(find(sum_segImg<5),2);
%     if tta>=(0.5*imgW)||(zeroLength/spikeLength>=0.75)  %%deviate to left or right
%         [clusterLabels, numOfAllClusters] = bwlabel(segImg);
%            if numOfAllClusters>=1  %%% always larger than or equal to 1
%                allClustersData = regionprops(clusterLabels,'basic');
%                clusterArea = 0;
%                clusterIndex = 0;
%                alargestSize = 0;
%                for xx=1:numOfAllClusters
%                    tmpArea = allClustersData(xx).Area;
%                    if tmpArea>alargestSize
%                        clusterIndex = xx;
%                        alargestSize = tmpArea;
%                    end
%                end
%            end
%        segImg = (clusterLabels==clusterIndex);
%        aTempImg = segImg;
%        
%        
       
       
       
       
       
       
       
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           votingArea = zeros(imCopyH,imCopyW);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%% voting area
           
            


            %%%%% searching for the most possible road (or road border) direction
            tmpVal = zeros(31,1);
            tempDist1 = zeros(31,1);
            

            angleInterval = 5;
            
            for i=3:33
                    tmpAngle = (i*angleInterval)*pi/180;
                    tmpAngle1 = i*angleInterval;
                    tmpCounter = 0;
                    tmpCounter1 = 0;
                    if tmpAngle1~=90

                        if tmpAngle1<90
                            yJoint = round(tan(tmpAngle)*(imCopyW-c)+r);
                            if yJoint<=imCopyH
                                tempDist1(i-2) = sqrt((imCopyW-c)*(imCopyW-c) + (yJoint-r)*(yJoint-r));
                            else
                                xJoint = (imCopyH-1 - r)/(tan(tmpAngle)) + c;
                                tempDist1(i-2) = sqrt((xJoint-c)*(xJoint-c) + (imCopyH-r)*(imCopyH-r));
                            end
                        else
                            yJoint = round(tan(tmpAngle)*(1-c)+r);
                            if yJoint<=imCopyH
                                tempDist1(i-2) = sqrt((1-c)*(1-c) + (yJoint-r)*(yJoint-r));
                            else
                                xJoint = (imCopyH-1 - r)/(tan(tmpAngle)) + c;
                                tempDist1(i-2) = sqrt((xJoint-c)*(xJoint-c) + (imCopyH-r)*(imCopyH-r));
                            end
                        end

                        for x=1:imCopyW
                            y = round(tan(tmpAngle)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)
                                tmpCounter1 = tmpCounter1+1;
                                ori = orientationMap(y,x);
                                if abs(ori-tmpAngle1)<=angleInterval*2;
                                    tmpCounter = tmpCounter+1;
                                end

                            end
                        end

                    else
                        tmpCounter1 = imCopyH-r;
                        tempDist1(i-2) = tmpCounter1;
                        for y=r:imCopyH
                            ori = orientationMap(y,c);
                            if abs(ori-90)<=angleInterval*2
                                tmpCounter = tmpCounter+1;
                            end
                        end

                    end

                    tmpVal(i-2) = tmpCounter/tmpCounter1;
                    
            end
            
            %%%% select the top 8 edges
            tempDist1_binary = double(tempDist1>=imCopyH*0.5);
            tmpVal = tmpVal.*tempDist1_binary;
            [sorted_tmpVal, order_tmpVal] = sort(tmpVal,'descend');
            tmpVal1 = tmpVal(order_tmpVal(1:8));
            tmpVal1 = tmpVal1/(max(tmpVal1));
            
            tmpTopAngles = 10 + order_tmpVal*5;
            [forClustering_sort,forClustering_order] = sort(tmpTopAngles,'descend');
            %%% judge how many clusters and find the largetst cluster
                            finalNum = 8;
                            clusterNum = 1;
                            
                            clusterSeg = [];
                                for xx=1:finalNum-1
                                    if abs(forClustering_sort(xx)-forClustering_sort(xx+1))>5
                                        clusterNum = clusterNum+1;
                                        clusterSeg = [clusterSeg xx];
                                    end
                                end
                                
                                numOfLinesWithinEachCluster = zeros(clusterNum,1);
                                if length(clusterSeg)>0
                                    numOfLinesWithinEachCluster(1) = clusterSeg(1);
                                    for xx=2:clusterNum-1
                                        numOfLinesWithinEachCluster(xx) = clusterSeg(xx)-clusterSeg(xx-1);
                                    end
                                    numOfLinesWithinEachCluster(clusterNum) = finalNum - clusterSeg(clusterNum-1);
                                    maxClusterIndex = find(numOfLinesWithinEachCluster==(max(numOfLinesWithinEachCluster)));
                                else
                                    clusterSeg = finalNum;
                                    maxClusterIndex = 1;
                                end
                                
                                if length(maxClusterIndex)==1
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                elseif length(maxClusterIndex)==2
                                    maxClusterIndex = maxClusterIndex(2);
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                else
                                    maxClusterIndex = maxClusterIndex(length(maxClusterIndex));
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                end
                                
                                
            sbswRatio = zeros(length(largestCluster),1);
            
            for i=1:length(largestCluster)
                emptyImg = zeros(imCopyH,imCopyW);
                angleCC = largestCluster(i); %%10 + order_tmpVal(i)*5;
                angleCCRad = angleCC*pi/180;
                if angleCC<60   
                    if angleCC<30
                        angleCC=30;
                    end
                    tmpAngle1 = angleCC-30;
                    tmpAngle2 = angleCC+30;
                    tmpAngle1Rad = tmpAngle1*pi/180;
                    tmpAngle2Rad = tmpAngle2*pi/180;
                    
                    tmp_y = round(tan(angleCCRad)*(imCopyW-c)+r);
                    tmp_y1 = round(tan(tmpAngle1Rad)*(imCopyW-c)+r);
                    tmp_y2 = round(tan(tmpAngle2Rad)*(imCopyW-c)+r);
                    if (tmp_y1<=imCopyH)&&(tmp_y2<=imCopyH)
                        for x=c:imCopyW
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                emptyImg(y1:y,x)=1;
                                emptyImg(y:y2,x)=2;
                            end
                        end
                    elseif (tmp_y1<=imCopyH)&&(tmp_y2>imCopyH)
                        tmp_x2 = round((imCopyH-r)/tan(tmpAngle2Rad)+c);
                        if tmp_x2>imCopyW 
                            tmp_x2=imCopyW; 
                        end
                        if tmp_y<=imCopyH
                            for x=c:imCopyW
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                    emptyImg(y1:y,x)=1;
                                end
                            end
                            for x=c:tmp_x2
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                    emptyImg(y:y2,x)=2;
                                end
                            end
                            for x=tmp_x2:imCopyW
                                y=round(tan(angleCCRad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)
                                    emptyImg(y:imCopyH,x)=2;
                                end
                            end
                        else
                            tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                            if tmp_x>imCopyW 
                                tmp_x=imCopyW; 
                            end
                            for x=c:tmp_x
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                    emptyImg(y1:y,x)=1;
                                end
                            end
                            for x=tmp_x:imCopyW
                                y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                if (y1>=r)&&(y1<=imCopyH)
                                    emptyImg(y1:imCopyH,x)=1;
                                end
                            end
                            for x=c:tmp_x2
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                    emptyImg(y:y2,x)=2;
                                end
                            end
                            for x=tmp_x2:tmp_x
                                y=round(tan(angleCCRad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)
                                    emptyImg(y:imCopyH,x)=2;
                                end
                            end
                        end
                    elseif (tmp_y1>imCopyH)&&(tmp_y2>imCopyH)
                        tmp_x1 = round((imCopyH-r)/tan(tmpAngle1Rad)+c);
                        tmp_x2 = round((imCopyH-r)/tan(tmpAngle2Rad)+c);
                        tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                        if tmp_x>imCopyW 
                            tmp_x=imCopyW; 
                        end
                        if tmp_x1>imCopyW 
                            tmp_x1=imCopyW; 
                        end
                        if tmp_x2>imCopyW 
                            tmp_x2=imCopyW; 
                        end
                        for x=c:tmp_x
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                emptyImg(y1:y,x)=1;
                            end
                        end
                        for x=tmp_x:tmp_x1
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            if (y1>=r)&&(y1<=imCopyH)
                                emptyImg(y1:imCopyH,x)=1;
                            end
                        end
                        for x=c:tmp_x2
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                emptyImg(y:y2,x)=2;
                            end
                        end
                        for x=tmp_x2:tmp_x
                            y=round(tan(angleCCRad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)
                                emptyImg(y:imCopyH,x)=2;
                            end
                        end
                        
                    end
                elseif angleCC>120
                    if angleCC>150
                        angleCC=150;
                    end
                    tmpAngle1 = angleCC-30;
                    tmpAngle2 = angleCC+30;
                    tmpAngle1Rad = tmpAngle1*pi/180;
                    tmpAngle2Rad = tmpAngle2*pi/180;
                    
                    tmp_y = round(tan(angleCCRad)*(1-c)+r);
                    tmp_y1 = round(tan(tmpAngle1Rad)*(1-c)+r);
                    tmp_y2 = round(tan(tmpAngle2Rad)*(1-c)+r);
                    if (tmp_y1<=imCopyH)&&(tmp_y2<=imCopyH)
                        for x=1:c
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                emptyImg(y2:y,x)=2;
                                emptyImg(y:y1,x)=1;
                            end
                        end
                    elseif (tmp_y1>imCopyH)&&(tmp_y2<=imCopyH)
                        tmp_x1 = round((imCopyH-r)/tan(tmpAngle1Rad)+c);
                        if tmp_x1<1 
                            tmp_x1=1; 
                        end
                        if tmp_y<=imCopyH
                            for x=1:c
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                    emptyImg(y2:y,x)=2;
                                end
                            end
                            for x=1:tmp_x1
                                y=round(tan(angleCCRad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)
                                    emptyImg(y:imCopyH,x)=1;
                                end
                            end
                            for x=tmp_x1:c
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                    emptyImg(y:y1,x)=1;
                                end
                            end
                        else
                            tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                            if tmp_x1<1 
                                tmp_x1=1; 
                            end
                            for x=1:tmp_x
                                y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                if (y2>=r)&&(y2<=imCopyH)
                                    emptyImg(y2:imCopyH,x)=2;
                                end
                            end
                            for x=tmp_x:c
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                    emptyImg(y2:y,x)=2;
                                end
                            end
                            for x=tmp_x:tmp_x1
                                y=round(tan(angleCCRad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)
                                    emptyImg(y:imCopyH,x)=1;
                                end
                            end
                            for x=tmp_x1:c
                                y=round(tan(angleCCRad)*(x-c)+r);
                                y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                    emptyImg(y:y1,x)=1;
                                end
                            end
                            
                        end
                    elseif (tmp_y1>imCopyH)&&(tmp_y2>imCopyH)
                        tmp_x1 = round((imCopyH-r)/tan(tmpAngle1Rad)+c);
                        tmp_x2 = round((imCopyH-r)/tan(tmpAngle2Rad)+c);
                        tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                        if tmp_x1<1 
                            tmp_x1=1; 
                        end
                        if tmp_x2<1 
                            tmp_x2=1; 
                        end
                        if tmp_x<1 
                            tmp_x=1; 
                        end
                        for x=tmp_x2:tmp_x
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y2>=r)&&(y2<=imCopyH)
                                emptyImg(y2:imCopyH,x)=2;
                            end
                        end
                        for x=tmp_x:c
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                emptyImg(y2:y,x)=2;
                            end
                        end
                        for x=tmp_x:tmp_x1
                            y=round(tan(angleCCRad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)
                                emptyImg(y:imCopyH,x)=1;
                            end
                        end
                        for x=tmp_x1:c
                            y=round(tan(angleCCRad)*(x-c)+r);
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                emptyImg(y:y1,x)=1;
                            end
                        end
                    end
                elseif (angleCC<=120)&&(angleCC>=60)
                    tmpAngle1 = angleCC-30;
                    tmpAngle2 = angleCC+30;
                    tmpAngle1Rad = tmpAngle1*pi/180;
                    tmpAngle2Rad = tmpAngle2*pi/180;
                    
%                     tmp_y = round(tan(angleCCRad)*(1-c)+r);
%                     tmp_y1 = round(tan(tmpAngle1Rad)*(1-c)+r);
%                     tmp_y2 = round(tan(tmpAngle2Rad)*(1-c)+r);
                    if angleCC==90
                        for x=1:c
                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                            if (y2>=r)&&(y2<=imCopyH)
                                emptyImg(y2:imCopyH,x)=2;
                            end
                        end
                        for x=c:imCopyW
                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                            if (y1>=r)&&(y1<=imCopyH)
                                emptyImg(y1:imCopyH,x)=1;
                            end
                        end
                    else
                        if angleCC<90
                            tmp_y = round(tan(angleCCRad)*(imCopyW-c)+r);
                            if tmpAngle2~=90
                                tmp_y1 = round(tan(tmpAngle1Rad)*(imCopyW-c)+r);
                                tmp_y2 = round(tan(tmpAngle2Rad)*(1-c)+r);
                                for x=1:c
                                    y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                    if (y2>=r)&&(y2<imCopyH)
                                        emptyImg(y2:imCopyH,x)=2;
                                    end
                                end
                                for x=c:imCopyW
                                    y=round(tan(angleCCRad)*(x-c)+r);
                                    if (y>=r)&&(y<imCopyH)
                                        emptyImg(y:imCopyH,x)=2;
                                    end
                                end
                                if tmp_y<=imCopyH
                                    for x=c:imCopyW
                                        y=round(tan(angleCCRad)*(x-c)+r);
                                        y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                        if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                            emptyImg(y1:y,x)=1;
                                        end
                                    end
                                else
                                    tmp_x1 = round((imCopyH-r)/tan(tmpAngle1Rad)+c);
                                    tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                                    if tmp_y1<=imCopyH
                                        for x=c:tmp_x
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:y,x)=1;
                                            end
                                        end
                                        for x=tmp_x:imCopyW
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=1;
                                            end
                                        end
                                    else
                                        for x=c:tmp_x
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:y,x)=1;
                                            end
                                        end
                                        for x=tmp_x:tmp_x1
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=1;
                                            end
                                        end
                                    end
                                end
                            else
                                tmp_y1 = round(tan(tmpAngle1Rad)*(imCopyW-c)+r);
                                for x=c:imCopyW
                                    y=round(tan(angleCCRad)*(x-c)+r);
                                    if (y>=r)&&(y<imCopyH)
                                        emptyImg(y:imCopyH,x)=2;
                                    end
                                end
                                if tmp_y<=imCopyH
                                    for x=c:imCopyW
                                        y=round(tan(angleCCRad)*(x-c)+r);
                                        y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                        if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                            emptyImg(y1:y,x)=1;
                                        end
                                    end
                                else
                                    tmp_x1 = round((imCopyH-r)/tan(tmpAngle1Rad)+c);
                                    tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                                    if tmp_y1<=imCopyH
                                        for x=c:tmp_x
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:y,x)=1;
                                            end
                                        end
                                        for x=tmp_x:imCopyW
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=1;
                                            end
                                        end
                                    else
                                        for x=c:tmp_x
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:y,x)=1;
                                            end
                                        end
                                        for x=tmp_x:tmp_x1
                                            y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                            if (y1>=r)&&(y1<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=1;
                                            end
                                        end
                                    end
                                end
                            end
                        elseif angleCC>90
                            tmp_y = round(tan(angleCCRad)*(1-c)+r);
                            if tmpAngle1~=90
                                tmp_y1 = round(tan(tmpAngle1Rad)*(imCopyW-c)+r);
                                tmp_y2 = round(tan(tmpAngle2Rad)*(1-c)+r);
                                for x=c:imCopyW
                                    y1=round(tan(tmpAngle1Rad)*(x-c)+r);
                                    if (y1>=r)&&(y1<imCopyH)
                                        emptyImg(y1:imCopyH,x)=1;
                                    end
                                end
                                for x=1:c
                                    y=round(tan(angleCCRad)*(x-c)+r);
                                    if (y>=r)&&(y<imCopyH)
                                        emptyImg(y:imCopyH,x)=1;
                                    end
                                end
                                if tmp_y<=imCopyH
                                    for x=1:c
                                        y=round(tan(angleCCRad)*(x-c)+r);
                                        y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                        if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                            emptyImg(y2:y,x)=2;
                                        end
                                    end
                                else
                                    tmp_x2 = round((imCopyH-r)/tan(tmpAngle2Rad)+c);
                                    tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                                    if tmp_y2<=imCopyH
                                        for x=tmp_x:c
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:y,x)=2;
                                            end
                                        end
                                        for x=1:tmp_x
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=2;
                                            end
                                        end
                                    else
                                        for x=tmp_x:c
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:y,x)=2;
                                            end
                                        end
                                        for x=tmp_x2:tmp_x
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:imCopyH,x)=2;
                                            end
                                        end
                                    end
                                end
                            else
                                tmp_y2 = round(tan(tmpAngle2Rad)*(1-c)+r);
                                for x=1:c
                                    y=round(tan(angleCCRad)*(x-c)+r);
                                    if (y>=r)&&(y<imCopyH)
                                        emptyImg(y:imCopyH,x)=1;
                                    end
                                end
                                if tmp_y<=imCopyH
                                    for x=1:c
                                        y=round(tan(angleCCRad)*(x-c)+r);
                                        y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                        if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                            emptyImg(y2:y,x)=2;
                                        end
                                    end
                                else
                                    tmp_x2 = round((imCopyH-r)/tan(tmpAngle2Rad)+c);
                                    tmp_x = round((imCopyH-r)/tan(angleCCRad)+c);
                                    if tmp_y2<=imCopyH
                                        for x=tmp_x:c
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:y,x)=2;
                                            end
                                        end
                                        for x=1:tmp_x
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y1:imCopyH,x)=2;
                                            end
                                        end
                                    else
                                        for x=tmp_x:c
                                            y=round(tan(angleCCRad)*(x-c)+r);
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y>=r)&&(y<=imCopyH)&&(y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:y,x)=2;
                                            end
                                        end
                                        for x=tmp_x2:tmp_x
                                            y2=round(tan(tmpAngle2Rad)*(x-c)+r);
                                            if (y2>=r)&&(y2<=imCopyH)
                                                emptyImg(y2:imCopyH,x)=2;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                size1 = sum(sum(double(emptyImg==1)));
                size2 = sum(sum(double(emptyImg==2)));
                meanR1 = mean(mean((double(emptyImg==1).*double(colorImgCopy(:,:,1)))));
                meanR2 = mean(mean((double(emptyImg==2).*double(colorImgCopy(:,:,1)))));
                tmpR1 = (double(emptyImg==1).*double(colorImgCopy(:,:,1))) + (double(emptyImg~=1)*meanR1);
                varR1 = sqrt(sum(sum((tmpR1-meanR1).*(tmpR1-meanR1)))/size1);
                tmpR2 = (double(emptyImg==2).*double(colorImgCopy(:,:,1))) + (double(emptyImg~=2)*meanR2);
                varR2 = sqrt(sum(sum((tmpR2-meanR2).*(tmpR2-meanR2)))/size1);
                Sw = min(varR1,varR2);%%(varR1+varR2)/2;
                Sb = sqrt((meanR1-meanR2)*(meanR1-meanR2));
                tmpSbSwRatio1 = Sb;%%/Sw;
                meanG1 = mean(mean((double(emptyImg==1).*double(colorImgCopy(:,:,2)))));
                meanG2 = mean(mean((double(emptyImg==2).*double(colorImgCopy(:,:,2)))));
                tmpG1 = (double(emptyImg==1).*double(colorImgCopy(:,:,2))) + (double(emptyImg~=1)*meanG1);
                varG1 = sqrt(sum(sum((tmpG1-meanG1).*(tmpG1-meanG1)))/size1);
                tmpG2 = (double(emptyImg==2).*double(colorImgCopy(:,:,2))) + (double(emptyImg~=2)*meanG2);
                varG2 = sqrt(sum(sum((tmpG2-meanG2).*(tmpG2-meanG2)))/size1);
                Sw = min(varG1,varG2); %%(varG1+varG2)/2;
                Sb = sqrt((meanG1-meanG2)*(meanG1-meanG2));
                tmpSbSwRatio2 = Sb;%%/Sw;
                meanB1 = mean(mean((double(emptyImg==1).*double(colorImgCopy(:,:,3)))));
                meanB2 = mean(mean((double(emptyImg==2).*double(colorImgCopy(:,:,3)))));
                tmpB1 = (double(emptyImg==1).*double(colorImgCopy(:,:,3))) + (double(emptyImg~=1)*meanB1);
                varB1 = sqrt(sum(sum((tmpB1-meanB1).*(tmpB1-meanB1)))/size1);
                tmpB2 = (double(emptyImg==2).*double(colorImgCopy(:,:,3))) + (double(emptyImg~=2)*meanB2);
                varB2 = sqrt(sum(sum((tmpB2-meanB2).*(tmpB2-meanB2)))/size1);
                Sw = min(varB1,varB2); %%(varB1+varB2)/2;
                Sb = sqrt((meanB1-meanB2)*(meanB1-meanB2));
                tmpSbSwRatio3 = Sb;%%/Sw;
                sbswRatio(i) = max([tmpSbSwRatio1,tmpSbSwRatio2,tmpSbSwRatio3]);
                
                
                
            end
            
            tmpVal1;
            sbswRatio;
            
            sbswRatio = sbswRatio/(max(sbswRatio));
            finalS = sbswRatio.*(tmpVal((largestCluster-10)/5)/max(tmpVal((largestCluster-10)/5))); %%.*tmpVal1;
            maxFinalS = max(finalS);
            para = zeros(2,1);
            angleTheta1 = largestCluster(find(finalS==maxFinalS)); %%order_tmpVal(find(finalS==maxFinalS))*5 + 10;
            
            if length(angleTheta1)>0
            para(1) = tan(angleTheta1*pi/180);
            
            
            
%             para = zeros(2,1);
%             para(1) = tan((tmpIndx*5)*pi/180);
%             angleTheta1 = tmpIndx*5;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update the second border
            angleInter = 5;
            numOfAngleChoices = 29;
            tanVal = zeros(numOfAngleChoices,1);
            for iii=1:numOfAngleChoices
                tanVal(iii) = tan((15+iii*5)*pi/180);
            end
            sumOftempCounterRatio_max = 0;
            selectedVp = zeros(2,1);
            selected_constructedLine = zeros(numOfAngleChoices,1);
            selected_tempCounterRatio = zeros(numOfAngleChoices,1);
            selected_tempCounter0 = zeros(numOfAngleChoices,1);
            selected_tempCounter = zeros(numOfAngleChoices,1);
            if para(1)>=0
                % endPoint = max(find(sum(aTempImg)>0));
                jointOfImh = (imCopyH-r)/(para(1)) + c;
                jointOfImhCopy = min(jointOfImh, imCopyW);
                for tempFixedVp_X= c-2:c+2; %2:imCopyW-3 %%c-1:c+1 %%3:5:endPoint
                    %%%%%% find the second dominant edge
                    %%%%%% the angle range of the second dominant edge
                    jointOfImh = jointOfImhCopy;
                    tempFixedVp_Y = round(para(1)*(tempFixedVp_X-c)+r);
                    
                    if (tempFixedVp_Y>1)&&(tempFixedVp_Y<imCopyH)
                        jointOfImh = jointOfImh - tempFixedVp_X;
                        constructedLine = zeros(numOfAngleChoices,1);
                        tempCounterRatio_array = zeros(numOfAngleChoices,1);
                        tempCounter0_array = zeros(numOfAngleChoices,1);
                        tempCounter_array = zeros(numOfAngleChoices,1);
                        for iii=1:numOfAngleChoices
                            tempCounter = 0;
                            tempCounter0 = 0;
                            tempCurrentAngle = 15+iii*5;
                            if (tempCurrentAngle~=90)%%&&(abs(tempCurrentAngle-angleTheta1)>=35)
                                if tempCurrentAngle<90
                                    yJoint = round(tanVal(iii)*(imCopyW-tempFixedVp_X)+tempFixedVp_Y);
                                    if yJoint<=imCopyH
                                        tempDist1 = sqrt((imCopyW-tempFixedVp_X)*(imCopyW-tempFixedVp_X) + (yJoint-tempFixedVp_Y)*(yJoint-tempFixedVp_Y));
                                    else
                                        xJoint = (imCopyH-1 - tempFixedVp_Y)/(tanVal(iii)) + tempFixedVp_X;
                                        tempDist1 = sqrt((xJoint-tempFixedVp_X)*(xJoint-tempFixedVp_X) + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                                    end
                                else
                                    yJoint = round(tanVal(iii)*(1-tempFixedVp_X)+tempFixedVp_Y);
                                    if yJoint<=imCopyH
                                        tempDist1 = sqrt((1-tempFixedVp_X)*(1-tempFixedVp_X) + (yJoint-tempFixedVp_Y)*(yJoint-tempFixedVp_Y));
                                    else
                                        xJoint = (imCopyH-1 - tempFixedVp_Y)/(tanVal(iii)) + tempFixedVp_X;
                                        tempDist1 = sqrt((xJoint-tempFixedVp_X)*(xJoint-tempFixedVp_X) + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                                    end
                                end
                                for jjj=1:imCopyW
                                    y = round(tanVal(iii)*(jjj-tempFixedVp_X)+tempFixedVp_Y);
                                    if (y>tempFixedVp_Y)&&(y<=imCopyH)
                                        tempCounter = tempCounter + 1;
                                        ori = orientationMap(y,jjj);
                                        if abs(ori-tempCurrentAngle)<=(angleInter+angleInter)
                                            tempCounter0 = tempCounter0+1;
                                        end
                                    end
                                end
    %                             if tempCounter>0
    %                                 if (tempCounter>=0.5*r)&&((tempCounter0/tempCounter)>0.15)
    %                                     constructedLine(iii)=1;
    %                                 end
    %                             end
    %                             tempCounterRatio(iii) = tempCounter0/tempCounter;
                            elseif (tempCurrentAngle==90)%%&&(abs(tempCurrentAngle-angleTheta1)>=35)
                                tempCounter = imCopyH-tempFixedVp_Y;
                                tempDist1 = tempCounter;
                                for jjj=tempFixedVp_Y:imCopyH
                                    ori = orientationMap(jjj,tempFixedVp_X);
                                    if abs(ori-tempCurrentAngle)<=(angleInter+angleInter)
                                        tempCounter0 = tempCounter0+1;
                                    end
                                end
    %                             if (tempCounter>=0.5*r)&&((tempCounter0/tempCounter)>0.15)
    %                                 constructedLine(iii)=1;
    %                             end
    %                             tempCounterRatio(iii) = tempCounter0/tempCounter;
                            end
                            tempCounterRatio = tempCounter0/tempCounter;
                            tempDist = sqrt(jointOfImh*jointOfImh + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                            if (tempCounterRatio>0.02)&&(abs(tempCurrentAngle-angleTheta1)>=20)&&(tempDist>0.35*imCopyH)&&(tempDist1>0.35*imCopyH)%%(tempCounterRatio>sumOftempCounterRatio_max)&&
                                tempCounterRatio_array(iii) = tempCounterRatio;
                                constructedLine(iii) = 1;
                            end
                        end
                        [tempCounterRatio_array_sort, tempCounterRatio_array_order] = sort(tempCounterRatio_array,'descend');
                        if (sum(tempCounterRatio_array_sort(1:6))>sumOftempCounterRatio_max)
                            selected_constructedLine = tempCounterRatio_array_order;
                            selectedVp = [tempFixedVp_X,tempFixedVp_Y]';
                            selected_tempCounterRatio = tempCounterRatio_array;
                            sumOftempCounterRatio_max = sum(tempCounterRatio_array_sort(1:6));
                        end
                                                        
                                
%                         tempCounterRatio_array = (tempCounter0_array./tempCounter_array).*(tempCounter_array/sum(tempCounter_array));
%                         %%%% sort
%                         [tempCounterRatio_array_sort, tempCounterRatio_array_order] = sort(tempCounterRatio_array,'descend');
% %                         aTempAngleDiff = sum(abs(angleTheta1-(20+tempCounterRatio_array_order(1:1)*5)));
%                         if (sum(tempCounterRatio_array_sort(1:5))>sumOftempCounterRatio_max)%%&&(aTempAngleDiff>=35)
%                             selectedVp = [tempFixedVp_X,tempFixedVp_Y];
%                             selected_tempCounterRatio = tempCounterRatio_array;
%                             selected_tempCounter0 = tempCounter0_array;
%                             selected_tempCounter = tempCounter_array;
%                             sumOftempCounterRatio_max = sum(tempCounterRatio_array_sort(1:5));
%                         end
                    end
                end
                
            else
                % endPoint = min(find(sum(aTempImg)>0));
                jointOfImh = (imCopyH-r)/(para(1)) + c;
                jointOfImhCopy = max(jointOfImh, 1);
                for tempFixedVp_X=c+2:-1:c-2 %%c-1:c+1 %%imCopyW-3:-5:endPoint
                    %%%%%% find the second dominant edge
                    %%%%%% the angle range of the second dominant edge
                    jointOfImh = jointOfImhCopy;
                    tempFixedVp_Y = round(para(1)*(tempFixedVp_X-c)+r);
                    
                    if (tempFixedVp_Y>1)&&(tempFixedVp_Y<imCopyH)
                        jointOfImh = tempFixedVp_X-jointOfImh;
                        constructedLine = zeros(numOfAngleChoices,1);
                        tempCounterRatio_array = zeros(numOfAngleChoices,1);
                        tempCounter0_array = zeros(numOfAngleChoices,1);
                        tempCounter_array = zeros(numOfAngleChoices,1);
                        for iii=1:numOfAngleChoices
                            tempCounter = 0;
                            tempCounter0 = 0;
                            tempCurrentAngle = 15+iii*5;
                            if (tempCurrentAngle~=90)%%&&(abs(tempCurrentAngle-angleTheta1)>=35)
                                if tempCurrentAngle<90
                                    yJoint = round(tanVal(iii)*(imCopyW-tempFixedVp_X)+tempFixedVp_Y);
                                    if yJoint<=imCopyH
                                        tempDist1 = sqrt((imCopyW-tempFixedVp_X)*(imCopyW-tempFixedVp_X) + (yJoint-tempFixedVp_Y)*(yJoint-tempFixedVp_Y));
                                    else
                                        xJoint = (imCopyH-1 - tempFixedVp_Y)/(tanVal(iii)) + tempFixedVp_X;
                                        tempDist1 = sqrt((xJoint-tempFixedVp_X)*(xJoint-tempFixedVp_X) + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                                    end
                                else
                                    yJoint = round(tanVal(iii)*(1-tempFixedVp_X)+tempFixedVp_Y);
                                    if yJoint<=imCopyH
                                        tempDist1 = sqrt((1-tempFixedVp_X)*(1-tempFixedVp_X) + (yJoint-tempFixedVp_Y)*(yJoint-tempFixedVp_Y));
                                    else
                                        xJoint = (imCopyH-1 - tempFixedVp_Y)/(tanVal(iii)) + tempFixedVp_X;
                                        tempDist1 = sqrt((xJoint-tempFixedVp_X)*(xJoint-tempFixedVp_X) + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                                    end
                                end
                                for jjj=1:imCopyW
                                    y = round(tanVal(iii)*(jjj-tempFixedVp_X)+tempFixedVp_Y);
                                    if (y>tempFixedVp_Y)&&(y<=imCopyH)
                                        tempCounter = tempCounter + 1;
                                        ori = orientationMap(y,jjj);
                                        if abs(ori-tempCurrentAngle)<=(angleInter+angleInter)
                                            tempCounter0 = tempCounter0+1;
                                        end
                                    end
                                end
    %                             if tempCounter>0
    %                                 if (tempCounter>=0.5*r)&&((tempCounter0/tempCounter)>0.15)
    %                                     constructedLine(iii)=1;
    %                                 end
    %                             end
    %                             tempCounterRatio(iii) = tempCounter0/tempCounter;
                            elseif (tempCurrentAngle==90)%%&&(abs(tempCurrentAngle-angleTheta1)>=35)
                                tempCounter = imCopyH-tempFixedVp_Y;
                                tempDist1 = tempCounter;
                                for jjj=tempFixedVp_Y:imCopyH
                                    ori = orientationMap(jjj,tempFixedVp_X);
                                    if abs(ori-tempCurrentAngle)<=(angleInter+angleInter)
                                        tempCounter0 = tempCounter0+1;
                                    end
                                end
    %                             if (tempCounter>=0.5*r)&&((tempCounter0/tempCounter)>0.15)
    %                                 constructedLine(iii)=1;
    %                             end
    %                             tempCounterRatio(iii) = tempCounter0/tempCounter;
                            end
                            tempCounterRatio = tempCounter0/tempCounter;
                            tempDist = sqrt(jointOfImh*jointOfImh + (imCopyH-tempFixedVp_Y)*(imCopyH-tempFixedVp_Y));
                            if (tempCounterRatio>0.02)&&(abs(tempCurrentAngle-angleTheta1)>=20)&&(tempDist>0.35*imCopyH)&&(tempDist1>0.35*imCopyH)%%(tempCounterRatio>sumOftempCounterRatio_max)&&
                                tempCounterRatio_array(iii) = tempCounterRatio;
                                constructedLine(iii) = 1;
                            end
                        end
                        [tempCounterRatio_array_sort, tempCounterRatio_array_order] = sort(tempCounterRatio_array,'descend');
                        if (sum(tempCounterRatio_array_sort(1:6))>sumOftempCounterRatio_max)
                            selected_constructedLine = tempCounterRatio_array_order;
                            selectedVp = [tempFixedVp_X,tempFixedVp_Y]';
                            selected_tempCounterRatio = tempCounterRatio_array;
                            sumOftempCounterRatio_max = sum(tempCounterRatio_array_sort(1:6));
                        end

%                         tempCounterRatio_array = (tempCounter0_array./tempCounter_array).*(tempCounter_array/sum(tempCounter_array));
%                         %%%% sort
%                         [tempCounterRatio_array_sort, tempCounterRatio_array_order] = sort(tempCounterRatio_array,'descend');
% %                         aTempAngleDiff = sum(abs(angleTheta1-(20+tempCounterRatio_array_order(1:1)*5)));
%                         
%                         if (sum(tempCounterRatio_array_sort(1:5))>sumOftempCounterRatio_max)%%%&&(aTempAngleDiff>=35)
%                             selectedVp = [tempFixedVp_X,tempFixedVp_Y];
%                             selected_tempCounterRatio = tempCounterRatio_array;
%                             selected_tempCounter0 = tempCounter0_array;
%                             selected_tempCounter = tempCounter_array;
%                             sumOftempCounterRatio_max = sum(tempCounterRatio_array_sort(1:5));
%                         end
                    end
                    
                end
            end
    %             tempCounterRatio
    %             constructedLine

%             selected_tempCounterRatio
             
            
            if sum(selectedVp)>0
                newC = selectedVp(1);
                newR = selectedVp(2);
            else
                newC = round(0.5*imCopyW);
                newR = round(0.5*imCopyH);
%                 selected_constructedLine = 10;
            end
            
%             [bbb,ccc] = sort(selected_tempCounterRatio,'descend');
%             
%             selected_constructedLine = ccc(1:1);

   %%%%%%%%%%%%%%%%%%% plot the possible dominant edge   
%         doim = zeros(imCopyH, imCopyW, 3);
%         doim(:,:,1) = grayImgCopy;
%         doim(:,:,2) = grayImgCopy;
%         doim(:,:,3) = grayImgCopy;
          doim = colorImgCopy;
                        if newC==1
                            newC=2;
                        end
                        if newR==1
                            newR=2;
                        end
            if (newC>4)&&(newC<imCopyW-4);
            doim(newR:newR+4,newC-4:newC+4,1) = 255;
            doim(newR:newR+4,newC-4:newC+4,2) = 0;
            doim(newR:newR+4,newC-4:newC+4,3) = 0;
            end
            for x=1:imCopyW
                y = round(para(1)*(x-newC)+newR);
                if (y>=newR)&&(y<=imCopyH-1)
                    doim(y-1:y+1,x,1) = 255;
                    doim(y-1:y+1,x,2) = 0;
                    doim(y-1:y+1,x,3) = 0;

                end
            end
            
            if sum(selectedVp)>0
            for iii=1:8
                if selected_tempCounterRatio(selected_constructedLine(iii))>0
                if (para(1)>=0)&&(tanVal(selected_constructedLine(iii))>=0)

                    for x=newC:imCopyW
                        y1 = min(round(para(1)*(x-newC)+newR),imCopyH);
                        y2 = min(round(tanVal(selected_constructedLine(iii))*(x-newC)+newR),imCopyH); %%round(tanVal(iii)*(x-newC)+newR);
                        if (y2>=newR)&&(y2<=imCopyH-1)
                            doim(y2-1:y2+1,x,1) = 0;
                            doim(y2-1:y2+1,x,2) = 255;
                            doim(y2-1:y2+1,x,3) = 0;
                        end
                        if (y1>=newR)&&(y1<=imCopyH)&&(y2>=newR)&&(y2<=imCopyH)
                            doim(y2-1:y2+1,x,1) = 0;
                            doim(y2-1:y2+1,x,2) = 255;
                            doim(y2-1:y2+1,x,3) = 0;
                            if y1>=y2
                                y1=min(y1+10,imCopyH);
                                y2 = max(y2-10,1);
                                votingArea(y2:y1,x)=1;
                            else
                                y2=min(y2+10,imCopyH);
                                y1 = max(y1-10,1);
                                votingArea(y1:y2,x)=1;
                            end
                        end
                    end
                elseif (para(1)>=0)&&(tanVal(selected_constructedLine(iii))<=0)
                    for x=1:newC
                        y2 = min(round(tanVal(selected_constructedLine(iii))*(x-newC)+newR),imCopyH); 
                        if (y2>=newR)&&(y2<=imCopyH-1)
                            y2a = max(y2-10,2);
                            doim(y2-1:y2+1,x,1) = 0;
                            doim(y2-1:y2+1,x,2) = 255;
                            doim(y2-1:y2+1,x,3) = 0;
                            votingArea(y2a:imCopyH,x)=1;
                        end
                    end
                    for x=newC:imCopyW
                        y1 = min(round(para(1)*(x-newC)+newR),imCopyH-1);
                        if (y1>=newR)&&(y1<=imCopyH-1)
                            y1 = max(y1-10,1);
                            votingArea(y1:imCopyH,x)=1;
                        end
                    end
                elseif (para(1)<=0)&&(tanVal(selected_constructedLine(iii))<=0)
                    for x=1:newC
                        y1 = min(round(para(1)*(x-newC)+newR),imCopyH);
                        y2 = min(round(tanVal(selected_constructedLine(iii))*(x-newC)+newR),imCopyH); %%round(tanVal(iii)*(x-newC)+newR);
                        if (y2>=newR)&&(y2<=imCopyH-1)
                            doim(y2-1:y2+1,x,1) = 0;
                            doim(y2-1:y2+1,x,2) = 255;
                            doim(y2-1:y2+1,x,3) = 0;
                        end
                        if (y1>=newR)&&(y1<=imCopyH)&&(y2>=newR)&&(y2<=imCopyH)
                            if y1>=y2
                                y1=min(y1+10,imCopyH);
                                y2 = max(y2-10,1);
                                votingArea(y2:y1,x)=1;
                            else
                                y2=min(y2+10,imCopyH);
                                y1 = max(y1-10,1);
                                votingArea(y1:y2,x)=1;
                            end
                        end
                    end
                elseif (para(1)<=0)&&(tanVal(selected_constructedLine(iii))>=0)
                    for x=1:newC
                        y1 = min(round(para(1)*(x-newC)+newR),imCopyH-1);
                        if (y1>=newR)&&(y1<=imCopyH-1)
                            y1 = max(y1-10,1);
                            votingArea(y1:imCopyH,x)=1;
                        end
                    end
                    for x=newC:imCopyW
                        y2 = min(round(tanVal(selected_constructedLine(iii))*(x-newC)+newR),imCopyH); 
                        if (y2>=newR)&&(y2<=imCopyH-1)
                            y2a = max(y2-10,2);
                            doim(y2-1:y2+1,x,1) = 0;
                            doim(y2-1:y2+1,x,2) = 255;
                            doim(y2-1:y2+1,x,3) = 0;
                            votingArea(y2a:imCopyH,x)=1;
                        end
                    end
                end
                end
            end
            end
            %%%% update the votingArea 
            imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vp1Map_Temp1.jpg'], 'jpg'); 
            
            
            
            
            
            
            
            
            
            %%%%%% rough road area 
                        
                        numOfAnglesLargerThan_angleTheta1=0;
                        numOfAnglesSmallerThan_angleTheta1=0;
                        
                        numOfAnglesLargerThan_90=0;
                        numOfAnglesSmallerThan_90=0;
                        
                        largestPossibleAngle = 0;
                        smallestPossibleAngle = 180;
                        
                        
                        doim = colorImgCopy;
                        doim1 = colorImgCopy;
                        roadBinaryImage = zeros(imCopyH, imCopyW);
                        
                        if sum(selectedVp)>0
                            numOfNonZeroRatios = sum(selected_tempCounterRatio>0);
                            numOfFinalLines = min(numOfNonZeroRatios,8);
                            forClustering = zeros(numOfFinalLines,1);
                            for iii=1:numOfFinalLines
                                if selected_tempCounterRatio(selected_constructedLine(iii))>0
                                    if (15+selected_constructedLine(iii)*5)>angleTheta1
                                        numOfAnglesLargerThan_angleTheta1 = numOfAnglesLargerThan_angleTheta1 + 1;
                                        if (15+selected_constructedLine(iii)*5)>largestPossibleAngle
                                            largestPossibleAngle = 15+selected_constructedLine(iii)*5;
                                        end
                                    else
                                        numOfAnglesSmallerThan_angleTheta1 = numOfAnglesSmallerThan_angleTheta1 + 1;
                                        if (15+selected_constructedLine(iii)*5)<smallestPossibleAngle
                                            smallestPossibleAngle = 15+selected_constructedLine(iii)*5;
                                        end
                                    end
                                    
                                    if angleTheta1>=90
                                        if (15+selected_constructedLine(iii)*5)<=90
                                            forClustering(iii)=15+selected_constructedLine(iii)*5;
                                        end
                                    else
                                        if (15+selected_constructedLine(iii)*5)>=90
                                            forClustering(iii)=15+selected_constructedLine(iii)*5;
                                        end
                                    end
                                end
                            end
                        
                            forClustering = forClustering((forClustering>0));
                            if angleTheta1>=90
                                [forClustering_sort, forClustering_order] = sort(forClustering,'descend');
                                forClustering_sort;
                            else
                                [forClustering_sort, forClustering_order] = sort(forClustering,'ascend');
                                forClustering_sort;
                            end
                            
                            
                            %%% judge how many clusters and find the
                            %%% largetst cluster
                            finalNum = length(forClustering>0);
                            clusterNum = 1;
                            
                            clusterSeg = [];
                            if angleTheta1>=90
                                for xx=1:finalNum-1
                                    if abs(forClustering_sort(xx)-forClustering_sort(xx+1))>5
                                        clusterNum = clusterNum+1;
                                        clusterSeg = [clusterSeg xx];
                                    end
                                end
                                
                                numOfLinesWithinEachCluster = zeros(clusterNum,1);
                                if length(clusterSeg)>0
                                    numOfLinesWithinEachCluster(1) = clusterSeg(1);
                                    for xx=2:clusterNum-1
                                        numOfLinesWithinEachCluster(xx) = clusterSeg(xx)-clusterSeg(xx-1);
                                    end
                                    numOfLinesWithinEachCluster(clusterNum) = finalNum - clusterSeg(clusterNum-1);
                                    maxClusterIndex = find(numOfLinesWithinEachCluster==(max(numOfLinesWithinEachCluster)));
                                else
                                    clusterSeg = finalNum;
                                    maxClusterIndex = 1;
                                end
                                
                                if length(maxClusterIndex)==1
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                elseif length(maxClusterIndex)==2
                                    maxClusterIndex = maxClusterIndex(2);
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                else
                                    maxClusterIndex = maxClusterIndex(length(maxClusterIndex));
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                end
                                
                            else
                                for xx=1:finalNum-1
                                    if abs(forClustering_sort(xx)-forClustering_sort(xx+1))>5
                                        clusterNum = clusterNum+1;
                                        clusterSeg = [clusterSeg xx];
                                    end
                                end
                                numOfLinesWithinEachCluster = zeros(clusterNum,1);
                                if length(clusterSeg)>0
                                    numOfLinesWithinEachCluster(1) = clusterSeg(1);
                                    for xx=2:clusterNum-1
                                        numOfLinesWithinEachCluster(xx) = clusterSeg(xx)-clusterSeg(xx-1);
                                    end
                                    numOfLinesWithinEachCluster(clusterNum) = finalNum - clusterSeg(clusterNum-1);
                                    maxClusterIndex = find(numOfLinesWithinEachCluster==(max(numOfLinesWithinEachCluster)));
                                else
                                    clusterSeg = finalNum;
                                    maxClusterIndex = 1;
                                end
                                
                                if length(maxClusterIndex)==1
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                elseif length(maxClusterIndex)==2
                                    maxClusterIndex = maxClusterIndex(2);
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                else
                                    maxClusterIndex = maxClusterIndex(length(maxClusterIndex));
                                    if maxClusterIndex==1
                                        largestCluster = forClustering_sort(1:clusterSeg(1));
                                    elseif maxClusterIndex<clusterNum
                                        largestCluster = forClustering_sort(clusterSeg(maxClusterIndex-1)+1:clusterSeg(maxClusterIndex));
                                    else
                                        largestCluster = forClustering_sort(clusterSeg(clusterNum-1)+1:finalNum);
                                    end
                                end
                            end

                            clusterSeg;
                            if clusterNum==1
                                mean_forClustering = mean(forClustering((forClustering>0)));
                            elseif clusterNum==2
                                if angleTheta1>=90
                                    if finalNum-clusterSeg(1)>=clusterSeg(1)
                                        mean_forClustering = mean(forClustering_sort(clusterSeg(1)+1:finalNum));
                                    else
                                        mean_forClustering = mean(forClustering_sort(1:clusterSeg(1)));
                                    end
                                else
                                    if finalNum-clusterSeg(1)>=clusterSeg(1)
                                        mean_forClustering = mean(forClustering_sort(clusterSeg(1)+1:finalNum));
                                    else
                                        mean_forClustering = mean(forClustering_sort(1:clusterSeg(1)));
                                    end
                                end
                            else
                                mean_forClustering = mean(forClustering((forClustering>0)));
                            end
                            
                              
                            if (angleTheta1>=90)&&(mean_forClustering<=90)
                                if angleTheta1>90
                                    if mean_forClustering==90
                                        for x=1:newC
                                            y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                    else
                                        for x=1:newC
                                            y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                        for x=newC:imCopyW
                                            y = round(tan(mean_forClustering*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                    end
                                else 
                                    for x=newC:imCopyW
                                        y = round(tan(mean_forClustering*pi/180)*(x-newC)+newR);
                                        if (y>=newR)&&(y<=imCopyH)
                                            doim1(y:imCopyH,x,1) = 255;
                                            roadBinaryImage(y:imCopyH,x) = 255;
                                        end
                                    end
                                end
                                
                            elseif (mean_forClustering>=90)&&(angleTheta1<=90)
                                if angleTheta1<90
                                    if mean_forClustering==90
                                        for x=newC:imCopyW
                                            y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                    else
                                        for x=newC:imCopyW
                                            y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                        for x=1:newC
                                            y = round(tan(mean_forClustering*pi/180)*(x-newC)+newR);
                                            if (y>=newR)&&(y<=imCopyH)
                                                doim1(y:imCopyH,x,1) = 255;
                                                roadBinaryImage(y:imCopyH,x) = 255;
                                            end
                                        end
                                    end
                                else 
                                    for x=1:newC
                                        y = round(tan(mean_forClustering*pi/180)*(x-newC)+newR);
                                        if (y>=newR)&&(y<=imCopyH)
                                            doim1(y:imCopyH,x,1) = 255;
                                            roadBinaryImage(y:imCopyH,x) = 255;
                                        end
                                    end
                                end
                            end
                            
                            for x=1:imCopyW
                                    y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                    if (y>=newR+1)&&(y<=imCopyH-2)
                                        doim1(y-2:y+2,x,1) = 0;
                                        doim1(y-2:y+2,x,2) = 0;
                                        doim1(y-2:y+2,x,3) = 255;
                                    end
                                    
                                    y = round(tan(mean_forClustering*pi/180)*(x-newC)+newR);
                                    if (y>=newR+1)&&(y<=imCopyH-2)
                                        doim1(y-2:y+2,x,1) = 0;
                                        doim1(y-2:y+2,x,2) = 0;
                                        doim1(y-2:y+2,x,3) = 255;
                                    end
                            end
                                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            roadBinaryImage = zeros(imCopyH, imCopyW);
                            if numOfAnglesLargerThan_angleTheta1>=numOfAnglesSmallerThan_angleTheta1+3
                                for x=1:imCopyW
                                    y = round(tan(largestPossibleAngle*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                                for x=1:imCopyW
                                    y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                            elseif numOfAnglesSmallerThan_angleTheta1>=numOfAnglesLargerThan_angleTheta1+3
                                for x=1:imCopyW
                                    y = round(tan(smallestPossibleAngle*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                                for x=1:imCopyW
                                    y = round(tan(angleTheta1*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                            else
                                for x=1:imCopyW
                                    y = round(tan(largestPossibleAngle*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                                for x=1:imCopyW
                                    y = round(tan(smallestPossibleAngle*pi/180)*(x-newC)+newR);
                                    if (y>=newR)&&(y<=imCopyH-1)
                                        doim(y-1:y+1,x,1) = 0;
                                        doim(y-1:y+1,x,2) = 0;
                                        doim(y-1:y+1,x,3) = 255;
                                        doim(y:imCopyH,x,1) = 255;
                                        roadBinaryImage(y:imCopyH,x) = 255;
                                    end
                                end
                            end
                        end
                        
                        imwrite(uint8(roadBinaryImage), [outputPath,int2str(numOfValidFiles),'roadBinary.jpg'], 'jpg');
                        imwrite(uint8(doim), [outputPath,int2str(numOfValidFiles),'vpMap_Road.jpg'], 'jpg');
                        % imwrite(uint8(doim1), [outputPath,int2str(numOfValidFiles),'vpMap_Road1.jpg'], 'jpg');
            
            
            end