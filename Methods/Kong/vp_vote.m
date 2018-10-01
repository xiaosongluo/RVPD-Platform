function [votingMap, vpX, vpY] = vp_vote(orientation, voter, halfKerSize)

[imgH, imgW] = size(voter);

[rowInd, colInd] = find(double(voter) == 1);
numOfEdgePixels = sum(sum(double(voter)));

votingMap = zeros(imgH, imgW);
interval = 1;
uppPercent = 0.9;
borderWidth = halfKerSize;

diagonal = sqrt(imgH * imgH + imgW * imgW);

radius = round(diagonal*0.35);

for i = 1 : interval : round(imgH * uppPercent)
    for j = borderWidth + 1 : interval : imgW - borderWidth
        for ind = 1 : numOfEdgePixels
            tmpI = rowInd(ind) - i;
            tmpJ = colInd(ind) - j;
            tempDist = sqrt(tmpJ * tmpJ + tmpI * tmpI);
            if (rowInd(ind) > i) && (tempDist <= radius)
                alpha = acos(tmpJ / tempDist) * 180 / pi;
                theta = orientation(rowInd(ind), colInd(ind));
                angleDiffer = abs(alpha - theta);
                dist = tempDist / diagonal;
                if angleDiffer <= (5 /(1 + 2 * dist))
                    votingMap(i,j) = votingMap(i, j) + 1/(1 + (dist * angleDiffer)^2);
                end
            end
        end
    end
end

max_votingMap = max(max(votingMap(1:interval:round(imgH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
min_votingMap = min(min(votingMap(1:interval:round(imgH*uppPercent),borderWidth+1:interval:imgW-borderWidth)));
votingMap1 = (votingMap(1:interval:round(imgH*uppPercent),borderWidth+1:interval:imgW-borderWidth)-min_votingMap)*255/(max_votingMap-min_votingMap);
votingMap(1:interval:round(imgH*uppPercent),borderWidth+1:interval:imgW-borderWidth) = votingMap1;

votingMap(round(uppPercent*imgH):imgH,:) =0;

votingMap = double(uint8(votingMap));

max_votingMap = max(max(votingMap));
[r,c] = find(votingMap == max_votingMap);
if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
end

vpX=c;
vpY=r;
