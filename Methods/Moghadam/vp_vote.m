function [votingMap, vpX, vpY] = vp_vote(orientation, halfKerSize)

[imgH, imgW] = size(orientation);

votingMap = zeros(imgH, imgW);
interval = 1;
borderWidth = halfKerSize;
variance = 0.25;

[rowInd, colInd] = find(double(orientation) ~= 0);
numOfEdgePixels = sum(sum(double(orientation) ~= 0));

for ind = 1 : numOfEdgePixels
    theta = orientation(rowInd(ind), colInd(ind));
    sine = sin(theta / 180 * pi);
    tg = tan(theta / 180 * pi);
    maxD = 0;
    if (theta >= 45 && theta <= 135)
        for i = 1 : rowInd(ind)-1
            j = colInd(ind) - round((rowInd(ind) - i) / tg);
            if (1 <= j)&&(j <= imgW)
                tmpI = rowInd(ind) - i;
                tmpJ = colInd(ind) - j;
                tempDist = sqrt(tmpJ * tmpJ + tmpI * tmpI);
                if maxD == 0
                    maxD = tempDist;
                end
                dist = tempDist / maxD;
                votingMap(i,j) = votingMap(i, j) + sine * exp(- dist^2 / (2 * variance));
            end
        end
    else
        for j = 1 : imgW
            i = rowInd(ind) - round((colInd(ind) - j) * tg);
            if (1 <= i)&&(i <= rowInd(ind)-1)
                tmpI = rowInd(ind) - i;
                tmpJ = colInd(ind) - j;
                tempDist = sqrt(tmpJ * tmpJ + tmpI * tmpI);
                if maxD == 0
                    maxD = tempDist;
                end
                dist = tempDist / maxD;
                votingMap(i,j) = votingMap(i, j) + sine * exp(- dist^2 / (2 * variance));
            end
        end
    end
end

max_votingMap = max(max(votingMap(1 : interval : imgH, borderWidth+1 : interval : imgW - borderWidth)));
min_votingMap = min(min(votingMap(1 : interval : imgH, borderWidth+1 : interval : imgW - borderWidth)));
votingMap1 = (votingMap(1 : interval : imgH, borderWidth+1 : interval : imgW - borderWidth) - min_votingMap) * 255 / (max_votingMap - min_votingMap);
votingMap(1 : interval : imgH, borderWidth+1 : interval : imgW - borderWidth) = votingMap1;

votingMap = double(uint8(votingMap));

max_votingMap = max(max(votingMap));
[r,c] = find(votingMap == max_votingMap);
if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
end

vpX=c;
vpY=r;