function [votingMap, vpX, vpY] = vp_vote(orientation, voter)

[imgH, imgW] = size(voter);

[rowInd, colInd] = find(double(voter) == 1);
numOfEdgePixels = sum(sum(double(voter)));

votingMap = zeros(imgH, imgW);
for ind = 1 : numOfEdgePixels
    theta = orientation(rowInd(ind), colInd(ind));
    cotangent = cot(theta / 180 * pi);
    if (theta >= 45 && theta <= 135)
        for i = 1 : imgH
            j = colInd(ind) - round((rowInd(ind) - i) * cotangent);
            if (1 <= j)&&(j <= imgW)
                votingMap(i, j) = votingMap(i, j) + 1;
            end
        end
    else
        for j = 1 : imgW
            i = rowInd(ind) - round((colInd(ind) - j) / cotangent);
            if (1 <= i)&&(i <= imgH)
                votingMap(i, j) = votingMap(i, j) + 1;
            end
        end
    end
end

max_votingMap = max(max(votingMap));
min_votingMap = min(min(votingMap));
votingMap = (votingMap - min_votingMap) * 255 / (max_votingMap - min_votingMap);

votingMap = double(uint8(votingMap));

windows = 5;

sum_votingMap = zeros(imgH,imgW);
for i=1+windows:imgH-windows
    for j=1+windows:imgW-windows
        sum_votingMap(i,j)=sum(sum(votingMap(i-windows:i+windows,j-windows:j+windows)));
    end
end

max_votingMap = max(max(sum_votingMap));
[r,c] = find(sum_votingMap==max_votingMap);
if length(r)>1
    r = round(mean(r));
    c = round(mean(c));
end

vpX=c;
vpY=r;