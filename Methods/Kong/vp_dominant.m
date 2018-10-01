function [orientation, confidence] = vp_dominant(convolution, halfKerSize)

[imgH, imgW, norient] = size(convolution);

confidence = zeros(imgH,imgW);
orientation = zeros(imgH,imgW);

for i = 1 + halfKerSize : imgH - halfKerSize
    for j = 1 + halfKerSize : imgW - halfKerSize
        %Calculate confidence
        convolutionVector = convolution(i,j,:);
        [a,b] = sort(convolutionVector,'descend');
        if a(1) > 1
            confidence(i,j) = 1 - mean(a(5:15)) / a(1);
        end
        %Calculate orientation
        maxma = max(convolution(i,j,:));
        indx = find(convolution(i,j,:)==maxma);
        orientation(i,j) = mean(indx);
    end
end

angleRange = 180;
angleInterval = angleRange / norient;
orientation = (orientation - 1) * angleInterval;

% Normalize Confidence
maxV = max(max(confidence(1+halfKerSize:imgH-halfKerSize,1+halfKerSize:imgW-halfKerSize)));
minV = min(min(confidence(1+halfKerSize:imgH-halfKerSize,1+halfKerSize:imgW-halfKerSize)));

aaTmp = confidence(1+halfKerSize:imgH-halfKerSize,1+halfKerSize:imgW-halfKerSize);
aaTmp = (aaTmp-minV)*255/(maxV-minV);
confidence(1+halfKerSize:imgH-halfKerSize,1+halfKerSize:imgW-halfKerSize) = aaTmp;