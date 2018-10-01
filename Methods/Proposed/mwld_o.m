function orient = mwld_o(gray, norient, r)

[imgH, imgW] = size(gray);

% Initialization 
range = 180;
interval = range / norient;

k = 2 * r + 1;
delta = k / 9;
lamda = round(k * pi / 10);
tmpDelta = -1/(delta * delta * 8);
c = 2 * pi / lamda;
cc = 2.2;

odd = zeros(k, k, norient);
even = zeros(k, k, norient);

% Kernels Generation
for theta = 90 + 1 : interval : range - interval + 1 + 90
    tmpTheta = (theta - 1) * pi / 180;
    for y= -r : r
        ySinTheta = y * sin(tmpTheta);
        yCosTheta = y * cos(tmpTheta);
        for x = -r : r
            xCosTheta = x * cos(tmpTheta);
            xSinTheta = x * sin(tmpTheta);
            a = xCosTheta + ySinTheta;
            b = -xSinTheta + yCosTheta;
            odd(y+r+1, x+r+1, (theta-1-90)/interval+1) = exp(tmpDelta*(4*a*a+b*b))*(sin(c*a)-exp(-cc*cc/2));
            even(y+r+1, x+r+1, (theta-1-90)/interval+1) = exp(tmpDelta*(4*a*a+b*b))*(cos(c*a)-exp(-cc*cc/2));
        end
    end
end

% Normalize Kernels
nodd = zeros(k, k, norient);
neven = zeros(k, k, norient);
for i = 1 : norient
	tmp = odd(:,:,i) - mean(mean(odd(:,:,i)));
	nodd(:,:,i) = tmp/(norm(tmp));
	tmp = even(:,:,i) - mean(mean(even(:,:,i)));
	neven(:,:,i) = tmp/(norm(tmp));
end

% Image convolution with Gabor Filter Banks
fodd = zeros(imgH, imgW, norient);
feven = zeros(imgH, imgW, norient);
convolution = zeros(imgH, imgW, norient);
for i = 1 : norient
    fodd(:,:,i) = conv2(double(gray), nodd(:,:,i), 'same');
    feven(:,:,i) = conv2(double(gray), neven(:,:,i), 'same');
    convolution(:,:,i) = fodd(:,:,i) .* fodd(:,:,i) + feven(:,:,i) .* feven(:,:,i);
    convolution(:,:,i) = convolution(:,:,i) / (k * k);
end

%calculate orientation
orient = zeros(imgH,imgW);
for i = 1 + r : imgH - r
    for j = 1 + r : imgW - r
        maxma = max(convolution(i,j,:));
        indx = find(convolution(i,j,:) == maxma);
        orient(i,j) = mean(indx);
    end
end
orient = (orient - 1) * interval;