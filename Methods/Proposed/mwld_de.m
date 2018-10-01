function diffex = mwld_de(gray, r)

% Determine the dimensions of the input image.
image = double(gray);
[imgH, imgW] = size(gray);
kerSize = r * 2 + 1;

% filter for differential excitation
kernel = ones(kerSize, kerSize);
kernel(r + 1, r + 1) = 0;
kernel = kernel / (kerSize * kerSize - 1);


% compute differential excitation
convolve = conv2(image, kernel, 'same'); %convolve with kernel
scaling = atan((image - convolve) ./ convolve); %perform the tangent scaling
diffex = 2 * scaling / pi;
diffex(image==0) = 0;
diffex(diffex<0) = 0;
diffex = sqrt(diffex);

diffex(1:r, :) = 0;
diffex(imgH-r:imgH, :) = 0;
diffex(:, 1:r) = 0;
diffex(:, imgW-r:imgW) = 0;