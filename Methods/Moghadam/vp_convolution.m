function [convolution, halfKerSize] = vp_convolution(grayImg)

[imgH, imgW] = size(grayImg);

% Initialization 
angleRange = 180;
angleInterval = 45;

lamda = 4 * sqrt(2);
kerSize = floor(10 * lamda / pi);

% Set KerSize
if mod(kerSize, 2)==0
    kerSize = kerSize + 1;    
end
halfKerSize = floor(kerSize / 2);

oddKernel = zeros(kerSize, kerSize, 4);
evenKernel = zeros(kerSize, kerSize, 4);

w0 = 2 * pi / lamda;
K = pi / 2;
tmpDelta = -w0 * w0/(K * K * 8);

% Kernels Generation
for theta = 90 + 1 : angleInterval : angleRange - angleInterval + 90 + 1
    tmpTheta = theta * pi/180;
    for y= -halfKerSize:halfKerSize
        ySinTheta = y*sin(tmpTheta);
        yCosTheta = y*cos(tmpTheta);
        for x=-halfKerSize:halfKerSize
            xCosTheta = x*cos(tmpTheta);
            xSinTheta = x*sin(tmpTheta);
            a = xCosTheta+ySinTheta;
            b = -xSinTheta+yCosTheta;
            oddKernel(y+halfKerSize+1,x+halfKerSize+1,(theta-90-1)/angleInterval+1) = imag(w0/(sqrt(2*pi)*K)*exp(tmpDelta*(4*a*a+b*b))*(exp(i*w0*a)-exp(-K*K/2)));
            evenKernel(y+halfKerSize+1,x+halfKerSize+1,(theta-90-1)/angleInterval+1) = real(w0/(sqrt(2*pi)*K)*exp(tmpDelta*(4*a*a+b*b))*(exp(i*w0*a)-exp(-K*K/2)));
        end
    end
end

% Image Convolution with Gabor Filter Banks
filteredImgsOdd = zeros(imgH,imgW,4*1);
filteredImgsEven = zeros(imgH,imgW,4*1);
convolution = zeros(imgH,imgW,4*1);
for x=1:4
    filteredImgsOdd(:,:,x) = conv2(double(grayImg), oddKernel(:,:,x), 'same');
    filteredImgsEven(:,:,x) = conv2(double(grayImg), evenKernel(:,:,x), 'same');
    convolution(:,:,x) = sqrt(filteredImgsOdd(:,:,x).*filteredImgsOdd(:,:,x) + filteredImgsEven(:,:,x).*filteredImgsEven(:,:,x));
end