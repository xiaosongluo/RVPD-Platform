function [convolution, half] = vp_convolution(grayImg, norient)

[imgH, imgW] = size(grayImg);

% Initialization
rang = 180;
interval = rang / norient;
c = 2.2;
w0 = 2.1;


% Set KerSize
lamda = 2^(floor(log2(imgW) - 5));
K = floor(10 * lamda / pi);
if mod(K, 2)==0
    K = K + 1;
end
half = floor(K / 2);

times = 5;
Tconvolution = zeros(imgH, imgW, norient, times);

for it = 0 : times -1
    % Set radial frequency
    w = w0 * 2^it;
    
    odd = zeros(K, K, norient);
    even = zeros(K, K, norient);
    
    % Kernels Generation
    for theta = 90 + 1 : interval : rang - interval + 1 + 90
        tmpTheta = (theta - 1) * pi / 180;
        for y= -half : half
            ySinTheta = y * sin(tmpTheta);
            yCosTheta = y * cos(tmpTheta);
            for x=-half : half
                xCosTheta = x * cos(tmpTheta);
                xSinTheta = x * sin(tmpTheta);
                a = xCosTheta + ySinTheta;
                b = -xSinTheta + yCosTheta;
                odd(y+half+1,x+half+1,(theta - 1-90)/interval+1) = imag(w/(sqrt(2*pi)*c)*exp(-w*w*(4*a*a+b*b)/(8*c*c))*(exp(1i*a*w)-exp(-c*c/2)));
                even(y+half+1,x+half+1,(theta - 1-90)/interval+1) = real(w/(sqrt(2*pi)*c)*exp(-w*w*(4*a*a+b*b)/(8*c*c))*(exp(1i*a*w)-exp(-c*c/2)));
            end
        end
    end
    
    % Normalize Kernels
    nodd = zeros(K, K, norient);
    neven = zeros(K, K, norient);
    for i=1:norient
        tmpKernel = odd(:,:,i)-mean(mean(odd(:,:,i)));
        tmpKernel = tmpKernel/(norm(tmpKernel));
        nodd(:,:,i) = tmpKernel;
        
        tmpKernel = even(:,:,i)-mean(mean(even(:,:,i)));
        tmpKernel = tmpKernel/(norm(tmpKernel));
        neven(:,:,i) = tmpKernel;
    end
    
    % Image Convolution with Gabor Filter Banks
    ImgsOdd = zeros(imgH, imgW, norient);
    ImgsEven = zeros(imgH, imgW, norient);
    for theta=1:norient
        ImgsOdd(:, :, theta) = conv2(double(grayImg), nodd(:, :, theta), 'same');
        ImgsEven(:, :, theta) = conv2(double(grayImg), neven(:, :, theta), 'same');
        Tconvolution(:, :, theta, it + 1) = ImgsOdd(:, :, theta).* ImgsOdd(:, :, theta) + ImgsEven(:, :, theta).* ImgsEven(:, :, theta);
    end
end

Tconvolution(:, :, :, :) = Tconvolution(:, :, :, :) / (K * K);

convolution = zeros(imgH, imgW, norient);
convolution(:, :, :) = Tconvolution(:, :, :, 1) + Tconvolution(:, :, :, 2) +  Tconvolution(:, :, :, 3) +  Tconvolution(:, :, :, 4) +  Tconvolution(:, :, :, 5);
convolution = convolution / 5;