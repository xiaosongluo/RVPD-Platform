%%%%%%%% convoluation of an image with gLoG filters of the same
%%%%%%%% orientation 

function [responseMap] = convolve_gLoG(img, smallestSigma, largestSigma, theta)

imgH = size(img,1);
imgW = size(img,2);
responseMap = zeros(imgH,imgW);

sigmaStep = -1;
newKerSize = 3*largestSigma;
hsize2 = newKerSize/2;
alpha = 1;

for sigmaX = largestSigma : sigmaStep: smallestSigma;
    for sigmaY = sigmaX : sigmaStep: smallestSigma;
        if sigmaX~=sigmaY
           [h] = -elipLog([newKerSize+1,newKerSize+1], sigmaX, sigmaY, theta);
           filteredImg = real(ifft2(fft2(img) .* fft2(h,imgH,imgW)));
           tmpMap = filteredImg;
           tmpMap(1:imgH-hsize2, 1:imgW-hsize2) = filteredImg(1+hsize2:imgH, 1+hsize2:imgW);
           tmpMap(1:imgH-hsize2,imgW-hsize2+1:imgW) = filteredImg(hsize2+1:imgH,1:hsize2);
           tmpMap(imgH-hsize2+1:imgH,1:imgW-hsize2) = filteredImg(1:hsize2,hsize2+1:imgW);
           tmpMap(imgH-hsize2+1:imgH, imgW-hsize2+1:imgW) = filteredImg(1:hsize2,1:hsize2);
           
           tmpMap = (1+log(sigmaX^alpha))*(1+log(sigmaY^alpha))*tmpMap;
           responseMap = responseMap+tmpMap;
        end
    end      
end

