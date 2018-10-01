function voter = vp_voter(edgeImg, orientation, overlap)

[imgH, imgW] =size(edgeImg);

outlierBinary = ones(imgH,imgW);
barlen = 30; %% odd number
halfbarlen = floor(barlen/2);
for i = 1 + 5 : imgW - 5
    for j = halfbarlen + 1 : imgH - halfbarlen
        tmpSum = sum(edgeImg(j - halfbarlen : j + halfbarlen, i));
        if tmpSum >= 20
            outlierBinary(j - halfbarlen : j + halfbarlen, i - 5 : i + 5) = 0;
        end
    end
end

nonHorizontal_edges = ones(imgH,imgW);
nonHorizontal_edges((orientation==180))=0;

voter = ones(imgH,imgW).*outlierBinary.*nonHorizontal_edges.*overlap;