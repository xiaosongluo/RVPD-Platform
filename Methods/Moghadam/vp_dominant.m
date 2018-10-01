function orientation = vp_dominant(convolution, halfKerSize)

[imgH, imgW, norient] = size(convolution);

orientation = zeros(imgH,imgW);

for i = 1 + halfKerSize : imgH - halfKerSize
    for j = 1 + halfKerSize : imgW - halfKerSize
        %Calculate orientation
		complexResponseVector = convolution(i,j,:);
		[a,b] = sort(complexResponseVector,'descend');
		cx = a(1) - a(4); %* cos(abs((b(1)-b(4)))*45*pi/180);
		cy = a(2) - a(3); %* cos(abs((b(2)-b(3)))*45*pi/180);
        tx = (b(1) - 1) * 45 * pi / 180;
        ty = (b(2) - 1) * 45 * pi / 180;
        vx = cos(tx) * cx + cos(ty) * cy;
        vy = sin(tx) * cx + sin(ty) * cy;
        orientation(i,j) = atan(vy / vx) * 180 / pi;
        if orientation(i,j) < 0
            orientation(i,j) = orientation(i,j) + 180;
        end
    end
end