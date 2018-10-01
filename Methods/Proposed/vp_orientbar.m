function orientbar = vp_orientbar(grayImg, orientation)

[imgH, imgW] = size(grayImg);
dispaly = orientation;

doim = zeros(imgH, imgW, 3);
doim(:,:,1) = grayImg;
doim(:,:,2) = grayImg;
doim(:,:,3) = grayImg;
doim_resized = imresize(doim,2,'bilinear');
dispaly_resized = imresize(dispaly,2,'bilinear'); 

barLen = 4;
barWid = 1;

for i = 10 : 8 : imgH * 2 - 10
    for j = 10 : 8 : imgW * 2 - 10
        ori = dispaly_resized(i,j); 
        if (ori==90)
            yy = i;
            xx = j;
            doim_resized(yy:yy+barLen,xx-barWid:xx,1) = 255;
            doim_resized(yy:yy+barLen,xx-barWid:xx,2) = 0;
            doim_resized(yy:yy+barLen,xx-barWid:xx,3) = 0;
        elseif (ori==180)||(ori==0)
            doim_resized(i-barWid:i,j:j+barLen,1) = 255;
            doim_resized(i-barWid:i,j:j+barLen,2) = 0;
            doim_resized(i-barWid:i,j:j+barLen,3) = 0;
        else
            if (ori<=45)||(ori>=135)
                kk = tan(ori*pi/180);
                for xx=j:j+barLen
                    yy = round(kk*(xx-j) + i);
                    if (yy>=i-barLen)&&(yy<=i+barLen) 
                        doim_resized(yy,xx-barWid:xx,1) = 255;
                        doim_resized(yy,xx-barWid:xx,2) = 0;
                        doim_resized(yy,xx-barWid:xx,3) = 0;
                    end
                end
            elseif (ori>45)&&(ori<135)
                kk = tan(ori*pi/180);
                for yy=i:i+barLen
                    xx = round((yy-i)/kk + j);
                    if (xx>=j-barLen)&&(xx<=j+barLen) 
                        doim_resized(yy-barWid:yy,xx,1) = 255;
                        doim_resized(yy-barWid:yy,xx,2) = 0;
                        doim_resized(yy-barWid:yy,xx,3) = 0;
                    end
                end
            end
        end       
    end
end

orientbar = uint8(doim_resized);