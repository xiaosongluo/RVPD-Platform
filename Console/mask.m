function  img=mask(px, py, origin, r, g, b)

[imgH, imgW] = size(origin);

img = origin;

coord = [px, py];

if (coord(1)<= 7) 
    coord(1) = 8; 
end

if (coord(1)>=imgH-8)
    coord(1)=imgH-8;
end

if (coord(2)<=7)
    coord(2) = 8;
end

if (coord(2)>=imgW-8)
    coord(2) = imgW-8;
end

y = coord(1);
x = coord(2);

img(y-3:y+3,x-7:x+7,1) = r;
img(y-7:y+7,x-3:x+3,1) = r;
img(y-3:y+3,x-7:x+7,2) = g;
img(y-7:y+7,x-3:x+3,2) = g;
img(y-3:y+3,x-7:x+7,3) = b;
img(y-7:y+7,x-3:x+3,3) = b;