function [descriptor] = circularLFR(radius)

maskLFR = zeros(radius+1);

radius2 = round(radius/2);
deneme1 = radius2*radius2;
for x=1:radius+1
    for y=1:radius+1
        if(deneme1 >= ((x-(radius2+1))*(x-(radius2+1)) + (y-(radius2+1))*(y-(radius2+1))))
            maskLFR(x,y) = 1;
        end
    end
end
descriptor = maskLFR;