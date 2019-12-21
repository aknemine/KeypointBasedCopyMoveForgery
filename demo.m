close all;
clear all;clc

forgedImage = imread('image2_org.png');
forgedGrayImage = double(rgb2gray(forgedImage));		

%//=======================================================================
%// Superpixel segmentation
%//=======================================================================
%// numberOfSegments is the target number of superpixels.
numberOfSegments = 50;

[labels] = mex_ers(forgedGrayImage,numberOfSegments);

%//Each segment is given a uniqe label.

[grayImageHeight,grayImageWidth] = size(forgedGrayImage);

[numberOfLabels] = zeros(1,numberOfSegments);
%// numberOfLabels keeps how many of each label.
[minimum] = 255*ones(1,numberOfSegments);
%***********************************************
for i = 1:(grayImageHeight)
    for j = 1:(grayImageWidth)
        for k = 0:(numberOfSegments-1)
            temp = labels(i,j);
            if(temp == k)
                numberOfLabels(k+1) = numberOfLabels(k+1) + 1;
                if(minimum(k+1) >= forgedGrayImage(i,j))
                    minimum(k+1) = forgedGrayImage(i,j);
                end
            end
        end
    end
end
clear i;
clear j;
clear k;

%//=======================================================================
%//                     Superpixel Classification                      //%
%//=======================================================================
[Ek] = zeros(1,numberOfSegments);

for t = 1:numberOfSegments
    Ek(t) = -1*((numberOfLabels(t)/(grayImageHeight*grayImageWidth))*log2(numberOfLabels(t)/(grayImageHeight*grayImageWidth)));
end
clear t;

E = max(Ek) - min(Ek);

EA = min(Ek) + E/2;
EB = min(Ek) + 3*E/4;

%// Color process
coloring = labels;

Bk = zeros(1,numberOfSegments);
for k = 1:numberOfSegments
    if(Ek(k) < EA)
        Bk(k) = 0;
    elseif(EA < Ek(k) && Ek(k) < EB)
        Bk(k) = 1;
    elseif(Ek(k)> EB)
        Bk(k) = 2;
    end
    
    for i = 1:(grayImageHeight)
        for j = 1:(grayImageWidth)
            if(k == labels(i,j))
                coloring(i,j) = Bk(k);
            end
        end
    end
end
clear i;
clear j; 
clear k;

[Tk] = zeros(1,numberOfSegments);

for t = 1:numberOfSegments
    if(Bk(t) == 0)
        Tk(t) = (0.0002 * 0.01*0.01 * minimum(t));
    elseif (Bk(t) == 1)
        Tk(t) = (0.0002 * 0.01*0.01 * minimum(t));
    elseif (Bk(t) == 2)
        Tk(t) = (0.0002 * minimum(t));
    end
end
clear t;

blueIndex = find(coloring==0);
greenIndex = find(coloring==1);
redIndex = find(coloring==2);

coloredImage = zeros(grayImageHeight, grayImageWidth, 3, 'uint8');

timg1 = coloredImage(:,:,1);
timg1(redIndex) = 255;

timg2 = coloredImage(:,:,2);
timg2(greenIndex) = 255;

timg3 = coloredImage(:,:,3);
timg3(blueIndex) = 255;

coloredImage(:,:,1) = timg1;
coloredImage(:,:,2) = timg2;
coloredImage(:,:,3) = timg3;

[borderOfSegments] = seg2bmap(labels,grayImageWidth,grayImageHeight);
bmapOnImg = forgedImage;
idx = find(borderOfSegments>0);
timg = forgedImage(:,:,1);
timg(idx) = 255;
bmapOnImg(:,:,1) = timg;
bmapOnImg(:,:,2) = forgedImage(:,:,2);
bmapOnImg(:,:,3) = forgedImage(:,:,3);
%  
%  imshow(bmapOnImg);
%  hold on;

%//---------------------------------------------------------------------\\%
%                   SURF - Speeded-Up Robust Features                     %
%//---------------------------------------------------------------------\\%

imageForSurf = rgb2gray(mat2gray(forgedImage));
keyPoints = detectSURFFeatures(imageForSurf);

scaleOfPoint = keyPoints.Scale;

surfPointsSize = keyPoints.size;
[diameter] = zeros(surfPointsSize(1,1),1);

strongestKeyPoints = keyPoints.selectStrongest(surfPointsSize(1,1));
[locationOfKeyPoint] = round(strongestKeyPoints.Location);

To = 2;
diameter(:) = To * round(scaleOfPoint(:));
% 
% imshow(forgedImage);
% hold on;
% plot(keyPoints.selectStrongest(250));

%//---------------------------------------------------------------------\\%
%              Log Polar Transform for Keypoint Descriptor                %
%//---------------------------------------------------------------------\\%

pointOfSurf(surfPointsSize(1,1)).descriptors = zeros(1,64);
pointOfSurf(surfPointsSize(1,1)).points = zeros(1,2);

for i=1:surfPointsSize(1,1)
    maskLFROutput = circularLFR(diameter(i));
    circularOfPoint = zeros(diameter(i) + 1);
    counter1=1;
    for j = locationOfKeyPoint(i,2) - round(diameter(i)/2) : locationOfKeyPoint(i,2) + round(diameter(i)/2) 
        counter2=1;
        for k = locationOfKeyPoint(i,1) - round(diameter(i)/2) : locationOfKeyPoint(i,1) + round(diameter(i)/2)
            circularOfPoint(counter1,counter2) = imageForSurf(j,k) * maskLFROutput(counter1,counter2);
            counter2 = counter2 + 1;
        end
       counter1 = counter1 + 1;
    end
    
    rMax = round(length(circularOfPoint)/2);
    tempVektor = logsample(circularOfPoint,1,rMax,rMax,rMax,8,8);
    for t = 1:8
        for y = 1:8
            pointOfSurf(i).descriptors((t-1)*8+y) = tempVektor(t,y);
        end
    end
    pointOfSurf(i).points(1) = locationOfKeyPoint(i,2);
    pointOfSurf(i).points(2) = locationOfKeyPoint(i,1);
end
clear i;
clear j;
clear k;
clear t;
clear y;

emre1(1).descriptors = zeros(surfPointsSize(1,1),64);
emre1(1).points = zeros(surfPointsSize(1,1),2);

for yu=1:surfPointsSize(1,1)
    emre1(1).descriptors(yu,:) = pointOfSurf(yu).descriptors(:);
    emre1(1).points(yu,:) = pointOfSurf(yu).points(:);
end

%//---------------------------------------------------------------------\\%
%                             G2NN Matching                               %
%//---------------------------------------------------------------------\\%

g2nnPointsStruct = g2nn(emre1,emre1,0.1);
matchPoints = filter_small_matches(g2nnPointsStruct,100);

sourcePoints = matchPoints.source;
targetPoints = matchPoints.target;


g2nnPointsStruct1 = g2nn(emre1,emre1,0.8);
matchPoints1 = filter_small_matches(g2nnPointsStruct1,40);

sourcePoints1 = matchPoints1.source;
targetPoints1 = matchPoints1.target;

boyanacaklar = ones(grayImageHeight,grayImageWidth);

deneme = extension(sourcePoints,targetPoints,boyanacaklar,forgedGrayImage);
indisler = find(deneme == -1);
indisler2 = find(deneme == 0);
tempImage = forgedImage;
maske = imread('maske_image1.png');
maske1 = imread('image2_mask.png');
boya = zeros(length(indisler)+length(indisler2),2);
for i = 1 : length(indisler)
    boya(i,1) = round(indisler(i)/300);
    boya(i,2) = mod(indisler(i),300);
end

for i = 1 : length(indisler2)
    boya(i+length(indisler),1) = round(indisler2(i)/300);
    boya(i+length(indisler),2) = mod(indisler2(i),300);
end

for i = 1 : length(boya)
    maske(boya(i,2)+1,boya(i,1)+1,1) = 255;
    maske(boya(i,2)+1,boya(i,1)+1,2) = 255;
    maske(boya(i,2)+1,boya(i,1)+1,3) = 255;
    tempImage(boya(i,2)+1,boya(i,1)+1,1) = 255;
    %tempImage(boya(i,2)+1,boya(i,1)+1,3) = 255;
end
%imwrite(maske,'deneme3.png');
beklenen = imread('deneme3.png');

gcf = figure(1);
subplot(2,4,1);
imshow(bmapOnImg,[]);
title('Bölütlenmiþ Görüntü');
subplot(2,4,2);
imshow(coloredImage,[]);
title('Class Görüntüsü');
subplot(2,4,3);
imshow(forgedImage,[]);
title('Bölütlenmiþ Görüntü');
hold on
plot(keyPoints.selectStrongest(50));
subplot(2,4,4);
imshow(forgedImage);
hold on
for i = 1 : length(targetPoints1)
     A(1,1) = sourcePoints1(i,1);
     A(1,2) = targetPoints1(i,1);
     B2(1,1) = sourcePoints1(i,2);
     B2(1,2) = targetPoints1(i,2);
     plot([B2(1,1), B2(1,2)],[A(1,1) ,A(1,2)],'Color','b','LineWidth',1);   
 end
title('Matching Görüntü');
subplot(2,4,5);
imshow(forgedImage);
hold on
for i = 1 : length(targetPoints)
     A(1,1) = sourcePoints(i,1);
     A(1,2) = targetPoints(i,1);
     B2(1,1) = sourcePoints(i,2);
     B2(1,2) = targetPoints(i,2);
     plot([B2(1,1), B2(1,2)],[A(1,1) ,A(1,2)],'Color','b','LineWidth',1);   
 end
title('Filtreli Matching Görüntü');
subplot(2,4,6);
imshow(tempImage,[]);
title('Renkli Görüntü');
subplot(2,4,7);
imshow(beklenen,[]);
title('Ýyileþtirme Öncesi Maske');
subplot(2,4,8);
imshow(maske,[]);
title('Ýyileþtirme Sonrasý Maske');
% subplot(3,3,9);
% imshow(maske,[]);
% title('Ýyileþtirme Sonrasý Maske');

gerceksonuc = getFmeasure(maske,maske1);

%// Compute the superpixel size histogram.
sizeOfSuperpixel = zeros(numberOfSegments,1);
for i=0:(numberOfSegments-1)
    sizeOfSuperpixel(i+1) = sum( labels(:)==i );
end
clear i;

[his, bins] = hist( sizeOfSuperpixel, 20 );
sizeOfKeyPoints = length(locationOfKeyPoint);
% figure , imshow(tempImage)
% hold on
% %plot(keyPoints.selectStrongest(sizeOfKeyPoints));
% for i = 1 : length(targetPoints)
%      A(1,1) = sourcePoints(i,1);
%      A(1,2) = targetPoints(i,1);
%      B2(1,1) = sourcePoints(i,2);
%      B2(1,2) = targetPoints(i,2);
%      plot([B2(1,1), B2(1,2)],[A(1,1) ,A(1,2)],'Color','b','LineWidth',1);   
%  end

%%
%//=======================================================================
%// Display 
%//=======================================================================
% gcf = figure(1);
% subplot(2,3,1);
% imshow(img,[]);
% title('input image.');
% subplot(2,3,2);
% imshow(bmapOnImg,[]);
% title('superpixel boundary map');
% subplot(2,3,3);
% imshow(out,[]);
% title('randomly-colored superpixels');
% subplot(2,3,4);
% imshow(coloredImage,[]);
% title('surffff image.');
% subplot(2,3,5);
% bar(bins,his,'b');
% title('the distribution of superpixel size');
% ylabel('# of superpixels');
% xlabel('superpixel sizes in pixel');
% scnsize = get(0,'ScreenSize');
% set(gcf,'OuterPosition',scnsize);