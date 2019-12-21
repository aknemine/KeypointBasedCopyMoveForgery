close all;
clear all;

disp('Entropy Rate Superpixel Segmentation Demo');

img = imread('forged1.png');
deneme3 = img;
%// Our implementation can take both color and grey scale images.
grey_img = double(rgb2gray(img));		

%%
%//=======================================================================
%// Superpixel segmentation
%//=======================================================================
%// nC is the target number of superpixels.
nC = 100;		
%// Call the mex function for superpixel segmentation\
%// !!! Note that the output label starts from 0 to nC-1.
t = cputime;

lambda_prime = 0.5;
sigma = 5.0; 	
conn8 = 1; 

%// flag for using 8 connected grid graph (default setting).
%[labels] = mex_ers(double(img),nC);		
%[labels] = mex_ers(double(img),nC,lambda_prime,sigma);
%[labels] = mex_ers(double(img),nC,lambda_prime,sigma,conn8);

% grey scale iamge
[labels] = mex_ers(grey_img,nC);
%[labels] = mex_ers(grey_img,nC,lambda_prime,sigma);
%[labels] = mex_ers(grey_img,nC,lambda_prime,sigma,conn8);

fprintf(1,'Use %f sec. \n',cputime-t);
fprintf(1,'\t to divide the image into %d superpixels.\n',nC);

%// You can also specify your preference parameters. The parameter values
%// (lambda_prime = 0.5, sigma = 5.0) are chosen based on the experiment
%// results in the Berkeley segmentation dataset.
%// lambda_prime = 0.5; sigma = 5.0;
%// [labels] = mex_ers(grey_img,nC,lambda_prime,sigma);
%// You can also use 4 connected-grid graph. The algorithm uses 8-connected 
%// graph as default setting. By setting conn8 = 0 and running
%// [labels] = mex_ers(grey_img,nC,lambda_prime,sigma,conn8),
%// the algorithm perform segmentation uses 4-connected graph. Note that 
%// 4 connected graph is faster.




%%
%//=======================================================================
%// Output
%//=======================================================================
[height,width] = size(grey_img);

[T] = zeros(1,nC);
[minimum] = 255*ones(1,nC);
%***********************************************
for i = 1:(height)
    for j = 1:(width)
        for k = 0:(nC-1)
            a = labels(i,j);
            if(a == k)
                T(k+1) = T(k+1) + 1;
                if(minimum(k+1) >= grey_img(i,j))
                    minimum(k+1) = grey_img(i,j);
                end
            end
        end
    end
end

%toplam = sum(T);
[E] = zeros(1,nC);

for t = 1:nC
    E(t) = -1*((T(t)/(height*width))*log2(T(t)/(height*width)));
end

ValuE = max(E) - min(E);

EA = min(E) + ValuE/2;
EB = min(E) + 3*ValuE/4;

boyama = labels;

B = zeros(1,nC);
for k = 1:nC
    if(E(k) < EA)
        B(k) = 0;
    elseif(EA < E(k) && E(k) < EB)
        B(k) = 1;
    elseif(E(k)> EB)
        B(k) = 2;
    end
    
    for i = 1:(height)
        for j = 1:(width)
            if(k == labels(i,j))
                boyama(i,j) = B(k);
            end
        end
    end
end

[Tk] = zeros(1,nC);

for t = 1:nC
    if(B(t) == 0)
        Tk(t) = (0.0002 * 0.01*0.01 * minimum(t));
    elseif (B(t) == 1)
        Tk(t) = (0.0002 * 0.01*0.01 * minimum(t));
    elseif (B(t) == 2)
        Tk(t) = (0.0002 * minimum(t));
    end
end

mavi_idx = find(boyama==0);
yesil_idx = find(boyama==1);
kirmizi_idx = find(boyama==2);

surfluimage = zeros(height, width, 3, 'uint8');

timg1 = surfluimage(:,:,1);
timg1(kirmizi_idx) = 255;

timg2 = surfluimage(:,:,2);
timg2(yesil_idx) = 255;

timg3 = surfluimage(:,:,3);
timg3(mavi_idx) = 255;

surfluimage(:,:,1) = timg1;
surfluimage(:,:,2) = timg2;
surfluimage(:,:,3) = timg3;


%// Compute the boundary map and superimpose it on the input image in the
%// green channel.
%// The seg2bmap function is directly duplicated from the Berkeley
%// Segmentation dataset which can be accessed via
%// http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/

[bmap] = seg2bmap(labels,width,height);
bmapOnImg = img;
idx = find(bmap>0);
timg = img(:,:,1);
timg(idx) = 255;
bmapOnImg(:,:,1) = timg;
bmapOnImg(:,:,2) = img(:,:,2);
bmapOnImg(:,:,3) = img(:,:,3);

deneme4 = rgb2gray(mat2gray(img));
points = detectSURFFeatures(deneme4);
%%[features,validPoints] = extractFeatures(deneme4,points);
noktalar = points.Scale;

surfPointsSize = points.size;
[radius] = zeros(surfPointsSize(1,1),1);

strongest = points.selectStrongest(surfPointsSize(1,1));
[locationPoint] = round(strongest.Location);

radius(:) = 2 * round(noktalar(:));

emre(surfPointsSize(1,1)).descriptors = zeros(1,64);
emre(surfPointsSize(1,1)).points = zeros(1,2);

for index1=1:surfPointsSize(1,1)
    maskLFROutput = circularLFR(radius(index1));
    circularOfPoint = zeros(radius(index1) + 1);
    counter1=1;
    for index2 = locationPoint(index1,2) - round(radius(index1)/2) : locationPoint(index1,2) + round(radius(index1)/2) 
        counter2=1;
        for index3 = locationPoint(index1,1) - round(radius(index1)/2) : locationPoint(index1,1) + round(radius(index1)/2)
            circularOfPoint(counter1,counter2) = deneme4(index2,index3) * maskLFROutput(counter1,counter2);
            counter2 = counter2 + 1;
        end
       counter1 = counter1 + 1;
    end
    [ab,ba] = size(circularOfPoint);
    rMax = round(ab/2);
    geciciVektor = logsample(circularOfPoint,1,rMax,rMax,rMax,8,8);
    for indis = 1:8
        for indis2 = 1:8
            emre(index1).descriptors((indis-1)*8+indis2) = geciciVektor(indis,indis2);
        end
    end
    emre(index1).points(1) = locationPoint(index1,2);
    emre(index1).points(2) = locationPoint(index1,1);
end

emre1(1).descriptors = zeros(surfPointsSize(1,1),64);
emre1(1).points = zeros(surfPointsSize(1,1),2);

for yu=1:surfPointsSize(1,1)
    emre1(1).descriptors(yu,:) = emre(yu).descriptors(:);
    emre1(1).points(yu,:) = emre(yu).points(:);
end
bizt = g2nn(emre1,emre1,0.10);
matches = filter_small_matches(bizt,40);
% 
% 
% if size(matches.source, 1) > 1
%         matches = unique_matches(matches);
%         matches = partition_matches(matches);
%         
%         %% RANSAC ESTIMATION
%         disp('RANSAC..');
%         RANSAC_ITERS = 15000;
%         start_time = cputime;
%         RANSAC_ERR_THRESH = 5;
%         r_matches = ransac_le_means(matches, RANSAC_ITERS, RANSAC_ERR_THRESH);
%     else
%         r_matches = matches;
% end

noktamiz = bizt.source;
hedefimiz = bizt.target;
i = size(noktamiz);
img45 = img;
for a = 1 : i(1,1)
    img45(noktamiz(a,1), noktamiz(a,2),1) = 255;
    img45(noktamiz(a,1), noktamiz(a,2),2) = 255;
    img45(noktamiz(a,1), noktamiz(a,2),3) = 255;

    img45(noktamiz(a,1)+1, noktamiz(a,2),1) = 255;
    img45(noktamiz(a,1)+1, noktamiz(a,2),2) = 255;
    img45(noktamiz(a,1)+1, noktamiz(a,2),3) = 255;

    img45(noktamiz(a,1)-1, noktamiz(a,2),1) = 255;
    img45(noktamiz(a,1)-1, noktamiz(a,2),2) = 255;
    img45(noktamiz(a,1)-1, noktamiz(a,2),3) = 255;

    img45(noktamiz(a,1), noktamiz(a,2)+1,1) = 255;
    img45(noktamiz(a,1), noktamiz(a,2)+1,2) = 255;
    img45(noktamiz(a,1), noktamiz(a,2)+1,3) = 255;

    img45(noktamiz(a,1), noktamiz(a,2)-1,1) = 255;
    img45(noktamiz(a,1), noktamiz(a,2)-1,2) = 255;
    img45(noktamiz(a,1), noktamiz(a,2)-1,3) = 255;

    img45(hedefimiz(a,1),hedefimiz(a,2),1) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2),2) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2),3) = 255;

    img45(hedefimiz(a,1)+1,hedefimiz(a,2),1) = 255;
    img45(hedefimiz(a,1)+1,hedefimiz(a,2),2) = 255;
    img45(hedefimiz(a,1)+1,hedefimiz(a,2),3) = 255;

    img45(hedefimiz(a,1)-1,hedefimiz(a,2),1) = 255;
    img45(hedefimiz(a,1)-1,hedefimiz(a,2),2) = 255;
    img45(hedefimiz(a,1)-1,hedefimiz(a,2),3) = 255;

    img45(hedefimiz(a,1),hedefimiz(a,2)+1,1) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2)+1,2) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2)+1,3) = 255;

    img45(hedefimiz(a,1),hedefimiz(a,2)-1,1) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2)-1,2) = 255;
    img45(hedefimiz(a,1),hedefimiz(a,2)-1,3) = 255;

end
%// Randomly color the superpixels
[out] = random_color( double(img) ,labels,nC);
% 
% imshow(img45);
% hold on;
%plot(selectStrongest(points, surfPointsSize(1,1)));

%// Compute the superpixel size histogram.
siz = zeros(nC,1);
for i=0:(nC-1)
    siz(i+1) = sum( labels(:)==i );
end
[his, bins] = hist( siz, 20 );
  
figure , imshow(img45)
hold on
plot(points.selectStrongest(200));
for i = 1 : length(hedefimiz)
 A(1,1) = noktamiz(i,1);
  A(1,2) = hedefimiz(i,1);
  B2(1,1) = noktamiz(i,2);
  B2(1,2) = hedefimiz(i,2);
    plot([B2(1,1), B2(1,2)],[A(1,1) ,A(1,2)],'Color','b','LineWidth',1);   
 end

%subplot(1,1,1);


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
% imshow(surfluimage,[]);
% title('surffff image.');
% subplot(2,3,5);
% bar(bins,his,'b');
% title('the distribution of superpixel size');
% ylabel('# of superpixels');
% xlabel('superpixel sizes in pixel');
% scnsize = get(0,'ScreenSize');
% set(gcf,'OuterPosition',scnsize);