function [bmap] = seg2bmap(seg,width,height)
% [bmap] = seg2bmap(seg)
%
% From a segmentation, compute a binary boundary map with 1 pixel wide
% boundaries.  The boundary pixels are offset by 1/2 pixel towards the
% origin from the actual segment boundary.
%
% INPUTS
%	seg		Segments labeled from 1..k.
%	[width]		Width of desired bmap, <= size(seg,2)
%	[height]	Height of desired bmap, <= size(seg,1)
%
% OUTPUTS
%	bmap		Binary boundary map.
%
% David Martin <dmartin@eecs.berkeley.edu>
% January 2003

if nargin<3,   
  [height,width] = size(seg);
end
[h,w] = size(seg);    

% check width and height
ar1 = width / height;
ar2 = w / h;
if width>w | height>h | abs(ar1-ar2)>0.01,    % HEIGHT VE WIDTH DE�?ERLERİ KONTROL EDİLİR.
  error(sprintf('Can''t convert %dx%d seg to %dx%d bmap.',w,h,width,height));
end

e = zeros(size(seg));     % DE�?İ�?KENLERİN İSİMLERİNİN ANLAMI BİLİNMİYOR????
s = zeros(size(seg));     % LABELS IN BOYUTU KADAR BİR 0 MATRİSİ OLU�?TURULUYOR.
se = zeros(size(seg));    

e(:,1:end-1) = seg(:,2:end);    % LABELS IN 2.SÜTUNUNDAN BA�?LAYARAK, e MATRİSİNE DE�?ERLER KOPYALANIYOR.
s(1:end-1,:) = seg(2:end,:);    % LABELS IN 2.SATIRINDAN BA�?LAYARAK, s MATRİSİNE DE�?ERLER KOPYALANIYOR.
se(1:end-1,1:end-1) = seg(2:end,2:end);   % LABELS IN 2.SATIR 2. SÜTUNUNDAN BA�?LAYARAK, se MATRİSİNE DE�?ERLER KOPYALANIYOR.

b = (seg~=e | seg~=s | seg~=se);    % kontrol sonucuna göre boundary matrisi oluşturuluyor.
b(end,:) = (seg(end,:)~=e(end,:));
b(:,end) = (seg(:,end)~=s(:,end));
b(end,end) = 0;
  
if w==width & h==height,
  
  bmap = b;

else
  
  bmap = zeros(height,width);
  for x = 1:w,
    for y = 1:h,
      if b(y,x),
        j = 1+floor((y-1)*height/h);
        i = 1+floor((x-1)*width/w);
        bmap(j,i) = 1;
      end
    end
  end

end

