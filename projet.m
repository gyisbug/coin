clear all;
close all;
img = (imread('images/1/1.jpg'));
figure(1);subplot(3,3,1);imshow(img,'initialmagnification','fit');title('Original image ') 


% Convert to grayscale image
img_gray = rgb2gray(img);
figure(1);subplot(3,3,2);imshow(img_gray,'initialmagnification','fit');title('Grayscale ') 

% Correct background drift
SE=strel('disk',100);
I=imclose(img_gray,SE);
F=I-img_gray;
img_gray=255-F;
figure(1);subplot(3,3,3);imshow(img_gray,'initialmagnification','fit');title('Correct background ')



% Edge detection
img_edge = edge(img_gray, 'roberts', 0.04,'both');
figure(1),subplot(3,3,4), imshow(img_edge,'initialmagnification','fit'), title('edge detection roberts');

% [VSFAT Threshold]=edge(img_gray, 'sobel','vertical');
% img_edge =edge(img_gray,'sobel',Threshold);
% figure(2),subplot(2,2,2), imshow(img_edge ,'initialmagnification','fit'), title('edge detection sobel');

%Closed operation image to make the coin edges complete without holes
se_disk = strel('disk',20);
img_close = imclose(img_edge, se_disk);
figure(1),subplot(3,3,5), imshow(img_close ,'initialmagnification','fit'), title('close')



%Remove small objects (noise) and fill complete objects
img_fill = imfill(img_close, 'holes');
img1=bwareaopen (img_fill, 500);
figure(1), subplot(3,3,6),imshow(img1,'initialmagnification','fit'), title('fill holes')


%clear cut between objects
img1=tse_imsplitobjects(img1);
figure(1), subplot(3,3,7), imshow(img1,'initialmagnification','fit'), title('fill holes')

%Measure the parameters of the object
[B, L] = bwboundaries(img1);
figure(1), subplot(3,3,8), imshow(label2rgb(L, @jet, [.5 .5 .5]),'initialmagnification','fit'), title('boundaries')
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

stats = regionprops(L,img(:,:,1),...
    'Area','Centroid','Orientation','EquivDiameter','MeanIntensity');
threshold = 0.80; % For differentiating coins from matches based on an objects circularity

coinCentroids = [];
coinArea = [];
coinRatios = [];
coinEquivDiameter=[];

for k = 1:length(B)
    boundary = B{k};
    delta_sq = diff(boundary).^2;
    perimeter = sum(sqrt(sum(delta_sq,2)));
    area = stats(k).Area;
    metric = 4*pi*area/perimeter^2;
    metric_string = sprintf('%2.2f',metric);
    angle_string = sprintf('%2.2f',stats(k).Orientation);
    centroid = stats(k).Centroid;
%     if metric > threshold
%         % Object is round, therefore a coin
        coinCentroids = [coinCentroids; centroid];
        coinArea=[coinArea; area];
        coinRatios = [coinRatios; stats(k).EquivDiameter/area];
        coinEquivDiameter=[coinEquivDiameter; stats(k).EquivDiameter];

    plot(centroid(1),centroid(2),'ko');
%     text(boundary(1,2)-35,boundary(1,1)+13,angle_string,'Color','y',...
%       'FontSize',14,'FontWeight','bold');

end

%coins count
c1=0;
c2=0;
c3=0;
c4=0;
c5=0;
c6=0;
c7=0;
c8=0;

for k = 1:length(coinEquivDiameter)
    if((stats(k).EquivDiameter)>810&&(stats(k).EquivDiameter)<900)
        c8=c8+1;
    elseif((stats(k).EquivDiameter)>740&&(stats(k).EquivDiameter)<775)
        c7=c7+1;
    elseif((stats(k).EquivDiameter)>775&&(stats(k).EquivDiameter)<810)
        c6=c6+1;
    elseif((stats(k).EquivDiameter)>705&&(stats(k).EquivDiameter)<740)
        c5=c5+1;
    elseif((stats(k).EquivDiameter)>620&&(stats(k).EquivDiameter)<670)
        c4=c4+1;
    elseif((stats(k).EquivDiameter)>670&&(stats(k).EquivDiameter)<705)
        c3=c3+1;
    elseif((stats(k).EquivDiameter)>580&&(stats(k).EquivDiameter)<620)
        c2=c2+1;
    elseif((stats(k).EquivDiameter)>500&&(stats(k).EquivDiameter)<550)
        c1=c1+1;
    
    end

end
sum1=c8*2+c7+c6*0.5+c5*0.2+c4*0.1+c3*0.05+c2*0.02+c1*0.01
tabcoins=[c1;c2;c3;c4;c5;c6;c7;c8]