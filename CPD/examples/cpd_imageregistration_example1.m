% CPD 2D image registration example.
% The idea is to extract feature points from both images beforehand.
% Register the extracted point sets. And then apply the estimated
% transformation to the images.

% clear all; %close all; clc;
% load cpd_data_geometry.mat; % load two images
% im = imcrop(im);
% % Extract feature points. Here we use a simple edge detector.
% % Most likely you'll have to use a more advanced feature detector here.
% [j,i]=find(edge(im,'canny'));
% X=[i j]; % first point set
% 
% [j,i]=find(edge(refim,'canny'));
% Y=[i j]; % second point set
close all
res = 0.251e-6;
A2_W = 119624;
A2_H = 82267;
A3_W = 82267;
A3_H = 74675;
W =A2_W;
H = A2_H;
IW = 0.0381;
IH = 0.030184;
AW = res*W;
AH = res*H;
mul_factor = 15890; %round(480/0.030184)
% refim1 = imread ('C:\Users\Farah\Dropbox\EOSC\Placental Data\B-mode Images\786064\vol2\786064_vol_2_5.png');
refim = a(:,:,35);
% refim1 = sum1;
refim = imcrop(refim,[18,1,604,480]);
figure,imagesc(refim), colormap jet

size(refim)

refim = imresize (refim, [round((IH)*mul_factor),round((IW)*mul_factor)] );
% close all
figure, imagesc(refim),colormap gray
h= imfreehand;Y = wait(h);
% sum2 = imcrop(sum1,[18,1,604,480]);
 imagesc(refim),colormap gray
h1= imfreehand;
Y1 = wait(h1);
Y = [Y; Y1];

% imagesc(sum2),colormap parula
% h1= imfreehand;
% Y1 = wait(h1);
% Y = [Y; Y1];


patho = imread('C:\Farah\LFE\PP\786045_A2_P2.tif');
patho = imcrop(patho,[1,1,size(patho,2)-5, size(patho,1)]);

aa= im2bw(patho,0.9);
figure, imshow((patho))
patho = imcrop(patho,[min(find(aa(500,:)==0)),1, max(find(aa(500,:)==0))- min(find(aa(500,:)==0)), size(patho,1)]);
% patho = imcrop(patho);
figure, imshow(patho)


patho = imresize (patho, [round((AH)*mul_factor),round((AW)*mul_factor)] );


% imshow(refim)
% h2= imfreehand;
% Y1 = wait(h2);
% Y = [Y; Y1];

figure, imshow(patho)
h= imfreehand;X = wait(h);
imshow(patho)
h1= imfreehand;
X1 = wait(h1);
X = [X;X1];
% imshow(patho)
% h1= imfreehand;
% X1 = wait(h1);
% X = [X;X1];
% imshow(patho)
% h2= imfreehand;
% X1 = wait(h2);
% X = [X;X1];

img1 = im2double(rgb2gray(patho));
refim = im2double(refim);

% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;
opt.scale =0; 
opt.rot = 0 ;
opt.corresp=1;
% show every iteration
% registering Y to X
[Transform,C]=cpd_register(X,Y,opt);

% Now lets apply the found transformation to the image
% Create a dense grid of points
[M,N]=size(refim);
[x,y]=meshgrid(1:N,1:M); 
grid=[x(:) y(:)];

% Transform the grid according to the estimated transformation
T=cpd_transform(grid, Transform);

% Interpolate the image
Tx=reshape(T(:,1),[M N]);
Ty=reshape(T(:,2),[M N]);
result=interp2(img1,Tx,Ty);   % evaluate image Y at the new gird (through interpolation)

% Show the result
figure,imshow(img1); title('Source image');
figure,imshow(refim); title('Reference image');
figure,imagesc(result); title('Result image');
