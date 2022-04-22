close all
clear all
clc


% read RGB image
RGB = imread('D:\SWAVE2.0\SWAVE2.0_DATA_ORGANIZED\placenta021\43879034_324085928401947_3278674177162215424_n.jpg');

figure,imagesc(RGB), colormap jet


h= imfreehand;
Y = wait(h);
h1= imfreehand;
Y1 = wait(h1);
Y = [Y; Y1];
h1= imfreehand;
Y2 = wait(h1);
Y = [Y; Y2];

% imagesc(sum2),colormap parula
% h1= imfreehand;
% Y1 = wait(h1);
% Y = [Y; Y1];

MRI = dicomread('D:\SWAVE2.0\SWAVE2.0_DATA_ORGANIZED\placenta021\MRI\1.2.840.113619.2.408.4738430.15806453.23679.1539702061.802\1.2.840.113619.2.408.4738430.15806453.21877.1539702577.50');
MRI = imrotate(MRI,180);

MRIinfo = dicominfo('D:\SWAVE2.0\SWAVE2.0_DATA_ORGANIZED\placenta021\MRI\1.2.840.113619.2.408.4738430.15806453.23679.1539702061.802\1.2.840.113619.2.408.4738430.15806453.21877.1539702577.50');


figure, imagesc(MRI)

h= imfreehand;X = wait(h);

h1= imfreehand;
X1 = wait(h1);
X = [X;X1];
h1= imfreehand;
X2 = wait(h1);
X = [X;X2];




img1 = im2double((RGB));
refim = im2double(MRI);

% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;
opt.scale =0; 
opt.rot = 1 ;
opt.corresp=1;
opt.tol = 1e-7;
% show every iteration
% registering Y to X
[Transform,C]=cpd_register(Y,X,opt);

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
result(:,:,1)=interp2(img1(:,:,1),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)
result(:,:,2)=interp2(img1(:,:,2),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)
result(:,:,3)=interp2(img1(:,:,3),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)


% Show the result
figure,imshow(img1); title('Source image');
figure,imshow(refim); title('Reference image');
figure,imagesc(result); title('Result image');

dicomwrite(result, 'result.dcm', MRIinfo);

%% intrinsic and world coordinate

RA = imref2d(size(MRI),MRIinfo.PixelSpacing(2),MRIinfo.PixelSpacing(1));
figure
imshow(MRI,RA,'DisplayRange',[0 512])
figure
imshow(result,RA,'DisplayRange',[0 512])