%Partial Matlab code for "Bundled Camera Paths for Video Stabilization" (SIGGRAPH 2013)
%Implementation of motion model estimation.
%1. As-similar-as-possible warping.
%2. Local homography estimation on mesh quads.
%require vision tool box for detectSURFFeatures, or you may want to use
%your own features. (N x 2)
% clear all;
% clc;

addpath('mesh');
addpath('RANSAC');

I1 = imread(strcat('out',string(910),'.png')); %right image. Will be warped
I2 = imread('out906.png');
[height,width,~] = size(I1);
%I1_alt = imwarp(I1,tform)
Tr = [[1,0,(1)*width],[0,1,0]];
I1alt = imtranslate(I1,[(-1)*width,0],'FillValues',0,'OutputView','full');
I1 = imtranslate(I1alt,[1*width,0],'FillValues',0,'OutputView','full');

I2alt = imtranslate(I2,[(-1)*width,0],'FillValues',0,'OutputView','full');
I2 = imtranslate(I2alt,[(1)*width,0],'FillValues',0,'OutputView','full');

fprintf('detect surf features...');
[I1_features,I2_features]=SURF(I1,I2);
plot(I1_features(:,1),I1_features(:,2),'r*');
plot(I2_features(:,1),I2_features(:,2),'b*');
fprintf('[DONE]');

if length(I1_features) < 10
    error('not enough matched features');
    return;
end

[height,width,~] = size(I1);
%3x3 mesh
quadWidth = width/(2^3);
quadHeight = height/(2^3);

% %4x4 mesh
% quadWidth = width/(2^4);
% quadHeight = height/(2^4);

lamda = 1; %mesh more rigid if larger value. [0.2~5]
asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform
%fprintf(homos);
gap = 0;
I1warp = asap.Warp(I1,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
fprintf("DR#");
gap = 0;
I1warp = asap.Warp(I1,gap);                    

%access local homography
[h,w,~,~] = size(homos);
for i=1:h-1
    for j=1:w-1
       H(:,:) = homos(i,j,:,:);
       fprintf('Quad=[%d %d]\n',i,j);
       H
    end
end
imwrite(I1warp,"right_warped.png");
%fancy_func
I1warp_gray = rgb2gray(I1warp);
I2_gray = rgb2gray(I2);
I1warp_mask = im2uint8(imbinarize(I1warp_gray,(1/255.0)));
I2_mask = im2uint8(imbinarize(I2_gray,(1/255.0)));
common_mask = bitand(I2_mask,I1warp_mask);
red_I1warp_common = bitand(common_mask,I1warp(:,:,1));
gre_I1warp_common = bitand(common_mask,I1warp(:,:,2));
blu_I1warp_common = bitand(common_mask,I1warp(:,:,3));
I1warp_common = I1warp;
I1warp_common(:,:,1) = red_I1warp_common;
I1warp_common(:,:,2) = gre_I1warp_common;
I1warp_common(:,:,3) = blu_I1warp_common;
I1warp_extra = imsubtract(I1warp,I1warp_common);
combine = imadd(I2,I1warp_extra);
imshow(combine);
imwrite(I2_mask,"red_I2_common.png");
imwrite(combine,"detect.png");

