addpath('mesh');
addpath('RANSAC');

num_images = 50;
width = 2;
start = 550;  %205 to 1350 - 285 with skip of 4
skip = 2;
I1 = imread(strcat("../final_frames/out",string(start),".png"));
I2 = imread(strcat("../final_frames/out",string(start+skip),".png"));
[height,offset,~] = size(I1);
I1alt = imtranslate(I1,[(-1)*(width-1)*offset,0],'FillValues',0,'OutputView','full');
%I1altalt = imtranslate(I1alt,[1*offset,0],'FillValues',0,'OutputView','full');
I2alt = imtranslate(I2,[(-1)*(width-1)*offset,0],'FillValues',0,'OutputView','full');
%I2altalt = imtranslate(I2alt,[(1)*offset,0],'FillValues',0,'OutputView','full');
[result,asap] = stitch2(I2alt,I1alt);
cnt = 2;
check = 0;
for i = (start+2*skip):skip:(start+(skip)*num_images)
    if cnt > 15
        cnt = 0;
        width = width+1;
        result = imtranslate(result,[(-1)*offset,0],'FillValues',0,'OutputView','full');
        check = 1;
        end
    fprintf(string(i));
    fprintf(" ");
    fprintf(string(cnt));
    fprintf('\n');
    I = imread(strcat("../final_frames/out",string(i),".png"));
    Ialt = imtranslate(I,[(-1)*(width-2)*offset,0],'FillValues',0,'OutputView','full');
    %Ialtalt = imtranslate(Ialt,[1*offset,0],'FillValues',0,'OutputView','full');
    I_warp = asap.Warp(Ialt,0);
    if check > 0
        I_warp = imtranslate(I_warp,[(-1)*offset,0],'FillValues',0,'OutputView','full');
        check = 0;
    end
    [result_temp,asap_temp] = stitch2(I_warp,result);
    if strcmp(result_temp,"bypass")
        fprintf("bypass \n");
        i = i-skip+1;
        continue
    end
    result = result_temp;
    asap = asap_temp;
    cnt = cnt+1;
end
imshow(result);
imwrite(result,"result.png");

function [result,asap] = stitch2(I1,I2)
    [height,width,~] = size(I1);
    %I1_alt = imwarp(I1,tform)

    %fprintf('detect surf features...');
    [I1_features,I2_features]=SURF(I1,I2);
    if strcmp(I1_features,"bypass")
        result = "bypass";
        asap = "bypass";
        return
    end
    plot(I1_features(:,1),I1_features(:,2),'r*');
    plot(I2_features(:,1),I2_features(:,2),'b*');
    %fprintf('[DONE]');

    if length(I1_features) < 3
        error('not enough matched features');
        return;
    end
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
    I1warp = asap.Warp(I1,gap); 
    if strcmp(I1warp,"bypass")
        result = "bypass";
        asap = "bypass";
        return 
    end
    %warp source image to target image
    I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
    %fprintf("DR#");
    gap = 0;
    I1warp = asap.Warp(I1,gap);                    
    %access local homography
    %[h,w,~,~] = size(homos);
    %for i=1:h-1
    %    for j=1:w-1
    %       H(:,:) = homos(i,j,:,:);
    %       %fprintf('Quad=[%d %d]\n',i,j);
    %       H
    %    end
    %end
    result = fancy_func(I1warp,I2);
    return
end

function result = fancy_func(I1warp,I2)
    I1warp_gray = rgb2gray(I1warp); %I1warp is the right warped image which will be overlayed onto the left
    I2_gray = rgb2gray(I2);
    I1warp_mask = im2uint8(imbinarize(I1warp_gray,(1/255.0)));
    I2_mask = im2uint8(imbinarize(I2_gray,(1/255.0)));
    common_mask = bitand(I2_mask,I1warp_mask);
    red_I2_common = bitand(common_mask,I2(:,:,1));
    gre_I2_common = bitand(common_mask,I2(:,:,2));
    blu_I2_common = bitand(common_mask,I2(:,:,3));
    I2_common = I2;
    I2_common(:,:,1) = red_I2_common;
    I2_common(:,:,2) = gre_I2_common;
    I2_common(:,:,3) = blu_I2_common;
    I2_extra = imsubtract(I2,I2_common);
    result = imadd(I2_extra,I1warp);
end



