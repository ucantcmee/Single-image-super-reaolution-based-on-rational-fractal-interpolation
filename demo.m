

clear all;
close all;
clc;
[file,path]=uigetfile('*.*','input image')



%     out_dir            =    '../results';
%     im_dir             =    '../data/';
%     im_name            =    'bird_GT.bmp';
%     
    im        =   imread(file);
    [m n ch]  =   size(im);
    


%  im=imresize(im,[28 128]);
figure,imshow(im)



%  downsample
    scale=2;         %% Scaling factors: 2, 3, 4
    %select down-sampling method
     dmethod=3;
    if dmethod == 1
        % bicubic down-sampling
        LR        =    imresize(im, 1/scale, 'bicubic');
    end
    if dmethod ==2
        % average down-samping
         LR  =  downsample(im);
    end
    if dmethod==3
        % LR images are obtained by down-sampling the HR images directly along both the horizontal and vertical directions by a factor of 2, 3, or 4.
     LR        =    im(1:scale:end,1:scale:end,:);  
    end
    
    
    if  ch == 3
        % change color space, work on illuminance only
        im_l_ycbcr = rgb2ycbcr(LR);
        im_l_y = im_l_ycbcr(:, :, 1);
        im_l_cb = im_l_ycbcr(:, :, 2);
        im_l_cr = im_l_ycbcr(:, :, 3);
        im_l_y=double(im_l_y);
        im_l_cb=double(im_l_cb);
        im_l_cr=double(im_l_cr);
        %expand the metrix
        [m n]=size(im_l_y);
        II(1:m,1:n) = im_l_y;
        II(m+1,:) = 2.*II(m,:) - II(m-1,:);
        II(:,n+1) = 2.*II(:,n) - II(:,n-1);
        II(m+2,:) = 2.*II(m+1,:) - II(m,:);
        II(:,n+2) = 2.*II(:,n+1) - II(:,n);
        II(m+3,:) =2.*II(m+2,:)-II(m+1,:);
        II(:,n+3) =2.*II(:,n+2)-II(:,n+1);
        II(m+4,:) =2.*II(m+3,:)-II(m+2,:);
        II(:,n+4) =2.*II(:,n+3)-II(:,n+2);
    
        figure,imshow(mat2gray(II))
        
        
        
        % image super-resolution
        im_h_y = main_function(II,m,n,scale);
        
        
        
        
    % upscale the chrominance simply by "bicubic" 
    [nrow, ncol] = size(im_h_y);
    im_h_cb = imresize(im_l_cb, [nrow, ncol], 'bicubic');
    im_h_cr = imresize(im_l_cr, [nrow, ncol], 'bicubic');
    
    im_h_ycbcr = zeros([nrow, ncol, 3]);
    im_h_ycbcr(:, :, 1) = im_h_y;
    im_h_ycbcr(:, :, 2) = im_h_cb;
    im_h_ycbcr(:, :, 3) = im_h_cr;
    im_h = ycbcr2rgb(uint8(im_h_ycbcr));
    if dmethod == 1 || dmethod == 2
        Image=double(im_h);
        Image1=Image;
        Half_size=2;
        F_size=2*Half_size+1;
        G_Filter=fspecial('gaussian',F_size,F_size/6);
        Image_Filter = imfilter(Image1, G_Filter,'conv');
        Image_Diff=Image-Image_Filter;
        Image_out=Image_Diff+0;
        Image1=Image+Image_out;
        im_h=uint8(Image1);
    end
    % show the images
    figure, imshow(im_h);
    fname            =   strcat('SRRFL_', im_name);    
%     imwrite(im_h, fullfile(out_dir, fname));  
    else
    
        I=double(LR);
        %expand the metrix
        [m n]=size(I);
        II(1:m,1:n) = I;
        II(m+1,:) = 2.*II(m,:) - II(m-1,:);
        II(:,n+1) = 2.*II(:,n) - II(:,n-1);
        II(m+2,:) = 2.*II(m+1,:) - II(m,:);
        II(:,n+2) = 2.*II(:,n+1) - II(:,n);
        II(m+3,:) =2.*II(m+2,:)-II(m+1,:);
        II(:,n+3) =2.*II(:,n+2)-II(:,n+1);
        II(m+4,:) =2.*II(m+3,:)-II(m+2,:);
        II(:,n+4) =2.*II(:,n+3)-II(:,n+2);
    
        % image super-resolution
        im_h= main_function(II,m,n,scale);
        im_h=uint8(im_h);
        if dmethod == 1 || dmethod == 2
        Image=double(im_h);
        Image1=Image;
        Half_size=2;
        F_size=2*Half_size+1;
        G_Filter=fspecial('gaussian',F_size,F_size/6);
        Image_Filter = imfilter(Image1, G_Filter,'conv');
        Image_Diff=Image-Image_Filter;
        Image_out=Image_Diff+0;
        Image1=Image+Image_out;
        im_h=uint8(Image1);
        end
        % show the images
        figure, imshow(im_h);
        fname            =   strcat('SRRFL_', im_name);    
%     imwrite(im_h, fullfile(out_dir, fname));  
    end
   
