clear all;
 close all;
 clc;
 [file,path]=uigetfile('*.*','input image')
 im=imread(file);
%  input=imresize(im,[256 256]);

  figure,imshow(im)
 scale=2;  
    input=imresize(im,[128 128]);
% 
% 
% 
% %%%%%%%%%%%dividing image into 16x16 blocks%%%%%%%%%%%
% [rows ,columns ,numberOfColorBands] = size(input);
% blockSizeR = 16; % Rows in block.
% blockSizeC = 16; % Columns in block.       
% wholeBlockRows = rows / blockSizeR;
% blockVectorR = [blockSizeR * ones(1, wholeBlockRows)];
% wholeBlockCols = columns / blockSizeC;
% blockVectorC = [blockSizeC * ones(1, wholeBlockCols)];
% 
% if numberOfColorBands > 1
% 	% It's a color image.
% 	ca = mat2cell(input, blockVectorR, blockVectorC, numberOfColorBands);
%     M=mat2cell(input, blockVectorR, blockVectorC, numberOfColorBands);
% else
% 	ca = mat2cell(input, blockVectorR, blockVectorC);
%     M=mat2cell(input, blockVectorR, blockVectorC);
% end
% % Now display all the blocks.
% 
% 
% plotIndex = 1;
% numPlotsR = size(ca, 1);
% numPlotsC = size(ca, 2);
% 
% 
% 
% 
% for r = 1 : numPlotsR
% 	for c = 1 : numPlotsC
%         figure(2),
% 		subplot(numPlotsR, numPlotsC, plotIndex);
% 		rgbBlock = ca{r,c};
% 		imshow(rgbBlock);
%        
%         
%         
%          plotIndex = plotIndex + 1;
%      
%     end
%  
%  
% 
% end
% 
% 
% 
% k=1;
% for i=1:8
%  for j=1:8
%     for m=1:16
%         for n=1:16
%           
%         Y=ca{i,j}(m,n)/(16*16);
%         delx{i,j}(m,n)=(ca{i,j}(m,n)-Y);
%         
%       
%         end
%     end
%     
%       figure(3),
% 		subplot(8,8,k);
%         imshow(delx{i,j})
%         
%        
%  
%       
%     k=k+1; 
%      
%  end    
% end
% 
% 
% 
% 
% 
% 
%  image=cell2mat(delx);
%  figure,imshow(image)
% 
%  
%  ima1=255-image;
%  figure,imshow(ima1,[]) 
     
 LR=im(1:scale:end,1:scale:end,:);  
 
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
    figure, imshow(im_h);
    
    
   