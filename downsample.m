function [ output_args ] = downsample( input_args )
%LR can be obtained
[m n ch]  =   size(input_args);

if ch == 3
       im_r = input_args(:, :, 1);
       im_g = input_args(:, :, 2);
       im_b = input_args(:, :, 3);
 
       im_r_l=fdownsampling(im_r);
       im_g_l=fdownsampling(im_g);
       im_b_l=fdownsampling(im_b);
       
       [nrow, ncol] = size(im_r_l);


output_args = zeros([nrow, ncol, 3]);
output_args(:, :, 1) = im_r_l;
output_args(:, :, 2) = im_g_l;
output_args(:, :, 3) = im_b_l;
output_args=uint8(output_args);
else

       output_args=fdownsampling(input_args);
       output_args=uint8(output_args);
end

end

