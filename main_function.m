function [ B ] = main_function( II,m,n,scale)
%MAIN_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
if scale==2 
    for i = 1:2:m
    for j = 1:2:n
         block = double(II(i:i+4,j:j+4));
         T=blanket(block);
         if(T>2.12)
             %rational fractal interpolation
             c=fractal_interpolation2(block);
         else
             %rational interpolation
             c=fractal_interpolation2_s0(block);            
         end         
       if(j==1)
            A=c(1:4,1:4);
       else
            A =cat(2,A,c(1:4,1:4));
       end
   end
    if(i==1)
         B=A;
    else
     B=cat(1,B,A);
    end
   end
end   
if scale==3
   for i = 1:2:m
   for j = 1:2:n
         block = double(II(i:i+4,j:j+4));
         T=blanket(block);
         if(T>2.21)
             %rational fractal interpolation
             c=fractal_interpolation3(block);
         else
             %rational interpolation
             c=fractal_interpolation3_s0(block);            
         end         
       if(j==1)
            A=c(1:6,1:6);
       else
            A =cat(2,A,c(1:6,1:6));
       end
   end
    if(i==1)
         B=A;
    else
     B=cat(1,B,A);
    end
   end  
%    B=B(1:end-4,1:end-4);
end
 if scale==4
   for i = 1:2:m
   for j = 1:2:n
         block = double(II(i:i+4,j:j+4));
         T=blanket(block);
         if(T>2.26)
             c=fractal_interpolation4(block);
         else
             c=fractal_interpolation4_s0(block);            
         end         
       if(j==1)
            A=c(1:8,1:8);
       else
            A =cat(2,A,c(1:8,1:8));
       end
   end
    if(i==1)
         B=A;
    else
     B=cat(1,B,A);
    end
    end
 end

end




