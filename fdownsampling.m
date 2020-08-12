function [ LR ] = fdownsampling(input)

input=double(input);
[m n]=size(input);
LR=zeros(m/2,n/2);
ii=1;
jj=1;
for i=1:2:m-1
    for j=1:2:n-1
        LR(ii,jj)=(input(i,j)+input(i,j+1)+input(i+1,j)+input(i+1,j+1))/4;
        if(jj==n/2)
            ii=ii+1;
            jj=1;
        else
            jj=jj+1;
        end
    end
end

end

