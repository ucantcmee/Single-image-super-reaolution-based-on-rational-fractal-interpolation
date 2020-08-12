function [ average_final ] = fractal_interpolation3(z)
%FRACTAL_INTERPOLATION Summary of this function goes here

s=[  0.02  0.02 0.02  0.02;
     0.02  0.02 0.02  0.02;
     0.02  0.02 0.02  0.02;
     0.02  0.02 0.02  0.02;
   ];

alpha0=[
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
   ];                 

beta0=[
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
       0.8  0.8  0.8  0.8;
   ];
gama0=[
       1.2  1.2  1.2  1.2;
       1.2  1.2  1.2  1.2;
       1.2  1.2  1.2  1.2;
       1.2  1.2  1.2  1.2;
   ];
alpha=[
       0.2  0.2  0.6  0.6;
       0.2  0.2  0.6  0.6;
       0.6  0.6  0.6  0.6;
       0.6  0.6  0.6  0.6;
   ];                 

beta=[
       0.2  0.2  0.6  0.6;
       0.2  0.2  0.6  0.6;
       0.6  0.6  0.6  0.6;
       0.6  0.6  0.6  0.6;
   ];              

gama=[
       2  2  1.5  1.5;
       2  2  1.5  1.5;
       1.5  1.5  1.5  1.5;
       1.5  1.5  1.5  1.5;
   ];
kk=1;
M=3;
N=3;
x=1:0.5:M;
y=1:0.5:N;
 [m,n]=size(z);
 ppp=m;
 for p=1:kk
     for i=1:n-1
         for k=1:4
             for j=1:m-1
              a=(x(i+1)-x(i))/(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1)); 
              b=(x(4*floor((i-1)/4)+5)*x(i)-x(4*floor((i-1)/4)+1)*x(i+1))/(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1));
              c=(y(j+1)-y(j))/(y(4*floor((j-1)/4)+5)-y(4*floor((j-1)/4)+1)); 
              d=(y(4*floor((j-1)/4)+5)*y(j)-y(4*floor((j-1)/4)+1)*y(j+1))/(y(4*floor((j-1)/4)+5)-y(4*floor((j-1)/4)+1));
              for t=1:4
                  sita(k,x,n);
                  eta(t,y,m);
                  dx(i,j,x,ppp,z);
                  dy(i,j,y,ppp,z);
                 % sss=dx(k,t,x,ppp,z)-dy(k,t,y,ppp,z);
                 a00=((1-sita(k,x,n))^2 * (alpha0(i,j)+sita(k,x,n)*gama0(i,j))*(1-eta(t,y,m))^2 * (alpha(i,j)+eta(t,y,m)*gama(i,j))) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 a01=((1-sita(k,x,n))^2 * (alpha0(i,j)+sita(k,x,n)*gama0(i,j))*eta(t,y,m)^2 * (beta(i,j)+(1-eta(t,y,m))*gama(i,j)))/ (((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 a10=(sita(k,x,n)^2 * (beta0(i,j)+(1-sita(k,x,n))*gama0(i,j))*(1-eta(t,y,m))^2 *(alpha(i,j)+eta(t,y,m)*gama(i,j)))/(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 a11=(sita(k,x,n)^2 * (beta0(i,j)+(1-sita(k,x,n))*gama0(i,j))*(eta(t,y,m))^2 *(beta(i,j)+(1-eta(t,y,m))*gama(i,j)))/ (((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 
                 b00=(sita(k,x,n) *(1-sita(k,x,n))^2 *alpha0(i,j)*(1-eta(t,y,m))^2 * (alpha(i,j)+eta(t,y,m)*gama(i,j))) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 b01= (sita(k,x,n) * (1-sita(k,x,n))^2 *alpha0(i,j)*eta(t,y,m)^2 * (beta(i,j)+(1-eta(t,y,m))*gama(i,j))) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 b10=-(sita(k,x,n)^2 * (1-sita(k,x,n)) *beta0(i,j)*(1-eta(t,y,m))^2 * (alpha(i,j)+eta(t,y,m)*gama(i,j))) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 b11=-((sita(k,x,n))^2* (1-sita(k,x,n)) *beta0(i,j)* eta(t,y,m)^2 *(beta(i,j)+(1-eta(t,y,m))*gama(i,j))) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                  
                 c00=((1-sita(k,x,n))^2* (alpha0(i,j)+sita(k,x,n)*gama0(i,j)) *eta(t,y,m) *(1-eta(t,y,m))^2*alpha(i,j)) / (((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 c01=-((1-sita(k,x,n))^2*(alpha0(i,j)+sita(k,x,n)*gama0(i,j))* eta(t,y,m)^2* (1-eta(t,y,m))*beta(i,j)) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 c10=(sita(k,x,n)^2*(beta0(i,j)+(1-sita(k,x,n))*gama0(i,j)) * eta(t,y,m)*(1-eta(t,y,m))^2*alpha(i,j)) /(((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 c11=-(sita(k,x,n)^2*(beta0(i,j)+(1-sita(k,x,n))*gama0(i,j)) * eta(t,y,m)^2 *(1-eta(t,y,m))*beta(i,j)) / (((1-sita(k,x,n))^2*alpha0(i,j)+sita(k,x,n)*(1-sita(k,x,n))*gama0(i,j)+sita(k,x,n)^2*beta0(i,j))*((1-eta(t,y,m))^2*alpha(i,j)+eta(t,y,m)*(1-eta(t,y,m))*gama(i,j)+eta(t,y,m)^2*beta(i,j)));
                 
                 
                  xx((i-1)*4+t)=a*x(4*floor((i-1)/4)+t)+b; 
                  yy((j-1)*4+t)=c*y(4*floor((j-1)/4)+t)+d; 
                    aa((i-1)*4+k,(j-1)*4+t)=alpha(k,t);
                    bb((i-1)*4+k,(j-1)*4+t)=beta(k,t);     
                    ss((i-1)*4+k,(j-1)*4+t)=s(k,t);
                  zz((i-1)*4+k,(j-1)*4+t)=s(i,j)*z(k,t) + a00*(z(i,j)-s(i,j)*z(1,1)) + a10*(z(i+1,j)-s(i,j)*z(ppp,1)) + a01*(z(i,j+1)-s(i,j)*z(1,ppp)) + a11*(z(i+1,j+1)-s(i,j)*z(ppp,ppp)) ...
                    + b00*((x(i+1)-x(i))*dx(i,j,x,ppp,z)-s(i,j)*(x(ppp)-x(1))*dx(1,1,x,ppp,z))+ b10*((x(i+1)-x(i))*dx(i+1,j,x,ppp,z)-s(i,j)*(x(ppp)-x(1))*dx(ppp,1,x,ppp,z)) ...
                    + b01*((x(i+1)-x(i))*dx(i,j+1,x,ppp,z)-s(i,j)*(x(ppp)-x(1))*dx(1,ppp,x,ppp,z)) + b11*((x(i+1)-x(i))*dx(i+1,j+1,x,ppp,z)-s(i,j)*(x(ppp)-x(1))*dx(ppp,ppp,x,ppp,z)) ...
                    +c00*((y(j+1)-y(j))*dy(i,j,y,ppp,z)-s(i,j)*(y(ppp)-y(1))*dy(1,1,y,ppp,z)) + c10*((y(j+1)-y(j))*dy(i+1,j,y,ppp,z)-s(i,j)*(y(ppp)-y(1))*dy(ppp,1,y,ppp,z))...
                    + c01*((y(j+1)-y(j))*dy(i,j+1,y,ppp,z)-s(i,j)*(y(ppp)-y(1))*dy(1,ppp,y,ppp,z)) + c11*((y(j+1)-y(j))*dy(i+1,j+1,y,ppp,z)-s(i,j)*(y(ppp)-y(1))*dy(ppp,ppp,y,ppp,z));
              end
             end
         end
     end
     zx=z(n,1:m);
     zy=z(1:n,m);
     for i=1:n-1
         for j=1:4
           sita(j,x,n);
           eta(j,y,m);
                 aa00=((1-sita(j,x,n))^2*(alpha0(i,m-1)+sita(j,x,n)*gama0(i,m-1)))/((1-sita(j,x,n))^2*alpha0(i,m-1)+sita(j,x,n)*(1-sita(j,x,n))*gama0(i,m-1)+sita(j,x,n)^2*beta0(i,m-1));
                 aa10=(sita(j,x,n)^2*((1-sita(j,x,n))*gama0(i,m-1)+beta0(i,m-1)))/((1-sita(j,x,n))^2*alpha0(i,m-1)+sita(j,x,n)*(1-sita(j,x,n))*gama0(i,m-1)+sita(j,x,n)^2*beta0(i,m-1));
                 bb10=(sita(j,x,n)*(1-sita(j,x,n))^2*alpha0(i,m-1))/((1-sita(j,x,n))^2*alpha0(i,m-1)+sita(j,x,n)*(1-sita(j,x,n))*gama0(i,m-1)+sita(j,x,n)^2*beta0(i,m-1));
                 bb11=-((sita(j,x,n))^2*(1-sita(j,x,n))*beta(i,m-1))/((1-sita(j,x,n))^2*alpha0(i,m-1)+sita(j,x,n)*(1-sita(j,x,n))*gama0(i,m-1)+sita(j,x,n)^2*beta0(i,m-1));
                 
                 zzx((i-1)*4+j)= s(i,m-1)*zx(4*floor((i-1)/4)+j)+ aa00*(zx(i)-s(i,m-1)*zx(4*floor((i-1)/4)+1))+ aa10*(zx(i+1)-s(i,m-1)*zx(4*floor((i-1)/4)+5))+ bb10*((x(i+1)-x(i))*dx(i,m,x,m,z)-s(i,m-1)*(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1))*dx(1,m,x,m,z))...
                 +bb11*((x(i+1)-x(i))*dx(i+1,m,x,n,z)-s(i,m-1)*(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1))*dx(n,m,x,n,z)); 
             
                  a00=((1-eta(j,y,m))^2*(alpha(n-1,i)+eta(j,y,m)*gama(n-1,i)))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  a01=(eta(j,y,m)^2*((1-eta(j,y,m))*gama(n-1,i)+beta(n-1,i)))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  c10=(eta(j,y,m)*(1-eta(j,y,m))^2*alpha(n-1,i))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  c11=-((1-eta(j,y,m))*eta(j,y,m)^2*beta(n-1,i))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
            
                  zzy((i-1)*4+j)= s(n-1,i)*zy(4*floor((i-1)/4)+j)+ a00*(zy(i)-s(n-1,i)*zy(4*floor((i-1)/4)+1))+ a01*(zy(i+1)-s(n-1,i)*zy(4*floor((i-1)/4)+5))+ c10*((y(i+1)-y(i))*dy(n,i,y,m,z)-s(n-1,i)*(y(4*floor((i-1)/4)+5)-y(4*floor((i-1)/4)+1))*dy(n,1,y,m,z))...
                 +c11*((y(i+1)-y(i))*dy(n,i+1,y,n,z)-s(n-1,i)*(y(4*floor((i-1)/4)+5)-y(4*floor((i-1)/4)+1))*dy(n,m,y,n,z)); 
         end
     end
     zz=cat(2,zz,zzy');
    tmp=cat(2,zzx,z(n,m));
    zz=cat(1,zz,tmp);
     
    xx=cat(2,xx,xx(4*(m-1))+xx(2)-xx(1));
    yy=cat(2,yy,yy(4*(m-1))+yy(2)-yy(1));
   
 %=========================================
 x=xx;
 y=yy;
 z=zz;
 alpha=aa;
 beta=bb;
 s=ss;
 [m,n]=size(z);
 end
 %pixel mapping
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
average_final=zeros(12,12);
% temp1=zeros(1,3);
% temp2=zeros(1,3);
% temp3=zeros(1,3);
% temp4=zeros(1,3);
temp5=zeros(1,4);
temp6=zeros(1,4);
temp7=zeros(1,4);
temp8=zeros(1,4);
row=1;
col=1;
for i=1:4:16
    for j=1:4:16
        average_final(row,col)=zz(i,j);

         temp5(1,1)=zz(i,j);
         temp5(1,2)=zz(i,j+1);
         temp5(1,3)=zz(i,j+2);
         temp5(1,4)=zz(i,j+3);
         max5=max(temp5);
         min5=min(temp5);
         ave1=(max5+min5)/2;
         average_final(row,col+1)=ave1;
         
         temp6(1,1)=zz(i,j+1);
         temp6(1,2)=zz(i,j+2);
         temp6(1,3)=zz(i,j+3);
         temp6(1,4)=zz(i,j+4);
         max6=max(temp6);
         min6=min(temp6);
         ave2=(max6+min6)/2;
         average_final(row,col+2)=ave2;

         temp7(1,1)=zz(i,j);
         temp7(1,2)=zz(i+1,j);
         temp7(1,3)=zz(i+2,j);
         temp7(1,4)=zz(i+3,j);
         max7=max(temp7);
         min7=min(temp7);
         ave3=(max7+min7)/2;
         average_final(row+1,col)=ave3;
         
         temp8(1,1)=zz(i+1,j);
         temp8(1,2)=zz(i+2,j);
         temp8(1,3)=zz(i+3,j);
         temp8(1,4)=zz(i+4,j);
         max8=max(temp8);
         min8=min(temp8);
         ave4=(max8+min8)/2;
         average_final(row+2,col)=ave4;
         
    
         
         average_final(row+1,col+1)= (zz(i+1,j+1)+zz(i+1,j+2)+zz(i+2,j+1)+zz(i+2,j+2))/4;
         average_final(row+1,col+2)= (zz(i+1,j+2)+zz(i+1,j+3)+zz(i+2,j+2)+zz(i+2,j+3))/4;
         average_final(row+2,col+1)= (zz(i+2,j+1)+zz(i+2,j+2)+zz(i+3,j+1)+zz(i+3,j+2))/4;
         average_final(row+2,col+2)= (zz(i+2,j+2)+zz(i+2,j+3)+zz(i+3,j+2)+zz(i+3,j+3))/4;

        col=col+3;
        if(col==13)
            row=row+3;
            col=1;
        end
    end
end
%�ڶ��ַ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tem=demo_SME_ImageZoom(zz);
% average_final=zeros(6,6);
% average_final(1,1)=tem(1,1);
% average_final(1,4)=tem(1,9);
% average_final(4,1)=tem(9,1);
% average_final(4,4)=tem(9,9);
% 
% average_final(1,2)=tem(1,4);
% average_final(1,3)=tem(1,6);
% average_final(2,1)=tem(4,1);
% average_final(2,2)=tem(4,4);
% average_final(2,3)=tem(4,6);
% average_final(3,1)=tem(6,1);
% average_final(3,2)=tem(6,4);
% average_final(3,3)=tem(6,6);
% average_final(4,2)=tem(9,4);
% average_final(4,3)=tem(9,6);
% average_final(5,1)=tem(12,1);
% average_final(5,2)=tem(12,4);
% average_final(5,3)=tem(12,6);
% average_final(6,1)=tem(14,1);
% average_final(6,2)=tem(14,4);
% average_final(6,3)=tem(14,6);
% 
% average_final(1,5)=tem(1,12);
% average_final(1,6)=tem(1,14);
% average_final(2,4)=tem(4,9);
% average_final(2,5)=tem(4,12);
% average_final(2,6)=tem(4,14);
% average_final(3,4)=tem(6,9);
% average_final(3,5)=tem(6,12);
% average_final(3,6)=tem(6,14);
% average_final(4,5)=tem(9,12);
% average_final(4,6)=tem(9,14);
% average_final(5,4)=tem(12,9);
% average_final(5,5)=tem(12,12);
% average_final(5,6)=tem(12,14);
% average_final(6,4)=tem(14,9);
% average_final(6,5)=tem(14,12);
% average_final(6,6)=tem(14,14);
% row=1;
% col=1;
% j=1;
% k=1;
% for ii=2:1:7
%     if(ii/2==0)
%         for i=2:1:7
%             if(i/2==0)
%                 average_final(row,col)=tem(j,k);
%                 k=k+3;
%                 col=col+1;                             
%             else
%                average_final(row,col)=tem(j,k);
%                k=k+2;
%                col=col+1;
%         
%             end
%         end 
%         col=1;
%         row=row+1;
%         k=1;
%         j=j+3;
%     else
%         for i=2:1:7
%             if(i/2==0)
%                 average_final(row,col)=tem(j,k);
%                 k=k+3;
%                 col=col+1;
%               
%             else
%                average_final(row,col)=tem(j,k);
%                k=k+2;
%                col=col+1;
%               
%             end
%         end 
%         col=1;
%         row=row+1;
%         k=1;
%         j=j+2;
%     end
% end
 %�����ַ���
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  average_final=zeros(12,12);
% temp1=zeros(1,3);
% temp2=zeros(1,3);
% temp3=zeros(1,3);
% temp4=zeros(1,3);
% 
% temp5=zeros(1,4);
% temp6=zeros(1,4);
% temp7=zeros(1,4);
% temp8=zeros(1,4);
% 
% temp9=zeros(1,3);
% temp10=zeros(1,3);
% temp11=zeros(1,3);
% temp12=zeros(1,3);
% 
% temp13=zeros(1,3);
% temp14=zeros(1,3);
% temp15=zeros(1,3);
% temp16=zeros(1,3);
% 
% temp17=zeros(1,3);
% temp18=zeros(1,3);
% temp19=zeros(1,3);
% temp20=zeros(1,3);
% row=1;
% col=1;
% for i=1:4:16
%     for j=1:4:16
%         average_final(row,col)=zz(i,j);
% 
%          temp5(1,1)=zz(i,j);
%          temp5(1,2)=zz(i,j+1);
%          temp5(1,3)=zz(i,j+2);
%          temp5(1,4)=zz(i,j+3);
%          max5=max(temp5);
%          min5=min(temp5);
%          ave1=(max5+min5)/2;
%          average_final(row,col+1)=ave1;
%          
%          temp6(1,1)=zz(i,j+1);
%          temp6(1,2)=zz(i,j+2);
%          temp6(1,3)=zz(i,j+3);
%          temp6(1,4)=zz(i,j+4);
%          max6=max(temp6);
%          min6=min(temp6);
%          ave2=(max6+min6)/2;
%          average_final(row,col+2)=ave2;
% 
%          temp7(1,1)=zz(i,j);
%          temp7(1,2)=zz(i+1,j);
%          temp7(1,3)=zz(i+2,j);
%          temp7(1,4)=zz(i+3,j);
%          max7=max(temp7);
%          min7=min(temp7);
%          ave3=(max7+min7)/2;
%          average_final(row+1,col)=ave3;
%          
%          temp8(1,1)=zz(i+1,j);
%          temp8(1,2)=zz(i+2,j);
%          temp8(1,3)=zz(i+3,j);
%          temp8(1,4)=zz(i+4,j);
%          max8=max(temp8);
%          min8=min(temp8);
%          ave4=(max8+min8)/2;
%          average_final(row+2,col)=ave4;
%          
%     
%          %%�м�4�������
%         temp1(1,1)=zz(i+1,j);
%         temp1(1,2)=zz(i+1,j+1);
%         temp1(1,3)=zz(i+1,j+2);
%         hor=var(temp1);
%        % hor=abs(zz(i+2,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+2,j+3));
%         temp2(1,1)=zz(i,j+1);
%         temp2(1,2)=zz(i+1,j+1);
%         temp2(1,3)=zz(i+2,j+1);
%         ver=var(temp2);
% %         ver=abs(zz(i+1,j+2)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+2));
%         temp3(1,1)=zz(i+2,j);
%         temp3(1,2)=zz(i+1,j+1);
%         temp3(1,3)=zz(i,j+2);
%         dia_45=var(temp3);
% %         dia_45=abs(zz(i+1,j+3)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+1));
%         temp4(1,1)=zz(i,j);
%         temp4(1,2)=zz(i+1,j+1);
%         temp4(1,3)=zz(i+2,j+2);
%         dia_145=var(temp4);
% %         dia_145=abs(zz(i+1,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+3));
%           %�������С�ķ���
%         m=min(hor,ver);
%         m1=min(m,dia_45);
%         m2=min(m1,dia_145);
%         if(m2==hor)
%               average_final(row+1,col+1)=(zz(i+1,j)+zz(i+1,j+1)+zz(i+1,j+2))/3;
%         else if(m2==ver)
%               average_final(row+1,col+1)=(zz(i,j+1)+zz(i+1,j+1)+zz(i+2,j+1))/3;
%             else if(m2==dia_45)
%                     average_final(row+1,col+1)=(zz(i+2,j)+zz(i+1,j+1)+zz(i,j+2))/3;
%                 else
%                       average_final(row+1,col+1)=(zz(i,j)+zz(i+1,j+1)+zz(i+2,j+2))/3;
%                 end
%             end
%         end
%         
%         temp9(1,1)=zz(i+1,j+2);
%         temp9(1,2)=zz(i+1,j+3);
%         temp9(1,3)=zz(i+1,j+4);
%         hor1=var(temp9);
%        % hor=abs(zz(i+2,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+2,j+3));
%         temp10(1,1)=zz(i,j+3);
%         temp10(1,2)=zz(i+1,j+3);
%         temp10(1,3)=zz(i+2,j+3);
%         ver1=var(temp10);
% %         ver=abs(zz(i+1,j+2)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+2));
%         temp11(1,1)=zz(i+2,j+2);
%         temp11(1,2)=zz(i+1,j+3);
%         temp11(1,3)=zz(i,j+4);
%         dia_451=var(temp11);
% %         dia_45=abs(zz(i+1,j+3)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+1));
%         temp12(1,1)=zz(i,j+2);
%         temp12(1,2)=zz(i+1,j+3);
%         temp12(1,3)=zz(i+2,j+4);
%         dia_1451=var(temp12);
% %         dia_145=abs(zz(i+1,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+3));
%           %�������С�ķ���
%         mm=min(hor1,ver1);
%         mm1=min(mm,dia_451);
%         mm2=min(mm1,dia_1451);
%         if(mm2==hor1)
%                average_final(row+1,col+2)=(zz(i+1,j+2)+zz(i+1,j+3)+zz(i+1,j+4))/3;
%         else if(mm2==ver1)
%                 average_final(row+1,col+2)=(zz(i,j+3)+zz(i+1,j+3)+zz(i+2,j+3))/3;
%             else if(mm2==dia_451)
%                     average_final(row+1,col+2)=(zz(i+2,j+2)+zz(i+1,j+3)+zz(i,j+4))/3;
%                 else
%                        average_final(row+1,col+2)=(zz(i,j+2)+zz(i+1,j+3)+zz(i+2,j+4))/3;
%                 end
%             end
%         end
%         
%         temp13(1,1)=zz(i+3,j);
%         temp13(1,2)=zz(i+3,j+1);
%         temp13(1,3)=zz(i+3,j+2);
%         hor2=var(temp13);
%        % hor=abs(zz(i+2,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+2,j+3));
%         temp14(1,1)=zz(i+2,j+1);
%         temp14(1,2)=zz(i+3,j+1);
%         temp14(1,3)=zz(i+4,j+1);
%         ver2=var(temp14);
% %         ver=abs(zz(i+1,j+2)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+2));
%         temp15(1,1)=zz(i+4,j);
%         temp15(1,2)=zz(i+3,j+1);
%         temp15(1,3)=zz(i+2,j+2);
%         dia_452=var(temp15);
% %         dia_45=abs(zz(i+1,j+3)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+1));
%         temp16(1,1)=zz(i+2,j);
%         temp16(1,2)=zz(i+3,j+1);
%         temp16(1,3)=zz(i+4,j+2);
%         dia_1452=var(temp16);
% %         dia_145=abs(zz(i+1,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+3));
%           %�������С�ķ���
%         mmm=min(hor2,ver2);
%         mmm1=min(mmm,dia_452);
%         mmm2=min(mmm1,dia_1452);
%         if(mmm2==hor2)
%                average_final(row+2,col+1)=(zz(i+3,j)+zz(i+3,j+1)+zz(i+3,j+2))/3;
%         else if(mmm2==ver2)
%               average_final(row+2,col+1)=(zz(i+2,j+1)+zz(i+3,j+1)+zz(i+4,j+1))/3;
%             else if(mmm2==dia_452)
%                    average_final(row+2,col+1)=(zz(i+4,j)+zz(i+3,j+1)+zz(i+2,j+2))/3;
%                 else
%                     average_final(row+2,col+1)=(zz(i+2,j)+zz(i+3,j+1)+zz(i+4,j+2))/3;
%                 end
%             end
%         end
%         
%         temp17(1,1)=zz(i+3,j+2);
%         temp17(1,2)=zz(i+3,j+3);
%         temp17(1,3)=zz(i+3,j+4);
%         hor3=var(temp17);
%        % hor=abs(zz(i+2,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+2,j+3));
%         temp18(1,1)=zz(i+2,j+3);
%         temp18(1,2)=zz(i+3,j+3);
%         temp18(1,3)=zz(i+4,j+3);
%         ver3=var(temp18);
% %         ver=abs(zz(i+1,j+2)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+2));
%         temp19(1,1)=zz(i+4,j+2);
%         temp19(1,2)=zz(i+3,j+3);
%         temp19(1,3)=zz(i+2,j+4);
%         dia_453=var(temp19);
% %         dia_45=abs(zz(i+1,j+3)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+1));
%         temp20(1,1)=zz(i+2,j+2);
%         temp20(1,2)=zz(i+3,j+3);
%         temp20(1,3)=zz(i+4,j+4);
%         dia_1453=var(temp20);
% %         dia_145=abs(zz(i+1,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+3));
%           %�������С�ķ���
%         mmmm=min(hor3,ver3);
%         mmmm1=min(mmmm,dia_453);
%         mmmm2=min(mmmm1,dia_1453);
%         if(mmmm2==hor3)
%               average_final(row+2,col+2)=(zz(i+3,j+2)+zz(i+3,j+3)+zz(i+3,j+4))/3;
%         else if(mmmm2==ver3)
%              average_final(row+2,col+2)=(zz(i+2,j+3)+zz(i+3,j+3)+zz(i+4,j+3))/3;
%             else if(mmmm2==dia_453)
%                  average_final(row+2,col+2)=(zz(i+4,j+2)+zz(i+3,j+3)+zz(i+2,j+4))/3;
%                 else
%                     average_final(row+2,col+2)=(zz(i+2,j+2)+zz(i+3,j+3)+zz(i+4,j+4))/3;
%                 end
%             end
%         end
% %          average_final(row+1,col+1)= (zz(i+1,j+1)+zz(i+1,j+2)+zz(i+2,j+1)+zz(i+2,j+2))/4;
% %          average_final(row+1,col+2)= (zz(i+1,j+2)+zz(i+1,j+3)+zz(i+2,j+2)+zz(i+2,j+3))/4;
% %          average_final(row+2,col+1)= (zz(i+2,j+1)+zz(i+2,j+2)+zz(i+3,j+1)+zz(i+3,j+2))/4;
% %          average_final(row+2,col+2)= (zz(i+2,j+2)+zz(i+2,j+3)+zz(i+3,j+2)+zz(i+3,j+3))/4;
% 
%         col=col+3;
%         if(col==13)
%             row=row+3;
%             col=1;
%         end
%     end
% end
end

