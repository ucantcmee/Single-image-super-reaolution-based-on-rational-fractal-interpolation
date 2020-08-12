function [ average_final ] = fractal_interpolation2_s0(z)


s=[    0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
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
                 %�����һhang
                 zzx((i-1)*4+j)= s(i,m-1)*zx(4*floor((i-1)/4)+j)+ aa00*(zx(i)-s(i,m-1)*zx(4*floor((i-1)/4)+1))+ aa10*(zx(i+1)-s(i,m-1)*zx(4*floor((i-1)/4)+5))+ bb10*((x(i+1)-x(i))*dx(i,m,x,m,z)-s(i,m-1)*(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1))*dx(1,m,x,m,z))...
                 +bb11*((x(i+1)-x(i))*dx(i+1,m,x,n,z)-s(i,m-1)*(x(4*floor((i-1)/4)+5)-x(4*floor((i-1)/4)+1))*dx(n,m,x,n,z)); 
             
                  a00=((1-eta(j,y,m))^2*(alpha(n-1,i)+eta(j,y,m)*gama(n-1,i)))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  a01=(eta(j,y,m)^2*((1-eta(j,y,m))*gama(n-1,i)+beta(n-1,i)))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  c10=(eta(j,y,m)*(1-eta(j,y,m))^2*alpha(n-1,i))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  c11=-((1-eta(j,y,m))*eta(j,y,m)^2*beta(n-1,i))/((1-eta(j,y,m))^2*alpha(n-1,i)+eta(j,y,m)*(1-eta(j,y,m))*gama(n-1,i)+eta(j,y,m)^2*beta(n-1,i));
                  %�����һ��
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
 %ȡƽ��ֵ+����
average_final=zeros(8,8);
temp1=zeros(1,3);
temp2=zeros(1,3);
temp3=zeros(1,3);
temp4=zeros(1,3);
temp5=zeros(1,5);
temp6=zeros(1,5);
row=1;
col=1;
for i=1:4:16
    for j=1:4:16
        average_final(row,col)=zz(i,j);
%          average_final(row,col+1)=(zz(i,j+1)+zz(i,j+2)+zz(i,j+3))/3;%�Ϸ�
         temp5(1,1)=zz(i,j);
         temp5(1,2)=zz(i,j+1);
         temp5(1,3)=zz(i,j+2);
         temp5(1,4)=zz(i,j+3);
         temp5(1,5)=zz(i,j+4);
         max5=max(temp5);
         min5=min(temp5);
         ave1=(max5+min5)/2;
         average_final(row,col+1)=ave1;
%         average_final(row+1,col)=(zz(i+1,j)+zz(i+2,j)+zz(i+3,j))/3;%���
         temp6(1,1)=zz(i,j);
         temp6(1,2)=zz(i+1,j);
         temp6(1,3)=zz(i+2,j);
         temp6(1,4)=zz(i+3,j);
         temp6(1,5)=zz(i+4,j);
         max6=max(temp6);
         min6=min(temp6);
         ave2=(max6+min6)/2;
         average_final(row+1,col)=ave2;
        %�����Ÿ���
        temp1(1,1)=zz(i+2,j+1);
        temp1(1,2)=zz(i+2,j+2);
        temp1(1,3)=zz(i+2,j+3);
        hor=var(temp1);
       % hor=abs(zz(i+2,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+2,j+3));
        temp2(1,1)=zz(i+1,j+2);
        temp2(1,2)=zz(i+2,j+2);
        temp2(1,3)=zz(i+3,j+2);
        ver=var(temp2);
%         ver=abs(zz(i+1,j+2)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+2));
        temp3(1,1)=zz(i+1,j+3);
        temp3(1,2)=zz(i+2,j+2);
        temp3(1,3)=zz(i+3,j+1);
        dia_45=var(temp3);
%         dia_45=abs(zz(i+1,j+3)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+1));
        temp4(1,1)=zz(i+1,j+1);
        temp4(1,2)=zz(i+2,j+2);
        temp4(1,3)=zz(i+3,j+3);
        dia_145=var(temp4);
%         dia_145=abs(zz(i+1,j+1)-zz(i+2,j+2))+abs(zz(i+2,j+2)-zz(i+3,j+3));
          %�������С�ķ���
        m=min(hor,ver);
        m1=min(m,dia_45);
        m2=min(m1,dia_145);
        if(m2==hor)
              average_final(row+1,col+1)=(zz(i+2,j+1)+zz(i+2,j+2)+zz(i+2,j+3))/3;
        else if(m2==ver)
              average_final(row+1,col+1)=(zz(i+1,j+2)+zz(i+2,j+2)+zz(i+3,j+2))/3;
            else if(m2==dia_45)
                    average_final(row+1,col+1)=(zz(i+1,j+3)+zz(i+2,j+2)+zz(i+3,j+1))/3;
                else
                      average_final(row+1,col+1)=(zz(i+1,j+1)+zz(i+2,j+2)+zz(i+3,j+3))/3;
                end
            end
        end
        col=col+2;
        if(col==9)
            row=row+2;
            col=1;
        end
    end
end

end

