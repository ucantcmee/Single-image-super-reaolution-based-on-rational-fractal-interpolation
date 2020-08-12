function values

x=0:1:4;
y=0:1:4;


%函数值
f=[ 4  3  6  5  4;
    3  5  3  4  3;
    5  3  4  3  5;
    3  3  3  4  5;
    4  3  4  5  3;
]


%插商值
for i=1:4;
    for j=1:5;
        deltax(i,j)=(f(i+1,j)-f(i,j))/(x(i+1)-x(i))
    end
end
for i=1:5;
    for j=1:4;
        deltay(i,j)=(f(i,j+1)-f(i,j))/(y(j+1)-y(j))
    end
end


%偏导数值
for j=1:5;
    dx(1,j)=deltax(1,j)-(deltax(2,j)-deltax(1,j))/2
end
for i=2:4;
    for j=1:5;
        dx(i,j)=(deltax(i,j)+deltax(i-1,j))/2
    end
end
for j=1:5;
    dx(5,j)=deltax(4,j)+(deltax(4,j)-deltax(4,j))/2
end

for i=1:5;
    dy(i,1)=deltay(i,1)-(deltay(i,2)-deltay(i,1))/2
end
for i=1:5;
    for j=2:4;
        dy(i,j)=(deltay(i,j)+deltay(i,j-1))/2
    end
end
for i=1:5;
    dy(i,5)=deltay(i,4)+(deltay(i,4)-deltay(i,3))/2
end


%尺度因子取值
for i=1:4
    for j=1:4
        hni=(x(5)-x(1))/(x(i+1)-x(i));
      lmj=(y(5)-y(1))/(y(j+1)-y(j));
      deltn1x=(f(5,1)-f(1,1))/(x(i+1)-x(i))
      deltnmx=(f(5,5)-f(1,5))/(x(i+1)-x(i));
      delt1my=(f(1,5)-f(1,1))/(y(j+1)-y(j));
      deltnmy=(f(5,5)-f(5,1))/(y(j+1)-y(j));
    end
end
N=5;
M=5;

for i=1:4;
    for j=1:4;
        s1(i,j)=min(dx(i,j)/(hni*dx(1,1)),dx(i+1,j)/(hni*dx(N,1)));
        s2(i,j)=min(dx(i,j+1)/(hni*dx(1,M)),dx(i+1,j+1)/(hni*dx(N,M)));
        s3(i,j)=min(deltax(i,j)/(deltn1x),deltax(i,j+1)/(deltnmx));
        
        s4(i,j)=min(dy(i,j)/(lmj*dy(1,1)),dy(i,j+1)/(lmj*dy(1,M)));
        s5(i,j)=min(dy(i+1,j)/(lmj*dy(N,1)),dy(i+1,j+1)/(lmj*dy(N,M)));
        s6(i,j)=min(deltay(i,j)/(delt1my),deltay(i+1,j)/(deltnmy));
        
        ss1(i,j)=min(s1(i,j),s2(i,j));
        ss2(i,j)=min(s3(i,j),s4(i,j));
        ss3(i,j)=min(s5(i,j),s6(i,j));
        
        sss1(i,j)=min(ss1(i,j),ss2(i,j));
        sss2(i,j)=min(0.2,ss3(i,j));
        
        s(i,j)=0.5*min(sss1(i,j),sss2(i,j))
    end
end

%形状参数取值

for i=1:4;
    for j=1:4;
        
        %gamma*/alpha*  取值
        
        xga1(i,j)=max(2*(dx(i,j)-s(i,j)*hni*dx(1,1))/(deltax(i,j)-s(i,j)*deltn1x),2*(dx(i,j+1)-s(i,j)*hni*dx(1,M))/(deltax(i,j+1)-s(i,j)*deltnmx));
        xga(i,j)=0.5+max(xga1(i,j),0)
        xga(i,5)=0.5+max(2*(dx(i,5)-s(i,4)*hni*dx(1,5))/(deltax(i,5)-s(i,4)*deltnmx),0)
        
        %gamma*/beta*  取值
        
        xgb1(i,j)=max(2*(dx(i+1,j)-s(i,j)*hni*dx(N,1))/(deltax(i,j)-s(i,j)*deltn1x),2*(dx(i+1,j+1)-s(i,j)*hni*dx(N,M))/(deltax(i,j+1)-s(i,j)*deltnmx));
        xgb(i,j)=0.5+max(xgb1(i,j),0)
        xgb(i,5)=0.5+max(2*(dx(i+1,5)-s(i,4)*hni*dx(N,5))/(deltax(i,5)-s(i,4)*deltnmx),0)
        
        %gamma/alpha  取值
        
        yga1(1,j)=max(2*(dy(1,j)-s(1,j)*lmj*dy(1,1))/(deltay(1,j)-s(1,j)*delt1my),2*(dy(2,j)-s(1,j)*lmj*dy(N,1))/(deltay(2,j)-s(1,j)*deltnmy));
        yga2(2,j)=max(2*(dy(2,j)-s(2,j)*lmj*dy(1,1))/(deltay(2,j)-s(2,j)*delt1my),2*(dy(3,j)-s(2,j)*lmj*dy(N,1))/(deltay(3,j)-s(2,j)*deltnmy));
        yga3(3,j)=max(2*(dy(3,j)-s(3,j)*lmj*dy(1,1))/(deltay(3,j)-s(3,j)*delt1my),2*(dy(4,j)-s(3,j)*lmj*dy(N,1))/(deltay(4,j)-s(3,j)*deltnmy));
        yga4(4,j)=max(2*(dy(4,j)-s(4,j)*lmj*dy(1,1))/(deltay(4,j)-s(4,j)*delt1my),2*(dy(5,j)-s(4,j)*lmj*dy(N,1))/(deltay(5,j)-s(4,j)*deltnmy));
        yyga1(i,j)=max(yga1(1,j),yga2(2,j));
        yyga2(i,j)=max(yga3(3,j),yga4(4,j));
        yga(i,j)=0.0+max(yyga1(i,j),yyga2(i,j))
        
        %gamma/beta  取值
        
        ygb1(1,j)=max(2*(dy(1,j+1)-s(1,j)*lmj*dy(1,M))/(deltay(1,j)-s(1,j)*delt1my),2*(dy(2,j+1)-s(1,j)*lmj*dy(N,M))/(deltay(2,j)-s(1,j)*deltnmy));
        ygb2(2,j)=max(2*(dy(2,j+1)-s(2,j)*lmj*dy(1,M))/(deltay(2,j)-s(2,j)*delt1my),2*(dy(3,j+1)-s(2,j)*lmj*dy(N,M))/(deltay(3,j)-s(2,j)*deltnmy));
        ygb3(3,j)=max(2*(dy(3,j+1)-s(3,j)*lmj*dy(1,M))/(deltay(3,j)-s(3,j)*delt1my),2*(dy(4,j+1)-s(3,j)*lmj*dy(N,M))/(deltay(4,j)-s(3,j)*deltnmy));
        ygb4(4,j)=max(2*(dy(4,j+1)-s(4,j)*lmj*dy(1,M))/(deltay(4,j)-s(4,j)*delt1my),2*(dy(5,j+1)-s(4,j)*lmj*dy(N,M))/(deltay(5,j)-s(4,j)*deltnmy));
        yygb1(i,j)=max(ygb1(1,j),ygb2(2,j));
        yygb2(i,j)=max(ygb3(3,j),ygb4(4,j));
        ygb(i,j)=0.0+max(yygb1(i,j),yygb2(i,j))
    end
end




