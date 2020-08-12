function fd=blanket(A)
M=double(A);
m=size(M);
e=101;
U=zeros(m(1),m(2),e);
B=zeros(m(1),m(2),e);

for i=1:1:m(1)
    for j=1:1:m(2)
        U(i,j,1)=M(i,j);
        B(i,j,1)=M(i,j);
    end
end

for p=2:1:e-1
    for i=1:1:m(1)
        for j=1:1:m(2)
            if i==1
                a=i;
            else a=i-1;
            end
            if j==1
                b=j;
            else b=j-1;
            end
            if i==m(1)
                c=i;
            else c=i+1;
            end
            if j==m(2)
                d=j;
            else d=j+1;
            end
  U(i,j,p)=max(U(i,j,p-1)+1,max(U(i,j,p-1),max(max(U(a,j,p-1),U(i,b,p-1)),max(U(c,j,p-1),U(i,d,p-1)))));
  B(i,j,p)=min(B(i,j,p-1)-1,min(B(i,j,p-1),min(min(B(a,j,p-1),B(i,b,p-1)),min(B(c,j,p-1),B(i,d,p-1)))));
  
        end
    end
end
V=zeros(51,1);
V1=zeros(m(1),m(2));
for p=50:1:e-1
    for i=1:1:m(1)
        for j=1:1:m(2)
            V1(i,j)=U(i,j,p)-B(i,j,p);
        end
    end
    V(p-49,1)=sum(sum(V1));
end
A=zeros(51,1); 
for p=50:1:e-1
    A(p-49,1)=V(p-49,1)/(2*p);
end


X=zeros(3,1);
Y=zeros(3,1);
 kk=1;fd=0;
for p=50:1:e-1   
     X(kk,1)=log(p);
     Y(kk,1)=log(A(p-49,1));
    if kk==3
        p1=polyfit(X,Y,1);
        fdtemp=2-p1(1);
        kk=0;
        if p==52
            fdt1=fdtemp;
        end
        fd=fd+fdtemp-fdt1;
    end
   kk=kk+1;
end




X=zeros(51,1);
Y=zeros(51,1);
for p=50:1:e-1
    X(p-49,1)=log(p);
    Y(p-49,1)=log(A(p-49,1));
end
p=polyfit(X,Y,1);
fd=2-p(1);
  