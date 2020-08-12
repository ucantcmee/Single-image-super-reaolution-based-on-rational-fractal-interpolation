function res = dx(i,j,x,ppp,z)
z=z;
if (i==1)                            
    res=(4*z(2,j)-3*z(1,j)-z(3,j))/(2*(x(i+1)-x(i)));
else if(i<=(ppp-1))
        res=(z(i+1,j)-z(i-1,j))/(2*(x(i+1)-x(i)));
    else
        res=(3*z(ppp,j)-4*z(ppp-1,j)+z(ppp-2,j))/(2*(x(i)-x(i-1)));
    end
end
   