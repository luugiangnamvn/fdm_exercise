function boundary  = g(x,y,i,h)
if(i==1)
    boundary= u_exact(x,0); % y=0
elseif(i==2)
    boundary= u_exact(0,y); % x=0
elseif(i==3)
    boundary = u_exact(x,1); % y=1
elseif(i==4)
    boundary= (3*u_exact(1,y)-4*u_exact(1-h,y)+u_exact(1-2*h,y))/(2*h); % dg/dx(y=1)
end
