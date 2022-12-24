function boundary  = g(x,y,i)
if(i==1)
    boundary= u_exact(0,y); % x=0
elseif(i==2)
    boundary= u_exact(1,y); % x=1
elseif(i==3)
    boundary = u_exact(x,0); % y=0
elseif(i==4)
    boundary= u_exact(x,1); % y=1
end
