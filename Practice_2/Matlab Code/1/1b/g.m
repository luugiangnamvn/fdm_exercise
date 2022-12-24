function boundary  = g(x,y,i)
if(i==1)
    boundary= y^2; % x=0
elseif(i==2)
    boundary= 1+y+y^2; % x=1
elseif(i==3)
    boundary = x^2; % y=0
elseif(i==4)
    boundary= 1+x+x^2; % y=1
end
