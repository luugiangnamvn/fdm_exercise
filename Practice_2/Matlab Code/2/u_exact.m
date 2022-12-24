function uex=u_exact(x,y);
% uex = -1/8*x*y*exp(-2*x^2*y^2);
uex=sin(2*pi*x^2)*sin(2*pi*y)+(x*y)/2 + 1;
end