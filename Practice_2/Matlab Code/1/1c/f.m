function f1=f(x,y);
% f1=2*x^3*y^5*exp(-2*x^2*y^2) - (3*x^3*y*exp(-2*x^2*y^2))/2 - (3*x*y^3*exp(-2*x^2*y^2))/2 + 2*x^5*y^3*exp(-2*x^2*y^2);
f1 = ((4*pi*x)^2+(2*pi)^2)*sin(2*pi*x^2)*sin(2*pi*y)-4*pi*cos(2*pi*x^2)*sin(2*pi*y);
end