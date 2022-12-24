function f = f0(a,b,x,t)
% f = - 2*exp(-2*t)*(x^2 + 1) - 2*a*exp(-2*t) - 2*b*x*exp(-2*t);
f = (25*pi^2*a*exp(-2*t)*cos((5*pi*x)/3))/9 - 2*exp(-2*t)*cos((5*pi*x)/3) + (5*pi*b*exp(-2*t)*sin((5*pi*x)/3))/3;
