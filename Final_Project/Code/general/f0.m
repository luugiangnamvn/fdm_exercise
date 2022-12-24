function f = f0(x,t)
% f = - 2*exp(-2*t)*(x^2 + 1) - 2*a1(x,t)*exp(-2*t) - 2*x*exp(-2*t)*(b1(x,t) + fa1(x,t)) - exp(-2*t)*(c1(x,t) + fb1(x,t))*(x^2 + 1);
f = (16*pi^2*a1(x,t)*exp(-2*t)*cos((4*pi*x)/3))/9 - exp(-2*t)*cos((4*pi*x)/3)*(c1(x,t) + fb1(x,t)) - 2*exp(-2*t)*cos((4*pi*x)/3) + (4*pi*exp(-2*t)*sin((4*pi*x)/3)*(b1(x,t) + fa1(x,t)))/3;