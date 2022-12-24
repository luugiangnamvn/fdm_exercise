function f1=f(x,y);
% f1= -(- y^12 + y^9)*(- 90*x^8 + 20*x^3) - (- x^10 + x^5)*(- 132*y^10 + 72*y^7);
f1 = ((4*pi*x)^2+(2*pi)^2)*sin(2*pi*x^2)*sin(2*pi*y)-4*pi*cos(2*pi*x^2)*sin(2*pi*y);
end