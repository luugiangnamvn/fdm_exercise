function u_ex=exact_solution(x,k);
if (k==1)
    u_ex=x^4+x^2+1;
end

if(k==2)
    u_ex=cos(9*x*pi/2)+3;
end