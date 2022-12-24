function u_ex=exact_solution(x,k);
if (k==1)
    u_ex=x^3-3/2*x^2 + 1/4;
end

if(k==2)
    u_ex=3*x^(10) -10*x^3 +49/22;
end