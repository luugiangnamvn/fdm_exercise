function u_ex=exact_solution(x,k);

if(k==1)
    u_ex=x^10+x^5-x+1;
end

if (k==2)
    u_ex=sin(17*x*pi/2)+1;
end
