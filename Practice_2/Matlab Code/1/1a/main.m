%% Solve equation -u''(x)=f(x) with the Dirichlet boundary condition 
clear all
close all
format long
clc
%% Initial informations
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;

%% Mesh
N=10;% number of mesh points of first mesh
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);


%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    number_mesh_point(inumber_mesh)=N*N;
    del_x=(bx-ax)/N;
    del_y=(by-ay)/N;
    
%% Create mesh point    
    x=zeros(N+1,1);
    for i=1:N+1
        x(i)=(i-1)*del_x;
    end
    y=zeros(N+1,1);
    for i=1:N+1
        y(i)=(i-1)*del_y;
    end
    %% Tao ma tran A
    
    A=sparse((N-1)*(N-1),(N-1)*(N-1));
    B=sparse(N-1,N-1);
    
    for i=1:N-1
        if (i==1)
            B(i,i) = 4;
            B(i,i+1) = -1;
        elseif (i==N-1)
            B(i,i) = 4;
            B(i,i-1)=-1;
        else
            B(i,i-1)=-1;
            B(i,i)=4;
            B(i,i+1)=-1;
        end
    end
    I=eye(N-1);
    for i=1:N-1
        if(i==1)
            A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1)) = B;
            A((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1)) = -I;
        elseif (i==N-1)
            A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1)) = B;
            A((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1)) = -I;
        else
            A((i-1)*(N-1)+1:i*(N-1),(i-2)*(N-1)+1:(i-1)*(N-1)) = -I;
            A((i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1)) = B;
            A((i-1)*(N-1)+1:i*(N-1),i*(N-1)+1:(i+1)*(N-1)) = -I;
        end
    end
    A = A/(del_x*del_y);
    %% Tao vector ve phai F
    F=zeros((N-1)*(N-1),1);
    for j=1:N-1
        for i = 1:N-1
            F((j-1)*(N-1)+i)=f(x(i+1),y(j+1));
        end
    end
   
    %% Tim nghiem xap xi u
    u=A\F;
    
    
     %% Tao ma tran cho ket qua chinh xac
    u_ex=zeros(N+1,N+1);
    for j=1:N+1
        for i=1:N+1
            u_ex(j,i)=u_exact(x(i),y(j));
        end
    end
    
    %% Tao ma tran chua ket qua xap xi va bien
    u_dis=zeros(N+1,N+1);
    for j=2:N
        for i=2:N
            u_dis(j,i)=u((j-2)*(N-1)+i-1);
        end
    end
    
    for i=2:N
       u_dis(1,i)= 0;
       u_dis(N+1,i)= 0;
       u_dis(1,N+1)= 0;
       u_dis(1,i)= 0;
    end
    u_dis(1,1) = 0;
    u_dis(N+1,N+1) = 0;
    
    %% Ve nghiem xap xi
    figure
    surf(x,y,u_dis)
    
    %% Vector c?a nghiem chinh xac
    vector_u_ex = zeros((N+1)^2,1);
    for j=1:N+1
        for i=1:N+1
            vector_u_ex((j-1)*(N+1)+i,1) = u_exact(x(i),y(j));
        end
    end
    
        %% Tao vector chua ket qua xap xi va bien
    vector_u_dis=zeros((N+1)^2,1);
    for i=1:N+1
        for j=1:N+1
            if(i==1)
                vector_u_dis((i-1)*(N+1)+j,1) =0;
            elseif(i==N+1)
                vector_u_dis((i-1)*(N+1)+j,1) =0;
            elseif(j==1)
                vector_u_dis((i-1)*(N+1)+j,1) =0;
            elseif(j==N+1)
                vector_u_dis((i-1)*(N+1)+j,1) =0;
            else
                vector_u_dis((i-1)*(N+1)+j,1)=u((i-2)*(N-1)+j-1);
            end
        end
    end
    
    %% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:(N+1)^2
        if (abs(vector_u_dis(i,1)-vector_u_ex(i,1)) > norm_max(inumber_mesh))
            norm_max(inumber_mesh)=abs(vector_u_dis(i,1)-vector_u_ex(i,1));
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:(N+1)^2
        norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(vector_u_dis(i)...
            -vector_u_ex(i))^2*del_x;
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:N^2
        if (abs(((vector_u_dis(i+1)-vector_u_ex(i+1))-(vector_u_dis(i)...
                -vector_u_ex(i)))/del_x) > norm_maxh1(inumber_mesh))
            norm_maxh1(inumber_mesh)=...
                abs(((vector_u_dis(i+1)-vector_u_ex(i+1))...
                -(vector_u_dis(i)-vector_u_ex(i)))/del_x);
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:N
        norm_h1(inumber_mesh)=...
            norm_h1(inumber_mesh)+(((vector_u_dis(i+1)-vector_u_ex(i+1))...
            -(vector_u_dis(i)-vector_u_ex(i)))/del_x)^2*del_x;
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions    
    figure
    plot(x,u_ex,'blue', x,u_dis,'red');
    xlabel('x');ylabel('value');
    title(['Comparison with N=', num2str(N)]);

    figure
    plot(y,u_ex,'blue', x,u_dis,'red');
    xlabel('y');ylabel('value');
    title(['Comparison with N=', num2str(N)]);

    %% Refine mesh
    N=2*N;
end

%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point),-log(norm_l2), 'red',...
    log(number_mesh_point), 1/2*log(number_mesh_point),'black',...
    log(number_mesh_point), log(number_mesh_point),'yellow');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_{l_2}','1/2x','x',...
    'Location','NorthEastOutside');  

