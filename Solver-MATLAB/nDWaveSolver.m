function [u,e] = nDWaveSolver(def,sigma,tf,nD,icase,oacc)
% n-Dimensional Wave Equation Solver
% Inputs:
%   definition - struct [a,b,c,N,f,g,B.C.s]
%   sigma      - CFL number
%   tf         - final time
%   nD         - number of dimensions
%   icase      - integer defining case to solve
%     Case 1:
%       B.C. u(a,t) = l(t), u_x(b,t) = 0
%       I.C. u(x,0) = f(x), u_t(x,0) = g(x)
%     Case 2:
%       B.C. u(a,t) = l(t), u(b,t) = r(t)
%       I.C. u(x,0) = f(x), u_t(x,0) = g(x)
%   oacc       - order of accuracy (2 or 4)
ys = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t)+ integral(def.g,x-def.c*t,x+def.c*t)/def.c);

% Setup
[sigma,x,dt,unm1,un,unp1,ja,jb] = setup(def,sigma,tf,nD,oacc);

% Initial Condition
unm1 = IC(def,x,unm1);

% First Time Step
un = firstStep(def,sigma,x,dt,ja,jb,unm1,un,nD,oacc);
% un = BC(def,sigma,x,1,dt,ja,jb,un,nD,icase,oacc);
un(ja-3) = ys(x(ja-3),dt);
un(ja-2) = ys(x(ja-2),dt);
un(ja-1) = ys(x(ja-1),dt);
un(jb+1) = ys(x(jb+1),dt);
un(jb+2) = ys(x(jb+2),dt);
un(jb+3) = ys(x(jb+3),dt);



% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc);
    
    % Apply Boundary Condition
%     unp1 = BC(def,sigma,x,n,dt,ja,jb,unp1,nD,icase,oacc);
    unp1(ja-3) = ys(x(ja-3),n*dt);
    unp1(ja-2) = ys(x(ja-2),n*dt);
    unp1(ja-1) = ys(x(ja-1),n*dt);
    unp1(jb+1) = ys(x(jb+1),n*dt);
    unp1(jb+2) = ys(x(jb+2),n*dt);
    unp1(jb+3) = ys(x(jb+3),n*dt);
    
    % Optional Animation
%     plot(x,unp1);
%     ylim([-1.1 1.1])
%     xlim([-.1 1.1])
%     pause(.1);
    
    % Update arrays
    unm1 = un;
    un = unp1;
    
    n = n+1;
end

e = 0;
if icase == 1
    % Error calc with exact solution for case 1
    ys = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t)+ integral(def.g,x-def.c*t,x+def.c*t)/def.c);

    for i = ja:jb
        e = max(e,abs(un(i)-ys(x(i),tf)));
    end
end

u = un(ja:jb);

end

