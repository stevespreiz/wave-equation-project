function [u,e] = nDWaveSolver(def,sigma,tf,nD,icase,oacc)
% n-Dimensional Wave Equation Solver
% Inputs:
%   definition - struct [a,b,c,N,f,g,B.C.s], each is vector for each
%   dimenison
%   sigma      - CFL number vector for x,y,z,...
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

% Setup
[sigma,x,y,dt,unm1,un,unp1,ja,jb] = setup(def,sigma,tf,nD,oacc);
% Initial Condition
unm1 = IC(def,x,y,unm1,nD);


% First Time Step
un = firstStep(def,sigma,x,y,dt,ja,jb,unm1,un,nD,oacc);
un = BC(def,sigma,x,y,1,dt,ja,jb,un,nD,icase,oacc);

% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc);
    
    % Apply Boundary Condition
    unp1 = BC(def,sigma,x,y,n,dt,ja,jb,unp1,nD,icase,oacc);
    
    % Optional Animation
%     plot(x,unp1);
%     ylim([-1.1 1.1])
%     xlim([-.1 1.1])
%     pause;

    
    % Update arrays
    unm1 = un;
    un = unp1;
    
    
    n = n+1;
end

e = 0;
if icase == 1
    % Error calc with exact solution for case 1
    y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t)+ integral(def.g,x-def.c*t,x+def.c*t)/def.c);

    for i = ja:jb
        e = max(e,abs(un(i)-y(x(i),tf)));
    end
end

u = un;

end

