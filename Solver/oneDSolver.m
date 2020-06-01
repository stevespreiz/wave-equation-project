function [u,e] = oneDSolver(def, tf, sigma, icase, oacc)
% 1-D Wave Equation FDM solver
% Inputs:
%   definition - struct [a,b,c,N,f,g,B.C.s]
%   sigma - CFL number
%   icase - integer defining case to solve
%     Case 1:
%       B.C. u(a,t) = l(t), u_x(b,t) = 0
%       I.C. u(x,0) = f(x), u_t(x,0) = g(x)
%     Case 2:
%       B.C. u(a,t) = l(t), u(b,t) = r(t)
%       I.C. u(x,0) = f(x), u_t(x,0) = g(x)
%   oacc - order of accuracy (2 or 4)

%% Setup
[x,dt,unm1,un,unp1,ja,jb] = setup(def,tf,sigma,oacc);

%% Initial Condition
unm1 = ICs(def,x,unm1);

%% First Time Step
un = firstStep(ja,jb,def,x,dt,sigma,unm1,un);
un = BC(ja,jb,def,1,dt,un,icase,oacc);

%% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = remainingSteps(ja,jb,sigma,unm1,un,unp1,oacc);
    
    % Update boundary conditions depending on case
    unp1 = BC(ja,jb,def,n,dt,unp1,icase,oacc);
    
    unm1 = un;
    un = unp1;
      
    n = n+1;
end

% Exact solution for case 1 
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

e = 0;
for i = ja:jb
   e = max(e, abs(un(i) - y(x(i),3.15))); 
end


u = un(ja:jb);

end

