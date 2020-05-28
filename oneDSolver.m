function [u,e] = oneDSolver_v2(def, tf, sigma, icase)
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

%% Setup
[x,dt,unm1,un,unp1] = setup(def,tf,sigma);

%% Initial Condition
unm1 = ICs(def,x,unm1);

%% First Time Step
un = firstStep(def,x,dt,sigma,unm1,un);
un = BC(def,1,dt,un,icase);

%% Remaining Time Steps
j = 2;
while j*dt <= tf
    % Update middle values
    unp1 = remainingSteps(def,sigma,unm1,un,unp1);
    
    % Update boundary conditions depending on case
    unp1 = BC(def,j,dt,unp1,icase);
    
    unm1 = un;
    un = unp1;
      
    j = j+1;
end

% Exact solution for case 1 
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

e = 0;
for i = 1:def.N+1
   e = max(e, abs(un(i) - y(x(i),3.15))); 
end

u = un;

end

