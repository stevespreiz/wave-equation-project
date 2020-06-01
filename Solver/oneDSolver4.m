function [u,e] = oneDSolver4(def, tf, sigma, icase, oacc)
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


% Exact solution for case 1 
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

%% Setup
[x,dt,unm1,un,unp1,ja,jb] = setup(def,tf,sigma,oacc);

%% Initial Condition
unm1 = ICs(def,x,unm1);
% plot(x(ja:jb),unm1(ja:jb));
% pause

%% First Time Step
% un = firstStep(ja,jb,def,x,dt,sigma,unm1,un);
% un = BC(ja,jb,def,1,dt,un,icase,oacc);
for i = 1:def.N+3
    un(i) = y(x(i),dt);
end

% plot(x(ja:jb),un(ja:jb));
% pause
%% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = remainingSteps(ja,jb,sigma,unm1,un,unp1,oacc);
    
    % Update boundary conditions depending on case
%     unp1 = BC(ja,jb,def,n,dt,unp1,icase,oacc);
    unp1(1) = y(x(1),n*dt);
    unp1(2) = y(x(2),n*dt);
    unp1(def.N+2) = y(x(def.N+2),n*dt);
    unp1(def.N+3) = y(x(def.N+3),n*dt);
    
    unm1 = un;
    un = unp1;
    
%     plot(x(ja:jb),un(ja:jb));
%     ylim([-1 1])
%     pause(.01)
    
    n = n+1;
end


e = 0;
for i = ja:jb
   e = max(e, abs(un(i) - y(x(i),3.15))); 
end


u = un(ja:jb);

end

