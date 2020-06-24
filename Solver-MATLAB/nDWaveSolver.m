function [u,e] = nDWaveSolver(def,sigma,tf,nD,icase,oacc,iTZ)
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
%   iTZ        - flag for Twilight Zone
%     Case 1: u_e = x^2*t^2
%   oacc       - order of accuracy (2 or 4)

% Setup
[sigma,x,dt,unm1,un,unp1,ja,jb] = setup(def,sigma,tf,nD,oacc);

% Initial Condition
unm1 = IC(def,x,unm1,iTZ);

% First Time Step
un = firstStep(def,sigma,x,dt,ja,jb,unm1,un,nD,oacc,iTZ);
un = BC(def,sigma,x,1,dt,ja,jb,un,nD,icase,oacc,iTZ);

% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = timeStep(def,sigma,x,n,dt,ja,jb,unm1,un,unp1,nD,oacc,iTZ);
    
    % Apply Boundary Condition
    unp1 = BC(def,sigma,x,n,dt,ja,jb,unp1,nD,icase,oacc,iTZ);
    
    % Optional Animation
    plot(x,unp1);
    ylim([-1.1 1.1])
    xlim([-.1 1.1])
    pause(.1);
    
    % Update arrays
    unm1 = un;
    un = unp1;
    
    n = n+1;
end

e = 0;
if iTZ == 0
    if icase == 1
        % Error calc with exact solution for case 1
        ys = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t)+ integral(def.g,x-def.c*t,x+def.c*t)/def.c);

        for i = ja:jb
            e = max(e,abs(un(i)-ys(x(i),tf)));
        end
    end
elseif iTZ == 1
    if icase == 1
        % Error calc with exact solution for case 1
        ys = @(x,t) x^2*t^2;

        for i = ja:jb
            e = max(e,abs(un(i)-ys(x(i),tf)));
        end
    end
end

u = un(ja:jb);

end

