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


ys = @(x,y,t) cos(2*pi*t)*sin(2*pi*x)*sin(2*pi*y); 


% First Time Step
% un = firstStep(def,sigma,x,y,dt,ja,jb,unm1,un,nD,oacc);
for i = ja(1):jb(1)
    for j = ja(2):jb(2)
        un(i,j) = ys(x(i),y(j),dt);
    end
end

% e = 0;
% if nD == 2
%     for i = ja(1):jb(1)
%         for j = ja(2):jb(2)
%             e = max( e , abs( un(i,j) - ys( x(i) , y(j) , dt ) ) );
%         end
%     end
% end


% un = BC(def,sigma,x,y,1,dt,ja,jb,un,nD,icase,oacc);
for i = ja(1):jb(1)
    un(i,ja(2) - 1) = ys(x(i),y(ja(2)-1),dt);
    un(i,jb(2) + 1) = ys(x(i),y(jb(2)+1),dt);
end
for j = ja(2):jb(2)
    un(ja(1)-1, j) = ys(x(ja(1)-1),y(i),dt);
    un(jb(1)+1, j) = ys(x(jb(1)+1),y(i),dt);
end
% figure(1)
% xlabel('x')
% ylabel('y')
% Remaining Time Steps
n = 2;
while n*dt <= tf
    % Update middle values
    unp1 = timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc);
%     for i = ja(1):jb(1)
%         for j = ja(2):jb(2)
%             unp1(i,j) = ys(x(i),y(j),n*dt);
%         end
%     end
    
    % Apply Boundary Condition
%     unp1 = BC(def,sigma,x,y,n,dt,ja,jb,unp1,nD,icase,oacc);
    for i = ja(1):jb(1)
        unp1(i,ja(2) - 1) = ys(x(i),y(ja(2)-1),n*dt);
        unp1(i,jb(2) + 1) = ys(x(i),y(jb(2)+1),n*dt);
    end
    for j = ja(2):jb(2)
        unp1(ja(1)-1, j) = ys(x(ja(1)-1),y(i),n*dt);
        unp1(jb(1)+1, j) = ys(x(jb(1)+1),y(i),n*dt);
    end

    % Optional Animation
%     surf(x,y,unp1);
%     zlim([-1 1]);
%     shading interp
%     pause(.05);
%     xlabel('x')
%     ylabel('y')

    
    % Update arrays
    unm1 = un;
    un = unp1;
   
    n = n+1;
end

% 2D error check with exact sol. on desk pad
e = 0;
if nD == 2
    for i = ja(1):jb(1)
        for j = ja(2):jb(2)
            e = max( e , abs( un(i,j) - ys( x(i) , y(j) , tf ) ) );
        end
    end
end

% e = 0;
% if icase == 1
%     % Error calc with exact solution for case 1
%     y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t)+ integral(def.g,x-def.c*t,x+def.c*t)/def.c);
% 
%     for i = ja:jb
%         e = max(e,abs(un(i)-y(x(i),tf)));
%     end
% end

u = un;

end

