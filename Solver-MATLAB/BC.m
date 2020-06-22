function [unp1] = BC(def,sigma,x,y,n,dt,ja,jb,unp1,nD,icase,oacc)
% Applies boundary conditions dependent on the case

if nD == 1
    %Dirchlet left, Neumann right
    if icase == 1
        if oacc == 2
            % currently for u_x = 0 and l_tt = 0
            unp1(jb+1) = unp1(jb-1);
            unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
        elseif oacc == 4
            dx = x(2) - x(1);
            % right hand Neumann - need to change to DDFA
            A = [2/3 -1/12; -2 1];
            b = A\[2/3*unp1(jb-1)-1/12*unp1(jb-2) + dx*def.r(n*dt)  ;
                -2*unp1(jb-1)+unp1(jb-2) + 2*dx/sigma^2*(def.r(dt*(n+1))-2*def.r(dt*n)+def.r(dt*(n-1)))];
            
            unp1(jb+1) = b(1);
            unp1(jb+2) = b(2);
            
            % left hand Dirchlet assuming l_tt = 0 - discrete delta
            % function approach
            f = @(u) [def.c^2/dx^2*(u(1)-2*unp1(ja)+unp1(ja+1) - (u(2)-4*u(1)+6*unp1(ja)-4*unp1(ja+1)+unp1(ja+2))/12);
                def.c^4/dx^4*(u(2)-4*u(1)+6*unp1(ja)-4*unp1(ja+1)+unp1(ja+2))];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0,f2-f0];
            b = -1.*f0;
            u = A\b;
            
            unp1(ja-1) = u(1);
            unp1(ja-2) = u(2);
        end
        % Dirchlet both sides
    elseif icase == 2
        if oacc == 2
            unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
            unp1(jb+1) = 2*unp1(jb)-unp1(jb-1);
        elseif oacc == 4
            dx = x(2) - x(1);
            % left hand Dirchlet assuming l_tt = 0 - discrete delta
            % function approach
            f = @(u) [def.c^2/dx^2*(u(1)-2*unp1(ja)+unp1(ja+1) - (u(2)-4*u(1)+6*unp1(ja)-4*unp1(ja+1)+unp1(ja+2))/12);
                def.c^4/dx^4*(u(2)-4*u(1)+6*unp1(ja)-4*unp1(ja+1)+unp1(ja+2))];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0,f2-f0];
            b = -1.*f0;
            ut = A\b;
            
            unp1(ja-1) = ut(1);
            unp1(ja-2) = ut(2);
            
            % Right hand Richlet with r_tt = 0
            f = @(u) [def.c^2/dx^2*(unp1(jb-1)-2*unp1(jb)+u(1) - (unp1(jb-2)-4*unp1(jb-1)+6*unp1(jb)-4*u(1)+u(2))/12);
                def.c^4/dx^4*(unp1(jb-2)-4*unp1(jb-1)+6*unp1(jb)-4*u(1)+u(2))];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0,f2-f0];
            b = -1.*f0;
            ut = A\b;
            
            unp1(jb+1) = ut(1);
            unp1(jb+2) = ut(2);
        end
    end
elseif nD == 2
    if icase == 1
        
    % All sides pinned at 0    
    elseif icase == 2
        if oacc == 2
            for i = ja(1):jb(1)
               unp1(i,ja(2) - 1) = 2*unp1(i,ja(2)) - unp1(i,ja(2)+1);
               unp1(i,jb(2) + 1) = 2*unp1(i,jb(2)) - unp1(i,jb(2)-1);
            end
            for j = ja(2):jb(2)
                unp1(ja(1)-1, j) = 2*unp1(ja(1), j) - unp1(ja(1)+1,j);
                unp1(jb(1)+1, j) = 2*unp1(jb(1), j) - unp1(jb(1)-1,j);
            end
        end
    end
end


end