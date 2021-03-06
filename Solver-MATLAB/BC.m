function [unp1] = BC(def,sigma,x,n,dt,ja,jb,unp1,nD,icase,oacc)
% Applies boundary conditions dependent on the case
dx = x(2)-x(1);
if nD == 1
    %Dirchlet left, Neumann right
    if icase == 1
        if oacc == 2
            % currently for l_tt = 0
            unp1(jb+1) = unp1(jb-1)+2*dx*def.r(n*dt);
            unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
        elseif oacc == 4
            % right hand Neumann 
            f = @(u) [u(1)-unp1(jb-1) - 1/3*(u(2)-2*u(1)+2*unp1(jb-1)-unp1(jb-2)) - def.r(n*dt)*2*dx ;
                      u(2) - 2*u(1) + 2*unp1(jb-1)-unp1(jb-2)];
            f0 = f([0;0]);
            f1 = f([1;0]);
            f2 = f([0;1]);
            A = [f1-f0,f2-f0];
            b = -1.*f0;
            u = A\b;
            
            unp1(jb+1) = u(1);
            unp1(jb+2) = u(2);
                  
                  
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
        elseif oacc == 6
            % Left side Dirchlet assuming l_tt = 0
            f = @(u) [1/90*unp1(ja+3)-3/20*unp1(ja+2)+3/2*unp1(ja+1)-49/18*unp1(ja)+3/2*u(1)-3/20*u(2)+1/90*u(3) ;
                -1/6*unp1(ja+3)+2*unp1(ja+2)-13/2*unp1(ja+1)+28/3*unp1(ja)-13/2*u(1)+2*u(2)-1/6*u(3);
                unp1(ja+3)-6*unp1(ja+2)+15*unp1(ja+1)-20*unp1(ja)+15*u(1)-6*u(2)+u(3)];
            f0 = f([0;0;0]);
            f1 = f([1;0;0]);
            f2 = f([0;1;0]);
            f3 = f([0;0;1]);
            A = [f1-f0,f2-f0,f3-f0];
            b = -1.*f0;
            u = A\b;
            unp1(ja-1) = u(1);
            unp1(ja-2) = u(2);
            unp1(ja-3) = u(3);
            
            % Right side Neumannn
            % assuming right hand condition is O(t) (only used in first row
            % of f
            q = 30; % pretty sure this is right
            f = @(u) [u(1)-unp1(jb-1) - 1/6*(u(2)-2*u(1)+2*unp1(jb-1)-unp1(jb-2)) + 1/q*(u(3)-4*u(2)+5*u(1)-5*unp1(jb-1)+4*unp1(jb-2)-unp1(jb-3))-def.r(n*dt)*dx*2;
                u(2)-2*u(1)+2*unp1(jb-1)-unp1(jb-2) + 1/q*(u(3)-4*u(2)+5*u(1)-5*unp1(jb-1)+4*unp1(jb-2)-unp1(jb-3));
                u(3)-4*u(2)+5*u(1)-5*unp1(jb-1)+4*unp1(jb-2)-unp1(jb-3)];
            f0 = f([0;0;0]);
            f1 = f([1;0;0]);
            f2 = f([0;1;0]);
            f3 = f([0;0;1]);
            A = [f1-f0,f2-f0,f3-f0];
            b = -1.*f0;
            u = A\b;
            unp1(jb+1) = u(1);
            unp1(jb+2) = u(2);
            unp1(jb+3) = u(3);
        end
        % Dirchlet both sides
    elseif icase == 2
        if oacc == 2
            unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
            unp1(jb+1) = 2*unp1(jb)-unp1(jb-1);
        elseif oacc == 4
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
            
            % Right hand Dirchlet with r_tt = 0
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
        elseif oacc == 6
            % Left side Dirchlet assuming l_tt = 0
            f = @(u) [1/90*unp1(ja+3)-3/20*unp1(ja+2)+3/2*unp1(ja+1)-49/18*unp1(ja)+3/2*u(1)-3/20*u(2)+1/90*u(3) ;
                -1/6*unp1(ja+3)+2*unp1(ja+2)-13/2*unp1(ja+1)+28/3*unp1(ja)-13/2*u(1)+2*u(2)-1/6*u(3);
                unp1(ja+3)-6*unp1(ja+2)+15*unp1(ja+1)-20*unp1(ja)+15*u(1)-6*u(2)+u(3)];
            f0 = f([0;0;0]);
            f1 = f([1;0;0]);
            f2 = f([0;1;0]);
            f3 = f([0;0;1]);
            A = [f1-f0,f2-f0,f3-f0];
            b = -1.*f0;
            u = A\b;
            unp1(ja-1) = u(1);
            unp1(ja-2) = u(2);
            unp1(ja-3) = u(3);
            
            % Right side Dirchlet assuming r_tt = 0
            f = @(u) [1/90*unp1(jb-3)-3/20*unp1(jb-2)+3/2*unp1(jb-1)-49/18*unp1(ja)+3/2*u(1)-3/20*u(2)+1/90*u(3) ;
                -1/6*unp1(jb-3)+2*unp1(jb-2)-13/2*unp1(jb-1)+28/3*unp1(ja)-13/2*u(1)+2*u(2)-1/6*u(3);
                unp1(jb-3)-6*unp1(jb-2)+15*unp1(jb-1)-20*unp1(ja)+15*u(1)-6*u(2)+u(3)];
            f0 = f([0;0;0]);
            f1 = f([1;0;0]);
            f2 = f([0;1;0]);
            f3 = f([0;0;1]);
            A = [f1-f0,f2-f0,f3-f0];
            b = -1.*f0;
            u = A\b;
            unp1(jb+1) = u(1);
            unp1(jb+2) = u(2);
            unp1(jb+3) = u(3);
        end
    end
end

end