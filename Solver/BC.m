function [unp1] = BC(ja,jb,def,sigma,n,x,dt,unp1,icase,oacc)
% Currently only second order accurate

if oacc == 2
    if icase == 1
        % Assuming l_tt = 0
        unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
        unp1(jb+1) = unp1(jb-1);
    else
        unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
        unp1(jb+1) = 2*unp1(jb)-unp1(jb-1);
    end
else
    if icase == 1
        dx = x(2) - x(1);
        A = [2/3 -1/12; -2 1];
        b = A\[2/3*unp1(jb-1)-1/12*unp1(jb-2) + dx*def.r(n*dt)  ;  
            -2*unp1(jb-1)+unp1(jb-2) + 2*dx/sigma^2*(def.r(dt*(n+1))-2*def.r(dt*n)+def.r(dt*(n-1)))];
        
        unp1(jb+1) = b(1);
        unp1(jb+2) = b(2);
    else
        
    end
end

end

