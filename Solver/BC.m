function [unp1] = BC(ja,jb,def,n,dt,unp1,icase,oacc)
% Currently only second order accurate

if icase == 1
    % Assuming l_tt = 0
    unp1(ja-1) = 2*unp1(ja)-unp1(ja+1);
    unp1(jb+1) = unp1(jb-1);
else
    unp1(ja-1) = def.l(n*dt);
    unp1(jb+1) = def.r(n*dt);
end

end

