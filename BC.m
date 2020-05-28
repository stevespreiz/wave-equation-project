function [unp1] = BC(def,j,dt,unp1,icase)

if icase == 1
    unp1(1) = def.l(j*dt);
    unp1(def.N+1) = (4*unp1(def.N)-unp1(def.N-1))/3;
else
    unp1(1) = def.l(j*dt);
    unp1(def.N+1) = def.r(j*dt);
end





end

