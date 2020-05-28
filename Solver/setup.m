function [x,dt,unm1,un,unp1] = setup(def,tf,sigma)

x = linspace(def.a,def.b,def.N+1);
dx = x(2)-x(1);

dttilde = sigma*dx/def.c;
nt = ceil(tf/dttilde);
dt = tf/nt;

unm1 = zeros(def.N+1,1);
un   = zeros(def.N+1,1);
unp1 = zeros(def.N+1,1);

end
