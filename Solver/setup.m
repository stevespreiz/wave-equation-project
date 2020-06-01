function [x,dt,unm1,un,unp1,ja,jb] = setup(def,tf,sigma, oacc)

dx = (def.b-def.a)/def.N;
x = linspace(def.a-dx*oacc/2,def.b+dx*oacc/2,def.N+1+oacc);

dttilde = sigma*dx/def.c;
nt = ceil(tf/dttilde);
dt = tf/nt;

unm1 = zeros(def.N+1+oacc,1);
un   = zeros(def.N+1+oacc,1);
unp1 = zeros(def.N+1+oacc,1);

ja = 1 + oacc/2;
jb = def.N + 1 + oacc/2;

end