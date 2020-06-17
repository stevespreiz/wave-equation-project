function [sigma,x,dt,unm1,un,unp1,ja,jb] = setup(def,sigma,tf,nD,oacc)
% Initialization of variables and preallocation of major memory blocks

% Steps in space
dx = (def.b-def.a)/def.N;
x  = linspace(def.a - dx*oacc/2, def.b + dx*oacc/2, def.N+1+oacc);

% Steps in time
dttilde = sigma*dx/def.c;
nt      = ceil(tf/dttilde);
dt      = tf/nt;

% Update sigma value
sigma = def.c*dt/dx;

% Initialize spacial memory for 3 time steps
if nD == 1
    unm1 = zeros(def.N+1+oacc,1);
    un   = zeros(def.N+1+oacc,1);
    unp1 = zeros(def.N+1+oacc,1);
end

% Set indexing variables
ja = 1 + oacc/2;
jb = def.N + 1 + oacc/2;


end

