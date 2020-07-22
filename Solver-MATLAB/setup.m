function [sigma,x,y,dt,unm1,un,unp1,ja,jb] = setup(def,sigma,tf,nD,oacc)
% Initialization of variables and preallocation of major memory blocks

if nD == 1
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
    
    unm1 = zeros(def.N+1+oacc,1);
    un   = zeros(def.N+1+oacc,1);
    unp1 = zeros(def.N+1+oacc,1);
    
    
    % Set indexing variables
    ja = 1 + oacc/2;
    jb = def.N + 1 + oacc/2;
    
    y = 0;
elseif nD == 2
    % Steps in space
    dx = (def.b(1)-def.a(1))/def.N(1);
    x  = linspace(def.a(1) - dx*oacc/2, def.b(1) + dx*oacc/2, def.N(1)+1+oacc);
    
    dy = (def.b(2)-def.a(2))/def.N(2);
    y  = linspace(def.a(2) - dx*oacc/2, def.b(2) + dx*oacc/2, def.N(2)+1+oacc);
    
    
    % Steps in time
    dttilde = sigma(1)*dx*dy/def.c*sqrt(1/dx^2+1/dy^2);
    nt      = ceil(tf/dttilde);
    dt      = tf/nt;
    
    % Update sigma value
    sigma(1) = def.c*dt/dx;
    sigma(2) = def.c*dt/dy;
    
    % Initialize spacial memory for 3 time steps
    
    unm1 = zeros(def.N(1)+1+oacc,def.N(2)+1+oacc);
    un   = zeros(def.N(1)+1+oacc,def.N(2)+1+oacc);
    unp1 = zeros(def.N(1)+1+oacc,def.N(2)+1+oacc);
    
    
    % Set indexing variables
    ja(1) = 1 + oacc/2;
    ja(2) = 1 + oacc/2;
    jb(1) = def.N(1) + 1 + oacc/2;
    jb(2) = def.N(2) + 1 + oacc/2;
end


end

