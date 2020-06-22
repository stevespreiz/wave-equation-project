% 4th order accurate method using exact soln for ghost points and first
% 2 time steps for now
%%
def.a = 0;
def.b = 1;
def.N = 10;
def.c = 1;
tf = 3;
def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0*x;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
sigma = .9;
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

nD = 1;
icase = 1;
oacc = 2;

[un,~] = nDWaveSolver(def,sigma,tf,nD,icase,oacc);
%%

% [u,e] = oneDSolver4(def,tf,sigma,1,4);
def.a = 0;
def.b = 1;
def.N = 10;
def.c = 1;
tf = 3;
def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0*x;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
sigma = .9;
for n = [10 100 1000]
   def.N = n;
   [u,er] = nDWaveSolver(def, sigma, 1, 1,1,2);
   e(log10(n)) = er;
end

h = [.1 .01 .001 ];
h2 = h.^2;
h3 = h.^3;
h4 = h.^4;
figure
loglog(h,e,'o',h,h,h, h2,h,h3,h,h4)
legend('error', 'h' ,'h^2','h^3', 'h^4')
%% 2-D Time
def.a = [0 0];
def.b = [pi pi];
def.N = [100 100];
def.c = 4;
tf = 1;
def.f = @(x,y) sin(x)*cos(y-pi/2);
def.g = @(x,y) 0*x*y;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
def.t = @(t) 0*t;
def.bot = @(t) 0*t;
sigma = [.4 .4];

nD = 2;
icase = 2;
oacc = 2;

[un,~] = nDWaveSolver(def,sigma,tf,nD,icase,oacc);