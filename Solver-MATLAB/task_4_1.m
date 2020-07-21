% 4th order accurate method using exact soln for ghost points and first
% 2 time steps for now
%%
def.a = 0;
def.b = 1;
def.N = 100;
def.c = 1;
tf = 2;
def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0*x;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
sigma = .5;

nD = 1;
icase = 1;
oacc = 6;

[un,~] = nDWaveSolver(def,sigma,tf,nD,icase,oacc);
%%
close all
% [u,e] = oneDSolver4(def,tf,sigma,1,4);

def.a = 0;
def.b = 1;
def.c = 1;

def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0*x;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
sigma = .2;
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

nD = 1;
icase = 1;
oacc = 4;
tf = 2;

for n = [10 20 50 100]
    def.N = n;
    [u,er] = nDWaveSolver(def, sigma, tf,nD,icase,oacc);
    if n == 10
        e = er;
    else
        e = [e er];
    end
end

h = [.1 .05 .02 .01 ];
h2 = h.^2;
h4 = h.^4;
h6 = h.^6;
figure
loglog(h,e,'o',h, h2,h,h4,h,h6)
legend('error','h^2', 'h^4', 'h^6')
