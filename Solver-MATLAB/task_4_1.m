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
tf = 5;
def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0*x;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
sigma = .9;
for n = [10 20 50 100]
    def.N = n;
    [u,er] = nDWaveSolver(def, sigma, 1, 1,1,2);
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



%% 2-D Time
% close all

def.a = [-1 -1];
def.b = [1 1];
def.N = [100 100];
def.c = sqrt(1/2);
tf = 1;
def.f = @(x,y) sin(2*pi*x)*sin(2*pi*y);
def.g = @(x,y) 0*x*y;
def.l = @(t) 0*t;
def.r = @(t) 0*t;
def.t = @(t) 0*t;
def.bot = @(t) 0*t;
sigma = [.3 .8];
nD = 2;
icase = 2;
oacc = 2;

% Steps in space
dx = (def.b(1)-def.a(1))/def.N(1);
x  = linspace(def.a(1) - dx*oacc/2, def.b(1) + dx*oacc/2, def.N(1)+1+oacc);

dy = (def.b(2)-def.a(2))/def.N(2);
y  = linspace(def.a(2) - dy*oacc/2, def.b(2) + dy*oacc/2, def.N(2)+1+oacc);


% [un,e] = nDWaveSolver(def,sigma,tf,nD,icase,oacc);

for n = [10 100 200]
    def.N(1) = n;
    def.N(2) = n;
    [u,er] = nDWaveSolver(def, sigma, tf, nD,icase,oacc);
    if n == 10
        e = er;
    else
        e = [e er]; 
    end
end

h = [.1 .01 .005 ];
h2 = h.^2*100;
h4 = h.^4;
h6 = h.^6;
figure(2)
loglog(h,e,'o',h, h2,h,h4)
legend('error','h^2', 'h^4')
