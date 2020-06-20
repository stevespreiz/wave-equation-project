% 4th order accurate method using exact soln for ghost points and first
% 2 time steps for now

def.a = 0;
def.b = 1;
def.N = 100;
def.c = 1;
tf = 3.15;
def.f = @(x) sin(3*pi/2*x);
def.g = @(x) 0;
def.l = @(t) 0;
def.r = @(t) 0;
sigma = .9;
y = @(x,t) 0.5*(def.f(x-def.c*t)+def.f(x+def.c*t));

% [u,e] = oneDSolver4(def,tf,sigma,1,4);

for n = [10 100 1000]
   def.N = n;
   [u,er] = oneDSolver4(def, tf, sigma, 1,4);
   e(log10(n)) = er;
end

h = [.1 .01 .001 ];
h2 = h.^2;
h3 = h.^3;
h4 = h.^4;
figure
loglog(h,e,'o',h,h,h, h2,h,h3,h,h4)
legend('error', 'h' ,'h^2','h^3', 'h^4')