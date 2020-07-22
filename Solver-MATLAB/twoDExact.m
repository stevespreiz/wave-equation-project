close all
%%
ys = @(x,y,t) cos(2*pi*t)*sin(2*pi*x)*sin(2*pi*y);


N = 100;
X = linspace(-1,1,N+1);
Y = linspace(-1,1,N+1);
u0 = zeros(N+1,N+1);

for t = 1
    for x = 1:N+1
        for y = 1:N+1
            u0(x,y) = ys(X(x),Y(y),t);
        end
    end
    figure(10)
    surf(X,Y,u0);
    xlabel('x')
    ylabel('y')
    shading interp
    zlim([-1 1])
    pause(.005)
end
