close all
u = @(x,y,t) cos(2*pi*t)*sin(2*pi*x)*sin(2*pi*y);


N = 100;
X = linspace(-1,1,N+1);
Y = linspace(-1,1,N+1);
u0 = zeros(N+1,N+1);

for t = linspace(0,2,1000)
    for x = 1:N+1
        for y = 1:N+1
            u0(x,y) = u(X(x),Y(y),t);
        end
    end

    surf(X,Y,u0);
    xlabel('x')
    ylabel('y')
    shading interp
    zlim([-1 1])
    pause(.005)
end
