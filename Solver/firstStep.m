function [un] = firstStep(ja,jb,def, x, dt, sigma, unm1, un,oacc)

if oacc == 2
    for i = ja:jb
        un(i) = (1-sigma^2)*unm1(i)+dt*def.g(x(i))+sigma^2/2*(unm1(i-1)+unm1(i+1));
    end
else
    for i = ja:jb
        un(i) = unm1(i) +  dt*def.g(x(i)) + sigma^2/2*( (unm1(i+1)-2*unm1(i)+unm1(i-1))/2 - (unm1(i+2)-4*unm1(i+1)+6*unm1(i)-4*unm1(i-1)+unm1(i-2))/12  );
        un(i) = un(i) + dt*sigma^2/6*( (def.g(x(i+1))-2*def.g(x(i))+def.g(x(i-1)))/2 -(def.g(x(i+2))-4*def.g(x(i+1))+6*def.g(x(i))-4*def.g(x(i-1))+def.g(x(i-2)))/12 );
    end
end

end