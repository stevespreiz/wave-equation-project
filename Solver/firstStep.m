function [un] = firstStep(ja,jb,def, x, dt, sigma, unm1, un)

for i = ja:jb
    un(i) = (1-sigma^2)*unm1(i)+dt*def.g(x(i))+sigma^2/2*(unm1(i-1)+unm1(i+1));
end

end