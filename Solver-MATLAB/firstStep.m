function [un] = firstStep(def,sigma,x,dt,ja,jb,unm1,un,nD,oacc)
% Executes first time step

if nD == 1
    if oacc == 2
        for i = ja:jb
            un(i) = (1-sigma^2)*unm1(i)+dt*def.g(x(i))+sigma^2/2*(unm1(i-1)+unm1(i+1));
        end
    elseif oacc == 4
        for i = ja:jb
            un(i) = unm1(i) + ...
                    dt*def.g(x(i)) + ...
                    sigma^2/2*( (unm1(i+1)-2*unm1(i)+unm1(i-1)) - (unm1(i+2)-4*unm1(i+1)+6*unm1(i)-4*unm1(i-1)+unm1(i-2))/12  ) + ...
                    dt*sigma^2/6*( (def.g(x(i+1))-2*def.g(x(i))+def.g(x(i-1))) -(def.g(x(i+2))-4*def.g(x(i+1))+6*def.g(x(i))-4*def.g(x(i-1))+def.g(x(i-2)))/12 ) + ...
                    sigma^4/24*(unm1(i+2) - 4*unm1(i+1) + 6*unm1(i) - 4*unm1(i-1) + unm1(i-2));
        end
    end


end
