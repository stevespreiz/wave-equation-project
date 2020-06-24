function [un] = firstStep(def,sigma,x,dt,ja,jb,unm1,un,nD,oacc,iTZ)
% Executes first time step

if iTZ == 0
    if nD == 1
        for i = ja:jb
            un(i) = unm1(i) + ...
                dt*def.g(x(i)) + ...
                sigma^2/2*( (unm1(i+1)-2*unm1(i)+unm1(i-1))) ;
            if oacc > 2
                un(i) = un(i)+...
                    -sigma^2/2*(unm1(i+2)-4*unm1(i+1)+6*unm1(i)-4*unm1(i-1)+unm1(i-2))/12  +...
                    dt*sigma^2/6*( (def.g(x(i+1))-2*def.g(x(i))+def.g(x(i-1))) -(def.g(x(i+2))-4*def.g(x(i+1))+6*def.g(x(i))-4*def.g(x(i-1))+def.g(x(i-2)))/12 ) + ...
                    sigma^4/24*(unm1(i+2) - 4*unm1(i+1) + 6*unm1(i) - 4*unm1(i-1) + unm1(i-2));
                if oacc > 4
                    un(i) = un(i)+...
                        sigma^2/2*(unm1(i+3)-6*unm1(i+2)+15*unm1(i+1)-20*unm1(i)+15*unm1(i-1)-6*unm1(i-2)+unm1(i-3))/90+...
                        dt*sigma^2/6*(def.g(x(i+3))-6*def.g(x(i+2))+15*def.g(x(i+1))-20*def.g(x(i))+15*def.g(x(i-1))-6*def.g(x(i-2))+def.g(x(i-3)))/90+...
                        -sigma^4/24*(unm1(i+3)-6*unm1(i+2)+15*unm1(i+1)-20*unm1(i)+15*unm1(i-1)-6*unm1(i-2)+unm1(i-3))/72+...
                        -dt*sigma^4/120*(def.g(x(i+3))-6*def.g(x(i+2))+15*def.g(x(i+1))-20*def.g(x(i))+15*def.g(x(i-1))-6*def.g(x(i-2))+def.g(x(i-3)))+...
                        sigma^6/720*(unm1(i+3)-6*unm1(i+2)+15*unm1(i+1)-20*unm1(i)+15*unm1(i-1)-6*unm1(i-2)+unm1(i-3));
                end
            end
        end
    end
    
elseif iTZ == 1
    if nD == 1
       for i = ja:jb
          un(i) = x(i)^2*dt^2; 
       end
    end
end

end
