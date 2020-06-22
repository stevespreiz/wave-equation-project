function [unm1] = IC(def,x,y,unm1,nD)
% Sets initial position of wave

[~,m] = size(x);
[~,n] = size(y);

if nD == 1
    for i = 1:m
        unm1(i) = def.f(x(i));
    end
elseif nD == 2
    for i = 1:m
       for j = 1:n
          unm1(i,j) = def.f(x(i),y(j)); 
       end
    end
end

end

