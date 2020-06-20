function [unm1] = IC(def,x,unm1)
% Sets initial position of wave

[~,m] = size(x);

for i = 1:m
    unm1(i) = def.f(x(i));
end

end

