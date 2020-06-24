function [unm1] = IC(def,x,unm1,iTZ)
% Sets initial position of wave

[~,m] = size(x);
if iTZ == 0
    for i = 1:m
        unm1(i) = def.f(x(i));
    end
elseif iTZ == 1
   for i = 1:m
       unm1(i) = 0;
   end
end

end

