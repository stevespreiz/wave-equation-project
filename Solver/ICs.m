function [unm1] = ICs(def, x, unm1)
[~,m] = size(x);
for i = 1:m
   unm1(i) = def.f(x(i)); 
end

end