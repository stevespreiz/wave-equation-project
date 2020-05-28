function [unm1] = ICs(def, x, unm1)

for i = 1:def.N+1
   unm1(i) = def.f(x(i)); 
end

end