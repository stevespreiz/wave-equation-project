function unp1 = remainingSteps(def,sigma,unm1,un,unp1)

for i = 2:def.N
    unp1(i) = 2*un(i)-unm1(i) + sigma^2 * (un(i+1) - 2*un(i) + un(i-1));
end
    
end