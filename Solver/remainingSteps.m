function unp1 = remainingSteps(ja,jb,sigma,unm1,un,unp1,oacc)

for i = ja:jb
    unp1(i) = 2*un(i)-unm1(i) + sigma^2 * (un(i+1) - 2*un(i) + un(i-1));
    if oacc == 4
       unp1(i) = unp1(i) - (sigma^2-sigma^4)/12*(un(i+2)-4*un(i+1)+6*un(i)-4*un(i-1)+un(i-2)); 
    end
end
    
end