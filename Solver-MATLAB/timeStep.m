function [unp1] = timeStep(sigma,ja,jb,unm1,un,unp1,nD,oacc)
% Main time step over middle values

if nD == 1
    for i = ja:jb
        unp1(i) =  2*un(i) - unm1(i) + sigma^2*(un(i+1)-2*un(i)+un(i-1));
        if oacc > 2
            unp1(i) = unp1(i) - (sigma^2-sigma^4)/12*(un(i+2)-4*un(i+1)+6*un(i)-4*un(i-1)+un(i-2));
            if oacc > 4
                
            end
        end
    end
elseif nD == 2
    for i = ja(1):jb(1)
        for j = ja(2):jb(2)
            unp1(i,j) = 2*un(i,j) - unm1(i,j) + ...
                sigma(1)^2*(un(i-1,j)-2*un(i,j)+un(i+1,j)) + ...
                sigma(2)^2*(un(i,j-1)-2*un(i,j)+un(i,j+1));
            
%             if oacc > 2
%                unp1(i,j) = unp1(i,j) - (sigma(1)^2-sigma(1)^4)/12*...
%                    (un(i+2,j)-4*un(i+1,j)+6*un(i,j)-4*un(i-1,j)+un(i-2,j))+...
%                    (sigma(2)^2-sigma(2)^4)/12*...
%                    (un(i,j+2)-4*un(i,j+1)+6*un(i,j)-4*un(i,j-1)+un(i,j-2));
%             end
        end
    end
    
end

end

